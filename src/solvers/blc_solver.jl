# Basic implementation of a Benders-like cuts solver. The big M must be provided by the user.

function solve_with_BLC!(inst::Instance, param::BLCparam)
    blcm::BlCMaster = inst.master

    # register solver
    @debug "Starting setup of the BlC solver."
    new_stat!(param.stats, "Solver", "BlCSolver")
    new_stat!(param.stats, "time_limit", param.runtime)
    new_stat!(param.stats, "NSub", length(inst.subproblems))
    new_stat!(param.stats, "Blcuts", 0)
    new_stat!(param.stats, "SepaTime", 0)  # time spend in separator

    # do some initial checks for master and sub solvers
    @debug "Doing some checks if master and sub were created correctly for Benders-like cuts solver."
    check(blcm, param)
    for sub in inst.subproblems
        check(sub, param)
    end

    # save HPR to lp file and set log file
    try
        write_to_file(blcm.hpr, param.output_folder_path * "/hpr.$(param.file_format_output)")
        set_optimizer_attribute(blcm.hpr, "LogFile", param.output_folder_path *"/blc_mip_log.txt")
    catch err 
        @error "Could not write model to file or set log file for folder $(param.output_folder_path). Error is $err"
    end

    # set runtime and number of threads
    set_time_limit_sec(blcm.hpr, param.runtime)
    set_attribute(blcm.hpr, MOI.NumberOfThreads(), param.threads_master)
    for sub in inst.subproblems
        set_nthreads(sub, param.threads_sub_con)
    end

    # add callback to master and solve 
    @debug "Finished model construction. Now proceeding to optimization process with GBC."
    msol_cuts_mapping = Dict()  # a mapping of master solution to found lazy constraints
    mtimeout = Dict()  # the master solutions for which we hit a timeout
    set_attribute(
        blcm.hpr,
        MOI.LazyConstraintCallback(),
        cb -> gbc_callback_function_blc(cb, inst, msol_cuts_mapping, param),
    )
    try 
        optimize!(blcm.hpr)
    catch e
        if (e isa TimeoutException)
            @warn "We run into a timeout when solving a Submodel with BlCSolver. Note that this implies that the Subproblem run for the time passed by timelimit and did not terminate. This indicates a bug or a very challenging second-level problem."
            new_stat!(param.stats, "BlCStatus", "Timeout_Submodel")
        else
            @error "BlCsolver suffered an error: $e"
            @error stacktrace(catch_backtrace())
            new_stat!(param.stats, "BlCStatus", "Terminate")
        end
    else
        # print solution to logger and collected data
        print_collected_cuts(param, msol_cuts_mapping)
        status = termination_status(blcm.hpr)
        if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED || status == MOI.TIME_LIMIT
            if primal_status(blcm.hpr) == MOI.FEASIBLE_POINT
                mobj = objective_value(blcm.hpr)
                xsol = Dict(a => value(blcm.link_vars[a]) for a in blcm.A)
                @debug "The master objective is $(mobj) and solution is $(xsol)."
                print_solution_to_file(mobj, xsol, param)
                new_stat!(param.stats, "Opt", mobj)

                # print MIP solution to file
                solution = Dict(JuMP.name(x) => JuMP.value(x) for x in all_variables(blcm.hpr))
                write(param.output_folder_path*"/solution.json", JSON.json(solution))
            end

            # set status 
            set_optimization_status_stats(status, param)

            # save run data to statistics
            new_stat!(param.stats, "runtime", solve_time(blcm.hpr))
            new_stat!(param.stats, "gap", JuMP.relative_gap(blcm.hpr))
            new_stat!(param.stats, "BNodes", MOI.get(blcm.hpr, MOI.NodeCount()))

            # TODO: find and print correct second level solutions?
        else
            @debug "The HPR MIP is infeasible with status: $status"
            set_optimization_status_stats(status, param)
        end
    end

end

function gbc_callback_function_blc(cb_data, inst::Instance, msol_cuts_mapping::Dict, parameter::BLCparam)
    blcm::BlCMaster = inst.master
    status = callback_node_status(cb_data, blcm.hpr)

    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        return
    elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # I would just like to mention that @belapsed does not work (as I cannot pass global vars into the scope...)
        sepatime = @elapsed begin

            # `callback_value(cb_data, x)` is integer (to some tolerance).
            x_vals = Dict(a => callback_value(cb_data, blcm.link_vars[a]) for a in blcm.A)
            round_master_solution(x_vals)  # round to integer
            @debug "Starting Benders-like cuts: Current values of the master linking variables are $(x_vals)."
            lazy = []

            # Because of multi thread (and some start trouble I cannot explain) we resolve the subproblems for same master solution multiple times. 
            ## To avoid this, we save the found cuts for each master solution and fall back on them before resolving the sub_problem
            msolkey = key_master_sol(x_vals, blcm.A)
            ## @debug "The current master_key=$(msolkey) and the saved solutions are $(msol_cuts_mapping).")
            if haskey(msol_cuts_mapping, msolkey)
                # if we already solved this sub_problem, we just recover the found lazy constraints 
                lazy = msol_cuts_mapping[msolkey]
                @debug "Finished Benders-like Cuts separation: We recovered the lazy cuts for solution $(x_vals) without resolving subproblems."
            else
                # solve each sub and add cut to list of cuts "lazy"
                @debug "We found no saved lazy cuts for current solution and solve subproblems now."
                for sub in inst.subproblems
                    # retrive current objective function value of sub_problem
                    subobjcurrent = callback_value(cb_data, blcm.sub_objectives[name(sub)])

                    # solve sub_problem 
                    ## there does not seem to be an easy way to obtain the remaining runtime in solver-independent callbacks. So we stick with a worst case approximation for the time limit.
                    feasible = false
                    subopt = -1
                    y_vals = Dict()
                    feasible, subopt, y_vals = solve_sub_for_x(sub, x_vals, parameter, parameter.runtime)
                    # in case of timeout within one subsolver, one of the Submodels run for timelimit time. We let this result in a timeout error 
                    
                    if !feasible
                        error("Terminate BlC solver: The passed first-level solution was not feasible for subsolver $(name(sub)). x=$x_vals")
                    end

                    # check if Benders-like cut is necessary and add it if this is the case
                    ## Attention: this leads to strange bugs with multi-threads settings. Just add Benders-like cuts whenever you find one
                    bigMterms = 0
                    for a in blcm.A
                        bigMterms +=
                            blcm.big_m(a, name(sub)) *
                            y_vals[a] *
                            (1 - blcm.link_vars[a])
                    end
                    cutopt = @build_constraint(
                        blcm.sub_objectives[name(sub)] <= subopt + bigMterms
                    )
                    @debug "Adding Benders-like cut $(cutopt) to the master problem for sub $(name(sub))."
                    add_stat!(parameter.stats, "Blcuts", 1)
                    push!(lazy, cutopt)
                end

                # save that we found lazy const. for this sol
                msol_cuts_mapping[msolkey] = lazy
            end

            # add lazy cuts to master model
            map(cu -> MOI.submit(blcm.hpr, MOI.LazyConstraint(cb_data), cu), lazy)

            # output cuts to file in case of debbug mode
            if should_debbug_print(parameter)
                outfilecut = parameter.output_folder_path * "/Blcuts.txt"
                append_constraintlist_to_file(lazy, outfilecut)
            end
        end
        add_stat!(parameter.stats, "SepaTime", sepatime)
    else
        @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        return
    end
end
