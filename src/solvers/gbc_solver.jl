
using JuMP
using BenchmarkTools
using CSV
import MathOptInterface as MOI


"""
    solve_with_GBC!(inst::Instance, param::GBCparam)

Solves the passed instance using generalized Benders cuts generation. Note that the master will be adjusted during the solution process. 
If you want to use this instance of master problem again, make a copy before. 
Also, master should **not** contain any variables name _subObj_ and both master and subproblems should be minimization problems. 
You should set the _optimizer_ for each JuMP model yourself beforehand.

# Arguments
- 'inst::Instance': The problem instance including the master and the subproblems.
- 'param::GBCparam': Parameters
"""
function solve_with_GBC!(inst::Instance, param::GBCparam)
    master = inst.master
    subs = inst.subproblems

    # register solver
    @debug "Starting setup of the GBC solver."
    new_stat!(param.stats, "Solver", "GBCSolver")
    new_stat!(param.stats, "time_limit", param.runtime)
    new_stat!(param.stats, "NSub", length(subs))
    new_stat!(param.stats, "NFeasCuts", 0)
    new_stat!(param.stats, "NOptCuts", 0)
    new_stat!(param.stats, "SepaTime", 0)  # time spend in separator
    new_stat!(param.stats, "SepaTimeCut", 0)  # time spend in separator for generating cuts only

    # do some initail checks for master and sub solvers
    @debug "Doing some checks if master and sub were created correctly."
    check(master, param)
    for sub in subs
        check(sub, param)
    end

    # add subObj variables, one for each sub_problem, to the master, and add them to objective
    @debug "Start with initialization of GBC solver."
    names = [name(sub) for sub in subs]
    subObj = @variable(master.model, subObj[names])
    @objective(master.model, Min, objective_function(master.model) + sum(subObj))

    # add partial decomposition to master if requested
    if !isnothing(master.partial_decomposition)
        master.partial_decomposition(master.model, subObj)
    end

    # debbug output
    try
        write_to_file(
            master.model,
            param.output_folder_path * "/master.$(param.file_format_output)",
        )
        set_optimizer_attribute(master.model, "LogFile", param.output_folder_path*"/gbc_mip_log.txt")
    catch err 
        @error "Could not write model to file or set log file for folder $(param.output_folder_path). Error is $err"
    end

    # build the sub LPs for Benders subroutine 
    clps, runtime_init = init_connectorLPs(subs, master.link_vars, subObj, param) 
    new_stat!(param.stats, "runtime_preprocessingGBC", runtime_init)

    # set time limit and number of threads
    true_runtime = param.runtime - runtime_init
    set_time_limit_sec(master.model, true_runtime)
    set_attribute(master.model, MOI.NumberOfThreads(), param.threads_master)
    for sub in subs
        set_nthreads(sub, param.threads_sub_con)
    end

    # add callback to master and solve 
    msol_cuts_mapping = Dict()  # a mapping of master solution to found lazy constraints
    if true_runtime > 0
        @debug "Finished model construction. Now proceeding to optimization process with GBC. Remaining runtime is $true_runtime"
        set_attribute(master.model, MOI.LazyConstraintCallback(), cb -> gbc_callback_function(cb, master, names, clps, subObj, msol_cuts_mapping, param))
    else
        @debug "We do not add any callbacks to GBCSolver as we run into a time out during the preprocessing."
        new_stat!(param.stats, "GBCStatus", "Timeout_Submodel")
    end

    try 
        optimize!(master.model)
    catch e
        if (e isa TimeoutException)
            @warn "We run into a timeout when solving a Submodel with GBCSolver. Note that this implies that the Subproblem run for the time passed by timelimit and did not terminate. "
            new_stat!(param.stats, "GBCStatus", "Timeout_Submodel")
        else
            @error "GBCsolver suffered an error: $e"
            @error stacktrace(catch_backtrace())
            #showerror(stdout, e, catch_backtrace())
            new_stat!(param.stats, "GBCStatus", "Terminate")
            #rethrow(e)
        end
    else
        # print solution and collected data
        print_collected_cuts(param, msol_cuts_mapping)
        if termination_status(master.model) == MOI.OPTIMAL || termination_status(master.model) == MOI.LOCALLY_SOLVED || termination_status(master.model) == MOI.TIME_LIMIT
            if primal_status(master.model) == MOI.FEASIBLE_POINT
                mobj = objective_value(master.model)
                xsol = Dict(a => value(master.link_vars[a]) for a in master.A)
                @debug "The master objective is $(mobj) and solution is $(xsol)."
                print_solution_to_file(mobj, xsol, param)
                new_stat!(param.stats, "Opt", mobj)

                # print full MIP solution to file
                solution = Dict(JuMP.name(x) => JuMP.value(x) for x in all_variables(master.model))
                write(param.output_folder_path*"/full_master_solution.json", JSON.json(solution))
            end

            # set status 
            status = termination_status(master.model)
            set_optimization_status_stats(status, param)
        else 
            @debug "The master MIP is infeasible with termination status: $(termination_status(master.model))"
            set_optimization_status_stats(termination_status(master.model), param)
        end
        
        # save run data to statistics
        new_stat!(param.stats, "runtime", solve_time(master.model))
        new_stat!(param.stats, "gap", JuMP.relative_gap(master.model))
        new_stat!(param.stats, "BNodes", MOI.get(master.model, MOI.NodeCount()))

        # TODO: find and print correct second level solutions?
    end

end


"""
    init_connectorLPs(subs, link_vars, subObjs, param)

Init all the ConnectorLP instances. Stop if time limit is reached.

# Return 
    - The build list of ConnectorLP objects. If we run into time out, the list is not complete
    - The overall time needed to init all ConnectorLP. In case of time out return 'param.runtime'.
"""
function init_connectorLPs(subs, link_vars, subObjs, param::GBCparam)
    @debug "Beginn building LPs for each subproblems required for Benders steps."
    timelimit_inner = param.runtime
    connectors = []
    try 
        for s in subs
            time_s = @elapsed begin
                con = build_connectorLP(s, link_vars, subObjs[name(s)], param, timelimit_inner)
                push!(connectors, con)
            end
            timelimit_inner = timelimit_inner - time_s

            if timelimit_inner <= 0
                return connectors, param.runtime # all time used up (but generaly the functions before throw TimeoutException)
            end
        end
        return connectors, param.runtime - timelimit_inner # return remaining time after preprocessing
    catch err
        if (err isa TimeoutException)
            @debug "Caught timeout exception while generating ConnectorLPs for GBCSolver"
        else
            rethrow(err)
        end
    end
    return connectors, param.runtime
end

"""
    build_connectorLP(sub::SubSolver, link_vars_master::Dict, subObj, parameter::GBCparam)

Generate the ConnectorLP objects that form the Benders subproblems within our hierarchical decomposition.

# Arguments
- 'sub': The subsolver that forms the original sub_problem behind the solver. 
- 'link_vars_master': The linking variables (in the master MIP) that are later used for the generelized Benders cuts.
- 'subObj': The master MIP objective variable representing the contribution of the sub_problem. We set a bound on it based on solving the original sub_problem with master objective function
- 'parameter': The solver parameters.
- 'timelimit': The time limit for this subroutine 

# Returns
- 'lp::ConnectorLP': The LP
"""
function build_connectorLP(sub::SubSolver, link_vars_master::Dict, subObjvar, parameter::GBCparam, timelimit)
    # Generate the ConnectorLP for the passed sub_problem
    @debug "Starting building of ConnectorLP for subpoblem $(name(sub))."

    # Construct LP (no objective or constraints)
    # Note that we use upper bounds to prevent unbounded solutions, see the ConnectorLP implimentation
    myLP = Model(() -> get_next_optimizer(parameter.solver))
    @variable(myLP, s <= parameter.infinity_num)
    @variable(myLP, k[sub.A] >= 0)
    @variable(myLP, 0 <= g <= parameter.infinity_num)

    # set number of threads
    set_attribute(myLP, MOI.NumberOfThreads(), parameter.threads_sub_con)

    # disable output of LP
    set_silent(myLP)

    # calculate a bound on the big M coefficients
    lbm = compute_lower_bound_master_contribution(sub, parameter, timelimit)
    @debug "For subpoblem $(name(sub)), the found lower bound to the master objective contribution is $(lbm). With this, construction of ConnectorLP finished."

    # set computed lower bound for sub_problem objective variables
    set_lower_bound(subObjvar, lbm)

    # build ConnectorLP obj
    return ConnectorLP(myLP, sub.A, link_vars_master, sub, lbm)
end



function gbc_callback_function(cb_data, master::Master, sub_names, clps, subObj, msol_cuts_mapping::Dict, parameter::GBCparam)
    # x are the linking variables and clps the connectors (one for each sub)
    # subObj are the obj. vars. in master (for each sub)

    status = callback_node_status(cb_data, master.model)
    if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
        return
    elseif status == MOI.CALLBACK_NODE_STATUS_INTEGER
        # I would just like to mention that @belapsed does not work (as I cannot pass global vars into the scope...)
        sepatime = @elapsed begin
            # `callback_value(cb_data, x)` is integer (to some tolerance).
            subObj_val = Dict(
                name => callback_value(cb_data, subObj[name]) for name in sub_names
            )
            @debug "Current values of the sub objectives are $(subObj_val)"
            x_vals = Dict(
                a => callback_value(cb_data, master.link_vars[a]) for a in master.A
            )
            round_master_solution(x_vals) # round to integer
            @debug "Current values of the master linking variables are $(x_vals)"
            lazy = []

            # Because of multi thread (and some start trouble I cannot explain) we resolve the subproblems for same master solution multiple times. 
            # To avoid this, we save the found cuts for each master solution and fall back on them before resolving the sub_problem
            msolkey = key_master_sol(x_vals, master.A)
            #@debug "The current master_key=$(msolkey) and the saved solutions are $(msol_cuts_mapping).")
            if haskey(msol_cuts_mapping, msolkey)
                # if we already solved this sub_problem, we just recover the found lazy constraints 
                lazy = msol_cuts_mapping[msolkey]
                @debug "We recovered the lazy cuts for solution $(x_vals) without resolving subproblems."
            else
                # solve each sub and add cut to list of cuts "lazy"
                @debug "We found no saved lazy cuts for current solution and solve subproblems now."
                for con in clps
                    cuttime = @elapsed begin
                        # solve corresponding connectorLP
                        ## there does not seem to be an easy way to obtain the remaining runtime in solver-independent callbacks. So we stick with a worst case approximation for the time limit.
                        feas = false
                        cut=nothing
                        pobj=-1
                        feas, cut, pobj = genBenders_cut!(con, x_vals, parameter, parameter.runtime)
                        # in case of timeout within one subsolver, one of the Submodels run for timelimit time. We let this result in a timeout error 

                        # generate cut from found solution
                        if feas
                            cutfeas = @build_constraint(cut >= 1)
                            @debug "Adding feasibility cut $(cutfeas) to the master problem for sub $(name(con.sub_solver))."
                            add_stat!(parameter.stats, "NFeasCuts", 1)
                            push!(lazy, cutfeas)
                        else
                            # add opt cut only if required, i.e., it truly cuts of solution 
                            if subObj_val[name(con)] + 1e-6 < pobj
                                cutopt = @build_constraint(cut <= subObj[name(con)])
                                @debug "Adding optimality cut $(cutopt) to the master problem for sub $(name(con.sub_solver))."
                                add_stat!(parameter.stats, "NOptCuts", 1)
                                push!(lazy, cutopt)
                            end
                        end
                    end
                    add_stat!(parameter.stats, "SepaTimeCut", cuttime)
                end

                # save that we found lazy const. for this sol
                msol_cuts_mapping[msolkey] = lazy
            end

            # add lazy cuts to master model
            map(cu -> MOI.submit(master.model, MOI.LazyConstraint(cb_data), cu), lazy)

            # output cuts to file in case of debbug mode
            if should_debbug_print(parameter)
                outfilecut = parameter.output_folder_path * "/mastercuts_last_update.txt"
                append_constraintlist_to_file(lazy, outfilecut)
            end
        end
        add_stat!(parameter.stats, "SepaTime", sepatime)
    else
        @assert status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
        return
    end
end


