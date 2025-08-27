using JuMP

"""
The LP Benders solver for a single sub_problem. 
"""
struct ConnectorLP{T}
    lp::Model  # The lp model. Variables should be registered as :s, :k[A], :g
    A::Vector{T}  # Iterable of resources
    link_vars::Dict{T,VariableRef}  # The linking variables (from master) as dict. Keys=A
    sub_solver::SubSolver  # The sub_solver that is used to solve the sub_problem
    lower_bound_obj_contribution::Number  # The lower bound on the contribution to the first level objective function
end

function name(clp::ConnectorLP)
    return name(clp.sub_solver)
end



"""
    genBenders_cut!(subLP::ConnectorLP, link_vals::Dict{Tuple{Int64, Int64}, Float64}, params::GBCparam, time_limit)

Generate an general Benders (feasibility or optimality) cut.

# Arguments
- 'subLP::ConnectorLP{T}': The LP that must be solved for the sub_problem. Note that the LP is modified during the solving process.
- 'link_vals::Dict{T, Float64}': The mapping of resource to values of the linking variables from the master.
- 'params::GBCparam': Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns

- 'feas::Bool': True iff we need a feasibility cut.
- 'cut': The left-hand (i.e., non-trivial) side of the cut.
- 'pobj': The contribution of the cut to the objective for the current solution (0 if feasibility cut).
"""
function genBenders_cut!(subLP::ConnectorLP{T}, link_vals::Dict{T,Float64}, params::GBCparam, time_limit) where T
    # adjust sub_problem by setting new objective
    @debug "We now solve for the found optimal master solution the ConnectorLP $(name(subLP.sub_solver))."

    foundfeas, optL2, y_vals = solve_sub_for_x(subLP.sub_solver, link_vals, params, time_limit)
    new_obj = subLP.lp[:s] - optL2 * subLP.lp[:g]
    new_obj +=
        -sum([
            round(capacity_linking(subLP.sub_solver, a, params) * link_vals[a]) * subLP.lp[:k][a]  # Todo: lets see if rounding helps numerics
            for a in subLP.A
        ])
    @objective(subLP.lp, Max, new_obj)

    # solve sub LP for new master
    @debug "Start iterative solution procedure for connector $(name(subLP.sub_solver))."

    # based if we found solution for current sub_problem, we need optimality or feasibility cut
    feas = !foundfeas  # a bit irreating, as we need feasibility cut iff sub is infeasible
    if feas
        # if infeasible, we do not have unbounded LP because of variable bounds added when initialized. This has some advantages:
        # 1. We can do the pareto-optimal cuts approach on this adjusted problem
        # 2. It seems to be numerical more stable for this LP that relying on dual information
        # 3. There seems to be only very little support for retrieving unbounded ray in JuMP

        # fix g to make sure that we obtain a ray with desired properties
        fix_g_constraint = @constraint(subLP.lp, subLP.lp[:g] == 0)

        @debug "Solving ConnectorLP $(name(subLP.sub_solver)) for feasibility cut generation."
        time_iterate = iterate_subsolver(subLP, params, time_limit)  # solve with new constraint

        
        # pareto optimality step for feasibility cuts
        if params.pareto == PARETO_OPTIMALITY_AND_FEASIBILITY
            time_limit_pareto = time_limit - time_iterate
            @debug "Start pareto-optimal Benders cut generation procedure for connector $(name(subLP.sub_solver)) for feasibility cut construction with remaining time limit $time_limit_pareto."
            pareto_optimal_decomposition(subLP, new_obj, optL2, params, time_limit_pareto)
        end

        # build feas cut
        cut = build_feas_cut(subLP)
        pobj = 0
    else
        # solve sub_problem 
        @debug "Solving ConnectorLP $(name(subLP.sub_solver)) for optimality cut generation."
        time_iterate = iterate_subsolver(subLP, params, time_limit)

        # pareto optimality step for optimality cuts
        if params.pareto == PARETO_OPTIMALITY_AND_FEASIBILITY || params.pareto == PARETO_OPTIMALITY_ONLY
            time_limit_pareto = time_limit - time_iterate
            @debug "Start pareto-optimal Benders cut generation procedure for connector $(name(subLP.sub_solver)) for optimality cut construction with remaining time limit $time_limit_pareto."
            pareto_optimal_decomposition(subLP, new_obj, optL2, params, time_limit_pareto)
        end

        #build opt cut
        cut = build_opt_cut(subLP, optL2, y_vals, params)
        pobj = value(new_obj)
    end

    # clean up and return
    if feas
        delete(subLP.lp, fix_g_constraint)
    end
    pareto_optimal_decomposition_cleanup(subLP)
    @debug "Finished solving connectLP $(name(subLP.sub_solver))."
    return feas, cut, pobj
end


function build_opt_cut(subLP::ConnectorLP, optL2, y_vals, params::GBCparam)
    master_vars = subLP.link_vars

    # s term of the cut
    sval = value(subLP.lp[:s])
    cut = sval

    # g term of the cut
    gval = value(subLP.lp[:g])
    bigMterm = sval - optL2 * gval - subLP.lower_bound_obj_contribution
    @debug "We use big_m=$(bigMterm) for this optimality cut"
    theta_xXg =
        optL2 * gval + sum(bigMterm * (1 - master_vars[a]) * y_vals[a] for a in subLP.A)  # The big M already contains the gval coefficient 
    cut -= theta_xXg


    # k term of the cut
    kvals = value.(subLP.lp[:k])
    cut -= sum(kvals[a] * master_vars[a] for a in subLP.A)
    return cut
end

function build_feas_cut(subLP::ConnectorLP)
    master_vars = subLP.link_vars
    sval = value(subLP.lp[:s])
    kvals = value.(subLP.lp[:k])
    fcut = sum(kvals[a] / sval * master_vars[a] for a in subLP.A)
    return fcut
end

"""
    check_solution_status_LP(lp::JuMP.Model)

Check if our connector LP model terminated with an optimal solution
"""
function check_solution_status_LP(me::ConnectorLP)
    status = termination_status(me.lp) 
    if status == MOI.TIME_LIMIT
        @debug "We hit the time limit while solving the ConnectorLP of user $(name(me))"
        throw(TimeoutException( "We hit the time limit while ConnectorLP of user $(name(me))" ))
    end


    if status == MOI.INFEASIBLE
        @debug "The ConnectorLP $(name(me)) is infeasible. Starting computations of IIS. But most of the time this is implied by numerical issues."
        compute_conflict!(me.lp)
        @debug "IIS computed, now printing it"
        iis_model, _ = copy_conflict(me.lp)
        if get_attribute(me.lp, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
            @debug iis_model # printing to file just causes error...
        else
            @debug "The IIS was not computed succesfully??? Most likely numerics??"
        end
        error("ConnectorLP $(name(me)) was infeasible. Computed IIS but stopping solution process (as it is clearly a bug). Most likely, it was caused by numerical issues.")
    elseif status == MOI.DUAL_INFEASIBLE
        error("The ConnectorLP $(name(me)) was unbounded what violates the way we handle its. Because we add significantly large bounds, we avoid unboundnes; so, this should not have happended.")
    end

    if !(status == MOI.OPTIMAL)
        error("ConnectorLP $(name(me)) could not be solver to optimality but terminated with status $(status). Connot continue as this is undefined behavior.")
    end
end


"""
    iterate_subsolver(subLP::ConnectorLP, params::GBCparam)

Solve the LP by a separation procedure. 
Uses an iterative while loop in the implementation, avoiding stack overflow exceptions compared to this function recursive version.

In case the separation takes longer that the set time limit, throw an exception. 

# Return
    The time it required to execute this function
"""
function iterate_subsolver(subLP::ConnectorLP, params::GBCparam, time_limit)
    # this part of the solver can run into quite some nasty infinity loops. To prevent the software just freezing (or seeming to freeze for the user),
    ## we stop the separation process in case we reach the timelimit set in the parameters (what still can take long, but at least the user is expected to wait this long) 
    current_time = time()

    violated_cut = true # true as long as violated constraint could exist in LP
    while violated_cut
        # solve the sub_problem iteratively (but first debbug output)
        if params.debbug_out
            write_to_file(
                subLP.lp,
                "$(params.output_folder_path)/conLP_$(name(subLP.sub_solver)).$(params.file_format_output)",
            )
        end

        # first, solve the LP
        #@debug subLP.lp # this debug output is not helpfull
        optimize!(subLP.lp) # Assumption: Solving the LP consumes neglectable time
        check_solution_status_LP(subLP)  

        # solve sub_problem for found solution
        kvals = Dict(a => value(subLP.lp[:k][a]) for a in subLP.A)
        @debug "The found sub_problem ConnectorLP solution is s=$(value(subLP.lp[:s])), g=$(value(subLP.lp[:g])), and k=$(kvals)"
        sub_solver = separation!(
            subLP.sub_solver,
            value(subLP.lp[:s]),
            value(subLP.lp[:g]),
            kvals,
            params,
            time_limit
        )

        # if we found a violated constraint, add it and resolve
        if sub_solver.vio
            # add constraints
            @debug "For connector $(name(subLP.sub_solver)), we found a sub_problem solution that violates current LP solution and uses resources $(sub_solver.A_sub)."
            new_const_left =
                subLP.lp[:s] - sum(subLP.lp[:k][a] for a in sub_solver.A_sub; init=0) -
                sub_solver.obj_second_level * subLP.lp[:g]
            c = @constraint(subLP.lp, new_const_left <= sub_solver.obj_first_level)
            @debug "We added an violated constraint $(c) for connector $(name(subLP.sub_solver)). Continue separation."

            # check for time limit
            if current_time + params.runtime < time() 
                @error "We reached the time limit when solving ConnectorLP $(name(subLP.sub_solver)). Terminating cut generation process."
                throw(TimeoutException("We reached the time limit when solving ConnectorLP $(name(subLP.sub_solver)). Terminating GBC solution procedure."))
            end

            # resolve
            violated_cut = true
        else
            # no violated constraint many more
            violated_cut = false
        end
    end

    return time() - current_time # return how long the process took
end


"""
    iterate_subsolver_recursive(subLP::ConnectorLP, params::GBCparam)

@Deprecated !

    Solves the current version of the LP using a recursive call of this function. 
Compared to the non-recursive 'iterate_subsolver' version of this function, it should be in my understanding of Julia.
However, it can easily lead to an 'StackOverflowError' if multiple constraints need to be added within one iteration. 
Therefore, this function DEPRECATED and should be avoided. 
"""
function iterate_subsolver_recursive(subLP::ConnectorLP, params::GBCparam)
    # solve the sub_problem iteratively (but first debbug output)
    if params.debbug_out
        write_to_file(
            subLP.lp,
            "$(params.output_folder_path)/conLP_$(name(subLP.sub_solver)).$(params.file_format_output)",
        )
    end

    # first, solve the LP
    #@debug subLP.lp # this debug output is not helpfull
    optimize!(subLP.lp)
    check_solution_status_LP(subLP)

    # solve sub_problem for found solution
    kvals = Dict(a => value(subLP.lp[:k][a]) for a in subLP.A)
    @debug "The found sub_problem ConnectorLP solution is s=$(value(subLP.lp[:s])), g=$(value(subLP.lp[:g])), and k=$(kvals)"
    sub_solver = separation!(
        subLP.sub_solver,
        value(subLP.lp[:s]),
        value(subLP.lp[:g]),
        kvals,
        params,
    )

    # if we found a violated constraint, add it and resolve
    if sub_solver.vio
        # add constraints
        @debug "For connector $(name(subLP.sub_solver)), we found a sub_problem solution that violates current LP solution and uses resources $(sub_solver.A_sub)."
        new_const_left =
            subLP.lp[:s] - sum(subLP.lp[:k][a] for a in sub_solver.A_sub; init=0) -
            sub_solver.obj_second_level * subLP.lp[:g]
        c = @constraint(subLP.lp, new_const_left <= sub_solver.obj_first_level)
        @debug "We added an violated constraint $(c) for connector $(name(subLP.sub_solver)). Continue separation."

        # resolve
        iterate_subsolver_recursive(subLP, params)
    end
end



"""
    pareto_optimal_decomposition(subLP::ConnectorLP, lp_obj, g_obj_coef, params::GBCparam)

When called after LP was solved to optimality, transforms the LP to generate a pareto-optimal cut and then resolves it.

# Arguments
- 'subLP::ConnectorLP': The LP
- 'lp_obj': The LP objective function (as JuMP expression)
- 'g_obj_coef': The coefficient of g in the objective function, i.e., the current solution value of the follower problem for selected master solution, or 0 if infeasible.
- 'params::GBCparam': The parameters
- 'time_limit': The runtime limit for this subroutine. In case of time out, throw an 'TimeoutException'
"""
function pareto_optimal_decomposition(subLP::ConnectorLP, lp_obj, g_obj_coef, params::GBCparam, time_limit)
    # resolve the subLP with steps necessary for pareto-optimal cuts
    cssol = value(subLP.lp[:s])
    cgsol = value(subLP.lp[:g])

    # The 0.1 term -> you still want to minimize the g
    bigMterm = (cssol - g_obj_coef * cgsol - subLP.lower_bound_obj_contribution) + 0.1 # this is just an educated guess 
    if bigMterm < 0
        error(
            "We got a negative big_m in ConnectorLP $(name(subLP)) because s=$(cssol), g=$(cgsol), gcoef=$(g_obj_coef), and bound=$(subLP.lower_bound_obj_contribution).",
        )
    end

    # build adjusted model
    current_obj = objective_value(subLP.lp)
    pBd_obj = sum(subLP.lp[:k]) + bigMterm * subLP.lp[:g]
    @objective(subLP.lp, Min, pBd_obj)
    @constraint(subLP.lp, pBd_const, lp_obj == current_obj)

    # print pareto adjusted model to file in debbug mode 
    if params.debbug_out
        write_to_file(
            subLP.lp,
            "$(params.output_folder_path)/conLP_ParetoInit_$(name(subLP.sub_solver)).$(params.file_format_output)",
        )
    end

    # resolve optimization problem
    iterate_subsolver(subLP, params, time_limit)
end


function pareto_optimal_decomposition_cleanup(subLP::ConnectorLP)
    # do the cleanup required for pareto-optimal Benders cuts
    if haskey(object_dictionary(subLP.lp), :pBd_const)
        delete(subLP.lp, subLP.lp[:pBd_const])
        unregister(subLP.lp, :pBd_const)
    end
    # no need to set new obj. This will happend when new master solution passed for next run
end





Base.show(io::IO, conLP::ConnectorLP) = print(io, "ConnectorLP_$(name(conLP.sub_solver))")





