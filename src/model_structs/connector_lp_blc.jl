## This is a version of the ConnectorLP that solves the Lagrangian dual that appears when approximating Benders-like cuts. 
##  Although we share many similarities with ConnectorLP, it essentially solves a different LP. Therefore, we decided to provide a new implementation.


### Constructor and basic functions ###
struct ConnectorLP_BlC{T}
    lp::Model  # The lp model. Variables should be registered as :s and :k[A]
    A::Vector{T}  # Iterable of resources
    link_vars::Dict{T,VariableRef}  # The linking variables (from master) as dict. Keys=A
    sub_solver::SubSolver  # The sub_solver that is used to solve the sub_problem
end

function (::Type{ConnectorLP_BlC{T}})(
    solver, infinity_num, A::Vector{T}, link_vars::Dict{T,VariableRef}, sub_solver::SubSolver
) where T
    # build LP
    myLP = Model(() -> get_next_optimizer(solver))
    @variable(myLP, s <= infinity_num)
    @variable(myLP, k[sub.A] >= 0)
    @variable(myLP, 0 <= g <= infinity_num)

    # build my struct
    ConnectorLP_BlC{T}(
      myLP, A, link_vars, sub_solver
    )
end

function name(clp::ConnectorLP_BlC)
    return name(clp.sub_solver)
end


### Cut Generation ###
"""
    genBenderslike_cut!(subLP::ConnectorLP, link_vals::Dict{Tuple{Int64, Int64}, Float64}, params::GBCparam, time_limit)

Generate a Benders-like optimality cut. 
We assume that this function will only be called for a first-level solution for which the second level is feasible. 
We solve the Lagrangian dual corresponding to the Benders-like cut and return the obtained cut coefficient. 

# Arguments
- 'subLP::ConnectorLP_BlC{T}': The LP(-wrapper) that must be solved for the sub_problem. Note that the LP is modified during the solving process.
- 'link_vals::Dict{T, Float64}': The mapping of resource to values of the linking variables from the master.
- 'params::GBCparam': Parameters passed down from the main solver. We respect the 'pareto' and 'warmstart' parameters.  
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns

- 'feas::Bool': True iff we need a feasibility cut.
- 'cut': The left-hand (i.e., non-trivial) side of the cut.
- 'pobj': The contribution of the cut to the objective for the current solution (0 if feasibility cut).
"""
function genBenderslike_cut!(subLP::ConnectorLP_BlC{T}, link_vals::Dict{T,Float64}, params::GBCparam, time_limit) where T
    # adjust sub_problem by setting new objective
    @debug "We now solve for the found optimal master solution the ConnectorLP_BlC $(name(subLP.sub_solver))."

    # Note: We assume that we have a feasible solution
    new_obj = subLP.lp[:s]
    new_obj +=
        sum([
            round(capacity_linking(subLP.sub_solver, a, params) * link_vals[a]) * subLP.lp[:k][a]  # TODO: lets see if rounding helps numerics
            for a in subLP.A
        ])
    @objective(subLP.lp, Min, new_obj)

    # solve sub LP for new master
    @debug "Start iterative solution procedure for ConnectorLP_BlC $(name(subLP.sub_solver))."
    time_iterate = iterate_subsolver_BlC(subLP, params, time_limit)

    # pareto optimality step for optimality cuts
    if params.pareto == PARETO_OPTIMALITY_AND_FEASIBILITY || params.pareto == PARETO_OPTIMALITY_ONLY
        time_limit_pareto = time_limit - time_iterate
        @debug "Start pareto-optimal Benders cut generation procedure for ConnectorLP_BlC $(name(subLP.sub_solver)) for optimality cut construction with remaining time limit $time_limit_pareto."
        pareto_optimal_decomposition_BlC(subLP, new_obj, optL2, params, time_limit_pareto)
    end

    #build opt cut
    cut = build_opt_cut_BlC(subLP, optL2, y_vals, params)
    pobj = value(new_obj)

    # clean up and return
    if feas
        delete(subLP.lp, fix_g_constraint)
    end
    pareto_optimal_decomposition_cleanup(subLP)
    if !params.warmstart
        # okay, I understand that these tests are kind of interesting, BUT I never felt so stupid implementing a feature
        @debug "As requested cleaning up ConnectorLP_BlC $(name(subLP.sub_solver)) by removing all generated constraints."
        for cref in all_constraints(subLP.lp; include_variable_in_set_constraints=false)
            @debug "Now deleting constraint $cref from ConnectorLP_BlC $(name(subLP.sub_solver))."
            delete(subLP.lp, cref)
        end
    end
    @debug "Finished solving ConnectorLP_BlC $(name(subLP.sub_solver))."
    return feas, cut, pobj
end

function build_opt_cut_BlC(subLP::ConnectorLP_BlC, optL2, y_vals, params::GBCparam)
    master_vars = subLP.link_vars

    # s term of the cut
    sval = value(subLP.lp[:s])
    cut = sval

    # k term of the cut
    kvals = value.(subLP.lp[:k])
    cut = sum(kvals[a] * master_vars[a] for a in subLP.A)
    return cut
end


"""
    iterate_subsolver(subLP::ConnectorLP_BlC, params::GBCparam)

Solve the LP by a separation procedure. 
Uses an iterative while loop in the implementation, avoiding stack overflow exceptions compared to this function recursive version.
In case the separation takes longer that the set time limit, throw an exception. 

# Return
    The time it required to execute this function
"""
function iterate_subsolver_BlC(subLP::ConnectorLP_BlC, params::GBCparam, time_limit)
    # this part of the solver can run into quite some nasty infinity loops. To prevent the software just freezing (or seeming to freeze for the user),
    ## we stop the separation process in case we reach the timelimit set in the parameters (what still can take long, but at least the user is expected to wait this long) 
    current_time = time()

    violated_cut = true # true as long as violated constraint could exist in LP
    while violated_cut
        # solve the sub_problem iteratively (but first debbug output)
        if params.debbug_out
            write_to_file(
                subLP.lp,
                "$(params.output_folder_path)/conLPBlC_$(name(subLP.sub_solver)).$(params.file_format_output)",
            )
        end

        # first, solve the LP
        #@debug subLP.lp # this debug output is not helpfull
        optimize!(subLP.lp) # Assumption: Solving the LP consumes neglectable time
        check_solution_status_LP(subLP)  

        # solve sub_problem for found solution
        kvals = Dict(a => value(subLP.lp[:k][a]) for a in subLP.A) 
        @debug "The found sub_problem ConnectorLP_BlC solution is s=$(value(subLP.lp[:s])) and k=$(kvals)"
        sub_solver = separation_BlC!(
            subLP.sub_solver,
            value(subLP.lp[:s]),
            kvals,
            params,
            time_limit
        )

        # if we found a violated constraint, add it and resolve
        if sub_solver.vio
            # add constraints
            @debug "For ConnectorLP_BlC $(name(subLP.sub_solver)), we found a sub_problem solution that violates current LP solution and uses resources $(sub_solver.A_sub)."
            new_const_left = subLP.lp[:s] + sum(subLP.lp[:k][a] for a in sub_solver.A_sub; init=0)
            c = @constraint(subLP.lp, new_const_left >= sub_solver.obj_second_level)
            @debug "We added an violated constraint $(c) for ConnectorLP_BlC $(name(subLP.sub_solver)). Continue separation."

            # check for time limit
            if current_time + params.runtime < time() 
                @error "We reached the time limit when solving ConnectorLP_BlC $(name(subLP.sub_solver)). Terminating cut generation process."
                throw(TimeoutException("We reached the time limit when solving ConnectorLP_BlC $(name(subLP.sub_solver)). Terminating GBC solution procedure."))
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
    pareto_optimal_decomposition_BlC(subLP::ConnectorLP_BlC, lp_obj, g_obj_coef, params::GBCparam)

When called after LP was solved to optimality, transforms the LP to generate a pareto-optimal cut and then resolves it.

# Arguments
- 'subLP::ConnectorLP_BlC': The LP
- 'lp_obj': The LP objective function (as JuMP expression)
- 'g_obj_coef': The coefficient of g in the objective function, i.e., the current solution value of the follower problem for selected master solution, or 0 if infeasible.
- 'params::GBCparam': The parameters
- 'time_limit': The runtime limit for this subroutine. In case of time out, throw an 'TimeoutException'
"""
function pareto_optimal_decomposition_BlC(subLP::ConnectorLP_BlC, lp_obj, g_obj_coef, params::GBCparam, time_limit)
    # resolve the subLP with steps necessary for pareto-optimal cuts
    cssol = value(subLP.lp[:s])

    # build adjusted model
    current_obj = objective_value(subLP.lp)
    pBd_obj = sum(subLP.lp[:k])
    @objective(subLP.lp, Min, pBd_obj)
    @constraint(subLP.lp, pBd_const, lp_obj == current_obj)

    # print pareto adjusted model to file in debbug mode 
    if params.debbug_out
        write_to_file(
            subLP.lp,
            "$(params.output_folder_path)/conLP_BlC_ParetoInit_$(name(subLP.sub_solver)).$(params.file_format_output)",
        )
    end

    # resolve optimization problem
    iterate_subsolver_BlC(subLP, params, time_limit)
end


function pareto_optimal_decomposition_cleanup(subLP::ConnectorLP_BlC)
    # do the cleanup required for pareto-optimal Benders cuts
    if haskey(object_dictionary(subLP.lp), :pBd_const)
        delete(subLP.lp, subLP.lp[:pBd_const])
        unregister(subLP.lp, :pBd_const)
    end
    # no need to set new obj. This will happend when new master solution passed for next run
end






### Auxiliary Functions ###
"""
    check_solution_status_LP(me::ConnectorLP_BlC)

Check if our connector LP model terminated with an optimal solution
"""
function check_solution_status_LP(me::ConnectorLP_BlC)
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
        error("The ConnectorLP_BlC $(name(me)) was unbounded what violates the way we handle its. Because we add significantly large bounds, we avoid unboundnes; so, this should not have happended.")
    end

    if !(status == MOI.OPTIMAL)
        error("ConnectorLP_BlC $(name(me)) could not be solver to optimality but terminated with status $(status). Connot continue as this is undefined behavior.")
    end
end

