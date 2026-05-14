using JuMP
import MathOptInterface as MOI

"""
A JuMP-based subsolver for bilevel follower problems that enforces bilevel
feasibility by adding persistent BlC-style no-good cuts over the copied linking
variables.

The working MIP is solved with the current separation objective. After each
optimal solve, an oracle function is used to solve the original second-level
problem for the found linking-variable pattern. If the candidate solution is not
optimal for the original second-level objective, a BlC-style cut is added to the
working MIP and the model is resolved from scratch. The generated cuts remain in
the working MIP across later subsolver calls.
"""
struct SubSolverBlCJuMP{T} <: SubSolver
    name::String
    mip_model::JuMP.Model
    A::Vector{T}
    link_varsC::Any
    y_vars::Any
    r_objterm::GenericAffExpr
    c_objterm::GenericAffExpr
    oracle_solve::Function
    big_m::Function
    cached_oracle_solutions::Dict{Any,NamedTuple}
end

function _copy_affexpr(expr::GenericAffExpr, ref_map)
    copied = zero(expr)
    JuMP.add_to_expression!(copied, JuMP.constant(expr))
    for (coef, var) in JuMP.linear_terms(expr)
        JuMP.add_to_expression!(copied, coef, ref_map[var])
    end
    return copied
end

function _build_default_oracle(mip_model, A, link_varsC, y_vars, c_objterm)
    oracle_model, oracle_ref_map = copy_model(mip_model)
    oracle_link_vars = Dict(a => oracle_ref_map[link_varsC[a]] for a in A)
    oracle_y_vars = Dict(a => oracle_ref_map[y_vars[a]] for a in A)
    oracle_c_objterm = _copy_affexpr(c_objterm, oracle_ref_map)
    set_silent(oracle_model)
    @objective(oracle_model, Min, oracle_c_objterm)

    return function (x_vals, params, time_limit)
        set_optimizer(oracle_model, () -> get_next_optimizer(params.solver))
        set_silent(oracle_model)
        set_time_limit_sec(oracle_model, time_limit)
        set_seed!(oracle_model, params.solver, get_seed(params))

        if haskey(params.stats.data, "threads_sub_con_used")
            oracle_threads =
                (
                    (params isa GBCparam && params.parallel_separation) ||
                    (params isa BlCLagparam && params.parallel_separation)
                ) ? 1 :
                used_nthreads(params.stats, "threads_sub_con")
            set_attribute(oracle_model, MOI.NumberOfThreads(), oracle_threads)
        end

        @constraint(oracle_model, fixc[a=A], oracle_link_vars[a] == round(x_vals[a]))
        try
            optimize!(oracle_model)

            status = termination_status(oracle_model)
            if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
                return false, 0, 0, Dict(a => 0.0 for a in A)
            elseif status == MOI.TIME_LIMIT
                throw(TimeoutException("We hit the time limit while solving the default oracle model for the bilevel subsolver."))
            elseif status != MOI.OPTIMAL
                error("The default oracle lower-level solve terminated with status $status")
            end

            y_vals = Dict(a => value(oracle_y_vars[a]) for a in A)
            return true, value(oracle_c_objterm), 0, y_vals
        finally
            for a in A
                delete(oracle_model, fixc[a])
            end
            unregister(oracle_model, :fixc)
        end
    end
end

function SubSolverBlCJuMP(
    name,
    mip_model,
    A::Vector{T},
    y_vars,
    r_objterm,
    c_objterm,
    oracle_solve::Function,
    big_m::Function,
) where T
    link_varsC = _create_link_varsC(mip_model, name, A, y_vars)
    return SubSolverBlCJuMP{T}(
        name,
        mip_model,
        A,
        link_varsC,
        y_vars,
        r_objterm,
        c_objterm,
        oracle_solve,
        big_m,
        Dict{Any,NamedTuple}(),
    )
end

function SubSolverBlCJuMP(
    name,
    mip_model,
    A::Vector{T},
    y_vars,
    r_objterm,
    c_objterm,
    big_m::Function,
) where T
    link_varsC = _create_link_varsC(mip_model, name, A, y_vars)
    oracle_solve = _build_default_oracle(mip_model, A, link_varsC, y_vars, c_objterm)
    return SubSolverBlCJuMP{T}(
        name,
        mip_model,
        A,
        link_varsC,
        y_vars,
        r_objterm,
        c_objterm,
        oracle_solve,
        big_m,
        Dict{Any,NamedTuple}(),
    )
end

function capacity_linking(sol::SubSolverBlCJuMP, a, params::SolverParam)
    return 1
end

function check(sol::SubSolverBlCJuMP, params::SolverParam)
    if !(objective_sense(sol.mip_model) == MIN_SENSE)
        error("The subsolver $(name(sol)) must be formulated as an minimization problem.")
    end

    if !(length(sol.link_varsC) == length(sol.y_vars))
        error(
            "In subsolver $(name(sol)), the number of linking variables does not match the number of second level variables in linking constraints! |link_vars| = $(length(sol.link_varsC)) and |y_vars| = $(length(sol.y_vars)).",
        )
    elseif !(length(sol.y_vars) == length(sol.A))
        error(
            "In subsolver $(name(sol)), the number of second level variables in linking constraints does not match the number of resources. |A| = $(length(sol.A)) and |y_vars| = $(length(sol.y_vars)).",
        )
    end

    for a in sol.A
        if !is_valid(sol.mip_model, sol.link_varsC[a])
            error("Linking variable for resource $(a) is not contained in subsolver $(name(sol))!")
        elseif !is_valid(sol.mip_model, sol.y_vars[a])
            error("Second level variable (y) for resource $(a) is not contained in subsolver $(name(sol))!")
        end
    end
end

function compute_lower_bound_master_contribution(sol::SubSolverBlCJuMP, params::SolverParam, time_limit)
    @objective(sol.mip_model, Min, sol.r_objterm)

    try
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/ubbound_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    solve_mip(sol, params, time_limit)
    check_solution_status(sol)

    @debug "The found lower bound solution for $(sol.name) is $(objective_value(sol.mip_model))."
    return objective_value(sol.mip_model)
end

function name(ss::SubSolverBlCJuMP)
    return ss.name
end

function separation!(sol::SubSolverBlCJuMP, sval, gvals, kvals::Dict, params::SolverParam, time_limit)
    obj_tolerance = 4
    var_non_zero_tolerance = 1e-3

    k_term = sum(kvals[a] * sol.link_varsC[a] for a in sol.A)
    k_penalty = sum(var_non_zero_tolerance * sol.link_varsC[a] for a in sol.A)
    new_obj = sol.c_objterm * gvals + sol.r_objterm + k_term + k_penalty
    @objective(sol.mip_model, Min, new_obj)

    try
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/sub_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    solve_mip(sol, params, time_limit)
    check_solution_status(sol)

    opt_obj = value(sol.c_objterm * gvals + sol.r_objterm + k_term)
    @debug "Optimal solution of sub $(sol.name) is $(opt_obj)"

    is_violate = Bool(!(sval < opt_obj + var_non_zero_tolerance))
    @debug "The violated status of sub $(sol.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    r = round(value(sol.r_objterm); digits=obj_tolerance)
    c = round(value(sol.c_objterm); digits=obj_tolerance)
    as = [a for a in sol.A if value(sol.link_varsC[a]) > 1 - var_non_zero_tolerance]
    return SubSolution(is_violate, r, c, opt_obj, as)
end

function separation_BlC!(sol::SubSolverBlCJuMP, sval, kvals::Dict, params::SolverParam, time_limit)
    obj_tolerance = 4
    var_non_zero_tolerance = 1e-3

    k_term = sum(kvals[a] * (1 - sol.link_varsC[a]) for a in sol.A)
    k_penalty = sum(var_non_zero_tolerance * (1 - sol.link_varsC[a]) for a in sol.A)
    new_obj = sol.c_objterm - k_term - k_penalty
    @objective(sol.mip_model, Max, new_obj)

    try
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/sub_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. Error message $err"
    end

    solve_mip(sol, params, time_limit)
    check_solution_status(sol)

    opt_obj = value(sol.c_objterm - k_term)
    @debug "Optimal solution of sub $(sol.name) is $(opt_obj)"

    is_violate = Bool(!(sval > opt_obj - var_non_zero_tolerance))
    @debug "The violated status of sub $(sol.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    r = round(value(sol.r_objterm); digits=obj_tolerance)
    c = round(value(sol.c_objterm); digits=obj_tolerance)
    as = [a for a in sol.A if value(sol.link_varsC[a]) < 1 - var_non_zero_tolerance]
    return SubSolution(is_violate, r, c, opt_obj, as)
end

function set_nthreads(sol::SubSolverBlCJuMP, n)
    set_attribute(sol.mip_model, MOI.NumberOfThreads(), capped_nthreads(n))
end

function set_singlethread(sol::SubSolverBlCJuMP)
    set_attribute(sol.mip_model, MOI.NumberOfThreads(), 1)
end

function supports_bilevel_subproblem_solver(sol::SubSolverBlCJuMP)
    return true
end

function solve_sub_for_x(sol::SubSolverBlCJuMP, xvals, params::SolverParam, time_limit)
    @objective(sol.mip_model, Min, sol.c_objterm)
    @constraint(sol.mip_model, fixc[a=sol.A], sol.link_varsC[a] == round(xvals[a]))

    try
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/subfix_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    @debug "Subproblem $(sol.name) MIP was adjusted. Start MIP solver with persistent BlC cuts."
    solve_mip(sol, params, time_limit)

    try
        status = termination_status(sol.mip_model)
        if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
            @debug "The Subproblem $(sol.name) was infeasible"
            return false, 0, 0, Dict()
        end

        check_solution_status(sol)

        y_vals = sol.y_vars isa AbstractDict ?
            Dict(a => value(sol.y_vars[a]) for a in keys(sol.y_vars)) :
            value.(sol.y_vars)
        osol = value(sol.c_objterm)
        osol_L1 = value(sol.r_objterm)
        @debug "We found a bilevel-feasible solution for subproblem $(sol.name) with value $(osol)."
        return true, osol, osol_L1, y_vals
    finally
        for a in sol.A
            delete(sol.mip_model, fixc[a])
        end
        unregister(sol.mip_model, :fixc)
    end
end

function verify_sub_for_x_optimistic(sol::SubSolverBlCJuMP, xvals, params::SolverParam, time_limit)
    start_time = time()

    @objective(sol.mip_model, Min, sol.c_objterm)
    @constraint(sol.mip_model, fixc[a=sol.A], sol.link_varsC[a] == round(xvals[a]))

    tie_constraint = nothing
    try
        solve_mip(sol, params, time_limit)

        status = termination_status(sol.mip_model)
        if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
            @debug "The Subproblem $(sol.name) was infeasible during optimistic verification."
            return false, 0, 0, Dict()
        end

        check_solution_status(sol)
        opt_cost = value(sol.c_objterm)

        remaining_time = time_limit - (time() - start_time)
        if remaining_time <= 0
            throw(TimeoutException("We hit the time limit while solving optimistic verification problem of bilevel subsolver $(sol.name)"))
        end

        tie_constraint = @constraint(sol.mip_model, sol.c_objterm <= opt_cost + 1e-6)
        @objective(sol.mip_model, Min, sol.r_objterm)
        solve_mip(sol, params, remaining_time)
        check_solution_status(sol)

        y_vals = sol.y_vars isa AbstractDict ?
            Dict(a => value(sol.y_vars[a]) for a in keys(sol.y_vars)) :
            value.(sol.y_vars)
        osol_L1 = value(sol.r_objterm)
        return true, opt_cost, osol_L1, y_vals
    finally
        if !isnothing(tie_constraint)
            delete(sol.mip_model, tie_constraint)
        end
        for a in sol.A
            delete(sol.mip_model, fixc[a])
        end
        unregister(sol.mip_model, :fixc)
    end
end

function check_solution_status(sol::SubSolverBlCJuMP)
    status = termination_status(sol.mip_model)
    if status == MOI.TIME_LIMIT
        @debug "We hit the time limit while solving MIP of bilevel subsolver $(sol.name)"
        throw(TimeoutException("We hit the time limit while solving MIP of bilevel subsolver $(sol.name)"))
    end

    if status == MOI.INFEASIBLE
        error("Submodel $(sol.name) was infeasible. This indicates inconsistent bilevel cuts or an invalid subproblem formulation.")
    elseif status == MOI.DUAL_INFEASIBLE
        error("The second level was unbounded what violates assumption on algorithm.")
    end

    if !(status == MOI.OPTIMAL)
        error("Submodel $(sol.name) could not be solved to optimality but terminated with status $status")
    end
end

function solve_mip(sol::SubSolverBlCJuMP, params::SolverParam, time_limit)
    inner_time_limit = time_limit
    var_non_zero_tolerance = 1e-3

    while true
        if inner_time_limit <= 0
            @debug "We hit the time limit while solving MIP of bilevel subsolver $(sol.name)."
            throw(TimeoutException("We hit the time limit while solving MIP of bilevel subsolver $(sol.name)."))
        end

        @debug "Setting time limit of bilevel subsolver $(sol.name) to $inner_time_limit"
        set_time_limit_sec(sol.mip_model, inner_time_limit)
        set_seed!(sol.mip_model, params.solver, get_seed(params))
        optimize!(sol.mip_model)
        inner_time_limit = inner_time_limit - solve_time(sol.mip_model)

        status = termination_status(sol.mip_model)
        @debug "Subproblem $(sol.name) MIP was solved with status $status"
        if status != MOI.OPTIMAL
            return
        end

        x_vals = Dict(a => round(value(sol.link_varsC[a])) for a in sol.A)
        x_key = key_master_sol(x_vals, sol.A)
        current_c = value(sol.c_objterm)
        oracle = _oracle_solution(sol, x_vals, x_key, params, inner_time_limit)

        if !oracle.feasible
            return
        end

        if current_c <= oracle.optL2 + var_non_zero_tolerance
            @debug "Current solution of $(sol.name) is bilevel-feasible for x=$(x_vals) with c=$(current_c) and oracle opt=$(oracle.optL2)."
            return
        end

        rhs = oracle.optL2 + sum(sol.big_m(a) * oracle.y_vals[a] * (1 - sol.link_varsC[a]) for a in sol.A)
        cut = @constraint(sol.mip_model, sol.c_objterm <= rhs)
        @debug "Added persistent BlC cut $(cut) for subsolver $(sol.name) at x=$(x_vals), candidate c=$(current_c), oracle c=$(oracle.optL2)."
    end
end

function _oracle_solution(sol::SubSolverBlCJuMP, x_vals, x_key, params::SolverParam, time_limit)
    if haskey(sol.cached_oracle_solutions, x_key)
        return sol.cached_oracle_solutions[x_key]
    end

    feasible, optL2, optL1, y_vals = sol.oracle_solve(x_vals, params, time_limit)
    result = (feasible=feasible, optL2=optL2, optL1=optL1, y_vals=y_vals)
    sol.cached_oracle_solutions[x_key] = result
    return result
end
