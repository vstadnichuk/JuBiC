
using JuMP
using BenchmarkTools
using CSV
import MathOptInterface as MOI
using Base.Threads


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
    if param.bigMwithLC
        new_stat!(param.stats, "Solver", "GBCLagSolver")  # Set other name as essentially a new solver
    else
        new_stat!(param.stats, "Solver", "GBCSolver")
    end
    new_stat!(param.stats, "time_limit", param.runtime)
    new_stat!(param.stats, "NSub", length(subs))
    new_stat!(param.stats, "NFeasCuts", 0)
    new_stat!(param.stats, "NOptCuts", 0)
    new_stat!(param.stats, "BlCLagCuts", 0)  # Number of added BlC constraints. Only added if their big M are generated in the subroutine
    new_stat!(param.stats, "NBigMlagCuts", 0)  # number of Lagrangian cuts computed to obtain better big M coef. 
    new_stat!(param.stats, "SepaTime", 0)  # time spend in separator
    new_stat!(param.stats, "SepaTimeCut", 0)  # time spend in separator for generating cuts only
    new_stat!(param.stats, "ConnectorLPTimeLP", 0.0)  # time spent solving the connector LP relaxation
    new_stat!(param.stats, "ConnectorLPTimePricing", 0.0)  # time spent in subsolver pricing / separation calls
    new_stat!(param.stats, "ConnectorLPTimePareto", 0.0)  # time spent in pareto refinement inside connectors
    new_stat!(param.stats, "ConnectorLPIterations", 0)  # number of connector LP resolve/pricing iterations
    new_stat!(param.stats, "ConnectorApproximation", string(param.connector_approximation))
    new_stat!(param.stats, "ConnectorPricingMaxError", 0.0)
    new_stat!(param.stats, "NInexactPricingCalls", 0)
    new_stat!(param.stats, "UsedInexactPricing", false)
    param.stats.data["ConnectorPricingMaxErrorBySub"] = Dict{String,Float64}()
    new_stat!(param.stats, "parallel_separation", param.parallel_separation)
    master_threads = resolve_nthreads!(param.stats, "threads_master", param.threads_master; context="the master MIP")
    sub_threads = resolve_nthreads!(param.stats, "threads_sub_con", param.threads_sub_con; context="the subproblem solvers")
    validate_parallel_subsolver_threads!(param.parallel_separation, sub_threads; context="Gurobi-backed subproblem solves")
    if param.parallel_separation
        parallel_workers = _resolve_parallel_workers!(param.stats, sub_threads)
        new_stat!(param.stats, "parallel_connector_workers_used", parallel_workers)
    end

    # do some initail checks for master and sub solvers
    @debug "Doing some checks if master and sub were created correctly."
    check(master, param)
    for sub in subs
        check(sub, param)
    end

    # Preserve the direct first-level objective so the final master incumbent
    # can be combined with exact fixed-x follower evaluations.
    base_master_objective = objective_function(master.model)

    # add subObj variables, one for each sub_problem, to the master, and add them to objective
    @debug "Start with initialization of GBC solver."
    names = [name(sub) for sub in subs]
    subObj = @variable(master.model, subObj[names])
    @objective(master.model, Min, objective_function(master.model) + sum(subObj))

    # add partial decomposition to master if requested
    if !isnothing(master.partial_decomposition)
        master.partial_decomposition(master.model, subObj)
    end

    # configure solver logging before preprocessing and final export
    if should_write_output_logs(param)
        try
            set_silent(master.model)
            set_optimizer_attribute(master.model, "LogFile", param.output_folder_path*"/gbc_mip_log.txt")
        catch err 
            @error "Could not set log file for folder $(param.output_folder_path). Error is $err"
        end
    else
        set_silent(master.model)
    end

    # build the sub LPs for Benders subroutine 
    clps, runtime_init = init_connectorLPs(subs, master.link_vars, subObj, param) 
    new_stat!(param.stats, "runtime_preprocessingGBC", runtime_init)

    # debug output: export the master only after preprocessing so stored subObj bounds match the live model
    if should_write_output_logs(param)
        try
            write_to_file(
                master.model,
                param.output_folder_path * "/master.$(param.file_format_output)",
            )
        catch err 
            @error "Could not write model to file for folder $(param.output_folder_path). Error is $err"
        end
    end

    # set time limit and number of threads
    true_runtime = param.runtime - runtime_init
    set_time_limit_sec(master.model, true_runtime)
    set_attribute(master.model, MOI.NumberOfThreads(), master_threads)
    set_seed!(master.model, param.solver, get_seed(param))
    for sub in subs
        if param.parallel_separation
            set_singlethread(sub)
        else
            set_nthreads(sub, sub_threads)
        end
    end

    # add callback to master and solve 
    msol_cuts_mapping = Dict()  # a mapping of master solution to found lazy constraints
    msol_cuts_mapping_blc = Dict()  # a mapping of master solution to found lazy blc constraints. They are only generated if BlC coef. are automatically computed in subroutine
    msol_subobj_mapping = Dict()  # for each master solution, store per subproblem which subObj values were already separated
    msol_certificate_mapping = Dict()  # final pricing/connector certificate per master solution and follower
    solve_start_time = time()
    if true_runtime > 0
        @debug "Finished model construction. Now proceeding to optimization process with GBC. Remaining runtime is $true_runtime"
        set_attribute(master.model, MOI.LazyConstraintCallback(), cb -> gbc_callback_function(cb, master, names, clps, subObj, msol_cuts_mapping, msol_cuts_mapping_blc, msol_subobj_mapping, msol_certificate_mapping, param, solve_start_time, true_runtime))
    else
        @debug "We do not add any callbacks to GBCSolver because preprocessing consumed the available runtime."
        new_stat!(param.stats, "GBCStatus", "Timelimit")
    end

    try 
        optimize!(master.model)
    catch e
        if (e isa TimeoutException)
            status = _gbc_timeout_status(param, true_runtime)
            if status == "Timeout_Submodel"
                @warn "GBC stopped because a submodel/connector timed out after receiving essentially the full remaining solve budget."
            else
                @warn "GBC reached the global runtime limit during submodel/connector separation. Reporting this as a normal timelimit."
            end
            param.stats.data["GBCStatus"] = status
            param.stats.data["Opt_status_override"] = status
        elseif (e isa NumericalIssueException)
            @error "GBCSolver stopped due to a detected numerical issue: $(e.message)"
            param.stats.data["GBCStatus"] = e.status
            param.stats.data["Opt_status_override"] = e.status
        else
            @error "GBCsolver suffered an error: $e"
            @error stacktrace(catch_backtrace())
            #showerror(stdout, e, catch_backtrace())
            param.stats.data["GBCStatus"] = "Terminate"
            param.stats.data["Opt_status_override"] = "Terminate"
            #rethrow(e)  
        end

        # collect whatever information is still available from the master after the interrupted solve
        try
            if primal_status(master.model) == MOI.FEASIBLE_POINT
                mobj = objective_value(master.model)
                xsol = round_master_solution(
                    Dict(a => value(master.link_vars[a]) for a in master.A),
                )
                print_solution_to_file(mobj, xsol, param)
                new_stat!(param.stats, "Opt", mobj)
            end
        catch err
            @warn "Could not recover a master incumbent after interrupted GBC solve: $(sprint(showerror, err))"
        end
        try
            set_optimization_status_stats(termination_status(master.model), param)
        catch
            set_optimization_status_stats(MOI.TIME_LIMIT, param)
        end
        try
            new_stat!(param.stats, "runtime", solve_time(master.model))
        catch
            new_stat!(param.stats, "runtime", param.runtime)
        end
        try
            new_stat!(param.stats, "gap", JuMP.relative_gap(master.model))
        catch
        end
        try
            new_stat!(param.stats, "BNodes", MOI.get(master.model, MOI.NodeCount()))
        catch
        end
    else
        # print solution and collected data
        print_collected_cuts(param, msol_cuts_mapping)
        print_collected_cuts(param, msol_cuts_mapping_blc; filename="mastercuts_blc.txt")
        if termination_status(master.model) == MOI.OPTIMAL || termination_status(master.model) == MOI.LOCALLY_SOLVED || termination_status(master.model) == MOI.TIME_LIMIT
            if primal_status(master.model) == MOI.FEASIBLE_POINT
                try
                    mobj = objective_value(master.model)
                    xsol = round_master_solution(
                        Dict(a => value(master.link_vars[a]) for a in master.A),
                    )
                    @debug "The master objective is $(mobj) and solution is $(xsol)."
                    print_solution_to_file(mobj, xsol, param)
                    new_stat!(param.stats, "Opt", mobj)

                    # print full MIP solution to file
                    if should_write_output_logs(param)
                        solution = Dict(JuMP.name(x) => JuMP.value(x) for x in all_variables(master.model))
                        write(param.output_folder_path*"/full_master_solution.json", JSON.json(solution))
                    end
                catch err
                    @warn "Could not read back the final GBC master incumbent after optimization: $(sprint(showerror, err))"
                    if haskey(param.stats.data, "GBCStatus")
                        param.stats.data["GBCStatus"] = "OptimizeNotCalled"
                    else
                        new_stat!(param.stats, "GBCStatus", "OptimizeNotCalled")
                    end
                end
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

    _finalize_gbc_objective_interval!(
        master,
        names,
        clps,
        base_master_objective,
        msol_certificate_mapping,
        param,
        solve_start_time,
        true_runtime,
    )

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
    set_attribute(myLP, MOI.NumberOfThreads(), _connector_thread_count(parameter))

    # TODO: This parameter combination seems to fix some numeric issues. Seems to have only necglectable impact on runtime
    if parameter.solver isa GurobiSolver
        set_optimizer_attribute(myLP, "NumericFocus", 3)
        set_optimizer_attribute(myLP, "CrossoverBasis", 1)
        set_optimizer_attribute(myLP, "Method", 2)
        set_optimizer_attribute(myLP, "DualReductions", 0)
        set_optimizer_attribute(myLP, "BarHomogeneous", 1)
    end

    # disable output of LP
    set_silent(myLP)

    # calculate a bound on the big M coefficients
    lbm = compute_lower_bound_master_contribution(sub, parameter, timelimit)
    @debug "For subpoblem $(name(sub)), the found lower bound to the master objective contribution is $(lbm). With this, construction of ConnectorLP finished."

    # set computed lower bound for sub_problem objective variables
    set_lower_bound(subObjvar, lbm)

    # add generator for BlC cuts if better cuts are requested
    blc_generator = nothing
    if parameter.bigMwithLC
        blc_generator = ConnectorLP_BlC(parameter, sub.A, link_vars_master, sub)
    end

    # build ConnectorLP obj
    return ConnectorLP(myLP, sub.A, link_vars_master, sub, lbm, blc_generator, Vector{ConSubsolCut}(), parameter.g_round_digit, Dict{Symbol,Any}())
end



function _normalize_subobj_cache_value(v::Real)
    return round(Float64(v), digits=6)
end

function _resolve_parallel_workers!(stats::RunStats, requested_workers::Integer)
    used = min(max(1, Int(requested_workers)), max(1, Threads.nthreads()))
    if used < requested_workers
        @warn "Requested $(requested_workers) parallel connector workers, but the Julia process only has $(Threads.nthreads()) thread(s). JuBiC will use $(used) worker(s) for parallel connector separation."
    end
    return used
end

function _connector_thread_count(params::GBCparam)
    return params.parallel_separation ? 1 : used_nthreads(params.stats, "threads_sub_con")
end

function _local_gbc_param(params::GBCparam)
    local_stats = RunStats()
    new_stat!(local_stats, "threads_master_used", get(params.stats.data, "threads_master_used", 1))
    new_stat!(local_stats, "threads_sub_con_used", get(params.stats.data, "threads_sub_con_used", 1))
    new_stat!(local_stats, "ConnectorLPTimeLP", 0.0)
    new_stat!(local_stats, "ConnectorLPTimePricing", 0.0)
    new_stat!(local_stats, "ConnectorLPTimePareto", 0.0)
    new_stat!(local_stats, "ConnectorLPIterations", 0)
    return GBCparam(
        params.solver,
        params.debbug_out,
        params.output_folder_path,
        params.file_format_output,
        local_stats,
        params.runtime,
        params.seed,
        params.threads_master,
        params.threads_sub_con,
        params.parallel_separation,
        params.pareto,
        params.warmstart,
        params.bigMwithLC,
        params.trim_coeff,
        params.infinity_num,
        params.g_round_digit,
        params.integer_obj,
        params.pareto_band_tolerance,
        params.blc_pareto_band_tolerance,
        params.connector_add_current_solution_cut,
        params.subsolver_numerical_preprocessing,
        params.connector_approximation,
    )
end

function _format_parallel_gbc_task_error(subname::AbstractString, err, bt)
    io = IOBuffer()
    print(io, "Parallel GBC separator failed for subproblem ", subname, ". ")
    showerror(io, err, bt)
    return String(take!(io))
end

function _merge_parallel_gbc_stats!(target::RunStats, local_stats::RunStats)
    for key in ("ConnectorLPTimeLP", "ConnectorLPTimePricing", "ConnectorLPTimePareto", "ConnectorLPIterations")
        if haskey(local_stats.data, key)
            add_stat!(target, key, local_stats.data[key])
        end
    end
    if haskey(local_stats.data, "NOptCutValidationWarnings")
        if haskey(target.data, "NOptCutValidationWarnings")
            add_stat!(target, "NOptCutValidationWarnings", local_stats.data["NOptCutValidationWarnings"])
        else
            new_stat!(target, "NOptCutValidationWarnings", local_stats.data["NOptCutValidationWarnings"])
        end
    end
    if haskey(local_stats.data, "OptCutValidationUsers")
        offenders = get!(target.data, "OptCutValidationUsers", String[])
        for uname in local_stats.data["OptCutValidationUsers"]
            uname in offenders || push!(offenders, uname)
        end
    end
    if get(local_stats.data, "Opt_status_override", nothing) == "Numerics"
        target.data["Opt_status_override"] = "Numerics"
        target.data["GBCStatus"] = "Numerics"
    end
    return nothing
end

function _cache_subobj_value!(mapping::Dict, msolkey, subname, value)
    if !haskey(mapping, msolkey)
        mapping[msolkey] = Dict{String, Set{Float64}}()
    end
    per_sub = mapping[msolkey]
    if !haskey(per_sub, subname)
        per_sub[subname] = Set{Float64}()
    end
    push!(per_sub[subname], _normalize_subobj_cache_value(value))
end

function _has_cached_subobj_value(mapping::Dict, msolkey, subname, value)
    haskey(mapping, msolkey) || return false
    per_sub = mapping[msolkey]
    haskey(per_sub, subname) || return false
    return _normalize_subobj_cache_value(value) in per_sub[subname]
end

function _remaining_gbc_callback_time(solve_start_time::Real, true_runtime::Real)
    return true_runtime - (time() - solve_start_time)
end

function _record_gbc_submodel_time_budget!(stats::RunStats, remaining_time::Real, true_runtime::Real)
    budget = max(0.0, Float64(remaining_time))
    total = max(Float64(true_runtime), eps(Float64))
    stats.data["GBCLastSubmodelTimeLimit"] = budget
    stats.data["GBCLastSubmodelTimeLimitShare"] = budget / total
    return nothing
end

function _gbc_timeout_status(param::GBCparam, true_runtime::Real)
    budget = get(param.stats.data, "GBCLastSubmodelTimeLimit", param.runtime)
    total = max(Float64(true_runtime), eps(Float64))
    full_budget_threshold = 0.9 * total
    return budget >= full_budget_threshold ? "Timeout_Submodel" : "Timelimit"
end

"""Return the serializable certificate left by the most recent connector solve."""
function _connector_certificate(con::ConnectorLP)
    state = con.numeric_state
    return (
        raw_value=Float64(get(state, :raw_connector_value, NaN)),
        cut_value=Float64(get(state, :local_cut_value, NaN)),
        pricing_upper=Float64(get(state, :pricing_upper_bound, NaN)),
        pricing_lower=Float64(get(state, :pricing_lower_bound, NaN)),
        error=Float64(get(state, :pricing_error, 0.0)),
        fixed_x_master_contribution=Float64(
            get(state, :fixed_x_master_contribution, NaN),
        ),
        pricing_is_exact=Bool(get(state, :pricing_is_exact, true)),
        fixed_x_feasible=Bool(get(state, :fixed_x_feasible, false)),
    )
end

"""
Compute the certified interval for the original bilevel optimum.

For safe underestimation, the master objective bound is already a lower bound
on the original optimum. For tight overestimation, each retained cut for
follower `k` may exceed its value function by at most `delta_kj`; subtracting
`sum(k, maximum(j, delta_kj))` from the heuristic-master bound restores a valid
lower bound. At the final binary incumbent, GBC attempts the optional
`solve_sub_for_x_optimistic` method for every follower. When all followers
support it, the interval width is bounded by the master MIP gap plus the
relevant connector pricing-error term. An unsupported follower falls back to
the exact but arbitrarily tie-broken fixed-`x` response cached during
separation. That fallback remains a valid incumbent upper bound, but its
additional follower tie-breaking error is not controlled by the pricing bound.
"""
function _finalize_gbc_objective_interval!(
    master::Master,
    sub_names,
    connectors,
    base_master_objective,
    certificate_mapping::Dict,
    param::GBCparam,
    solve_start_time::Real,
    solve_runtime::Real,
)
    model = master.model
    primal_status(model) == MOI.FEASIBLE_POINT || return nothing

    master_bound = try
        Float64(objective_bound(model))
    catch
        -Inf
    end
    master_incumbent = try
        Float64(objective_value(model))
    catch
        Inf
    end
    param.stats.data["MasterObjective"] = master_incumbent
    param.stats.data["MasterObjectiveBound"] = master_bound
    param.stats.data["MasterAbsoluteGap"] = master_incumbent - master_bound
    # Compatibility aliases for existing benchmark post-processing.
    param.stats.data["HeuristicMasterObjective"] = master_incumbent
    param.stats.data["HeuristicMasterObjectiveBound"] = master_bound
    param.stats.data["HeuristicMasterAbsoluteGap"] = master_incumbent - master_bound

    global_error = 0.0
    if param.connector_approximation == CONNECTOR_OVERESTIMATION
        by_sub = get(
            param.stats.data,
            "ConnectorPricingMaxErrorBySub",
            Dict{String,Float64}(),
        )
        global_error = sum(get(by_sub, String(name), 0.0) for name in sub_names)
    end
    interval_lower = master_bound - global_error

    # Normalize solver tolerances before the values are used as a cache key.
    xvals = round_master_solution(
        Dict(a => Float64(value(master.link_vars[a])) for a in master.A),
    )
    msolkey = key_master_sol(xvals, master.A)
    final_certificates = get(certificate_mapping, msolkey, Dict{String,Any}())
    final_error = 0.0
    connector_upper_sum = 0.0
    cache_complete = true
    for subname in sub_names
        cert = get(final_certificates, String(subname), nothing)
        if isnothing(cert) || !cert.fixed_x_feasible
            cache_complete = false
            break
        end
        connector_upper_sum += cert.raw_value
        final_error += cert.error
    end
    # Prefer a lexicographically optimistic fixed-x solution for the final
    # upper endpoint. This optional solve is intentionally isolated from cut
    # generation so custom subsolvers need not implement it.
    optimistic_risk_sum = 0.0
    optimistic_complete = true
    optimistic_used = false
    fallback_followers = String[]
    optimistic_by_sub = Dict{String,Bool}()
    connector_by_name = Dict(String(name(con.sub_solver)) => con.sub_solver for con in connectors)
    optimistic_start = time()
    for subname in sub_names
        subkey = String(subname)
        subsolver = get(connector_by_name, subkey, nothing)
        contribution = nothing
        if !isnothing(subsolver) && applicable(
            solve_sub_for_x_optimistic,
            subsolver,
            xvals,
            param,
            1.0,
        )
            remaining = solve_runtime - (time() - solve_start_time)
            if remaining > 0
                try
                    found, _, optimistic_contribution, _ = solve_sub_for_x_optimistic(
                        subsolver,
                        xvals,
                        param,
                        remaining,
                    )
                    found || error(
                        "Optimistic fixed-x evaluation declared final incumbent x infeasible for follower $(subkey).",
                    )
                    contribution = Float64(optimistic_contribution)
                    optimistic_used = true
                    optimistic_by_sub[subkey] = true
                catch err
                    if !(err isa TimeoutException)
                        rethrow()
                    end
                    @warn "Optimistic fixed-x evaluation timed out for follower $(subkey); using its cached arbitrary follower-optimal response."
                end
            end
        end

        if isnothing(contribution)
            optimistic_complete = false
            optimistic_by_sub[subkey] = false
            push!(fallback_followers, subkey)
            cert = get(final_certificates, subkey, nothing)
            if isnothing(cert) || !cert.fixed_x_feasible
                cache_complete = false
                continue
            end
            contribution = Float64(cert.fixed_x_master_contribution)
        end
        optimistic_risk_sum += contribution
    end
    optimistic_time = time() - optimistic_start

    first_level_value = try
        Float64(value(base_master_objective))
    catch
        NaN
    end
    upper_complete = optimistic_complete || cache_complete
    interval_upper = upper_complete ? first_level_value + optimistic_risk_sum : Inf
    connector_upper = cache_complete ? first_level_value + connector_upper_sum : Inf

    pricing_width_error = if param.connector_approximation == CONNECTOR_OVERESTIMATION
        global_error
    else
        cache_complete ? final_error : Inf
    end
    master_gap = max(0.0, master_incumbent - master_bound)
    width_bounded_by_pricing = optimistic_complete && cache_complete &&
        isfinite(master_gap) && isfinite(pricing_width_error)
    width_bound = width_bounded_by_pricing ? master_gap + pricing_width_error : Inf

    param.stats.data["ConnectorApproximationErrorBound"] = global_error
    param.stats.data["FinalConnectorErrorSum"] = cache_complete ? final_error : Inf
    param.stats.data["FinalConnectorUpperBound"] = connector_upper
    param.stats.data["FinalOptimisticEvaluationUsed"] = optimistic_used
    param.stats.data["FinalOptimisticEvaluationComplete"] = optimistic_complete
    param.stats.data["FinalOptimisticEvaluationBySub"] = optimistic_by_sub
    param.stats.data["FinalOptimisticEvaluationFallbackFollowers"] = fallback_followers
    param.stats.data["FinalOptimisticEvaluationTime"] = optimistic_time
    param.stats.data["IncumbentObjectiveUpperBound"] = interval_upper
    param.stats.data["OptimisticIncumbentObjective"] = optimistic_complete ? interval_upper : Inf
    # Historical compatibility alias for the selected exact fixed-x response.
    param.stats.data["ExactIncumbentObjective"] = interval_upper
    param.stats.data["ObjectiveIntervalLower"] = interval_lower
    param.stats.data["ObjectiveIntervalUpper"] = interval_upper
    param.stats.data["ObjectiveIntervalWidth"] = interval_upper - interval_lower
    param.stats.data["ObjectiveIntervalCertified"] = isfinite(interval_lower) && isfinite(interval_upper)
    param.stats.data["ObjectiveIntervalWidthBoundedByPricingError"] = width_bounded_by_pricing
    param.stats.data["ObjectiveIntervalWidthBound"] = width_bound

    if !upper_complete
        @warn "The final GBC master incumbent was not fully represented in the connector certificate cache. JuBiC does not perform a fallback fixed-x solve and reports an infinite objective-interval upper endpoint."
    end
    if !optimistic_complete
        @warn "Followers $(fallback_followers) do not provide a completed optimistic fixed-x evaluation. The reported objective interval remains valid, but its width is not bounded solely by the master gap and certified subsolver pricing errors because follower tie-breaking error is unknown."
    end

    used_inexact = Bool(get(param.stats.data, "UsedInexactPricing", false))
    solution_type = used_inexact ? "Heuristic" : "Exact"
    master_optimal = termination_status(model) in (MOI.OPTIMAL, MOI.LOCALLY_SOLVED)
    result_status = if master_optimal
        used_inexact ? "HeuristicOptimal" : "Optimal"
    else
        base_status = string(get(param.stats.data, "Opt_status", termination_status(model)))
        used_inexact ? "Heuristic_$(base_status)" : base_status
    end
    param.stats.data["GBCSolutionType"] = solution_type
    param.stats.data["GBCResultStatus"] = result_status
    param.stats.data["MasterSolvedToOptimality"] = master_optimal
    if master_optimal && !haskey(param.stats.data, "Opt_status_override")
        param.stats.data["GBCStatus"] = result_status
        param.stats.data["Opt_status"] = result_status
    end
    return nothing
end

function gbc_callback_function(cb_data, master::Master, sub_names, clps, subObj, msol_cuts_mapping::Dict, msol_cuts_mapping_blc::Dict, msol_subobj_mapping::Dict, msol_certificate_mapping::Dict, parameter::GBCparam, solve_start_time::Real, true_runtime::Real)
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
            x_vals = round_master_solution(
                Dict(
                    a => callback_value(cb_data, master.link_vars[a]) for a in master.A
                ),
            )
            @debug "Current values of the master linking variables are $(x_vals)"
            lazy = []
            lazy_blc = []

            # Because of multi thread (and some start trouble I cannot explain) we resolve the subproblems for same master solution multiple times. 
            # To avoid this, we save the found cuts for each master solution and fall back on them before resolving the sub_problem
            msolkey = key_master_sol(x_vals, master.A)
            #@debug "The current master_key=$(msolkey) and the saved solutions are $(msol_cuts_mapping).")
            if haskey(msol_cuts_mapping, msolkey)
                # if we already solved this sub_problem, we just recover the found lazy constraints 
                lazy = msol_cuts_mapping[msolkey]
                lazy_blc = msol_cuts_mapping_blc[msolkey]
                @debug "We recovered the lazy cuts for solution $(x_vals). We only rerun separators for subproblems whose current subObj value was not tested yet for this master solution."
            else
                # solve each sub and add cut to list of cuts "lazy"
                @debug "We found no saved lazy cuts for current solution and solve subproblems now."
                msol_cuts_mapping[msolkey] = lazy
                msol_cuts_mapping_blc[msolkey] = lazy_blc
            end

            pending = Tuple{Int,Any,String,Float64}[]
            for (idx, con) in enumerate(clps)
                subname = name(con)
                current_subobj = subObj_val[subname]
                if _has_cached_subobj_value(msol_subobj_mapping, msolkey, subname, current_subobj)
                    @debug "For master solution $(x_vals) and sub $(subname), the current subObj value $(current_subobj) was already checked. Skip resolving this separator."
                    continue
                end
                push!(pending, (idx, con, subname, current_subobj))
            end

            results = Vector{Any}(undef, length(pending))
            task_errors = Vector{Any}(undef, length(pending))
            fill!(task_errors, nothing)
            if parameter.parallel_separation && length(pending) > 1
                sem = Base.Semaphore(get(parameter.stats.data, "parallel_connector_workers_used", 1))
                @sync for (res_idx, (_, con, subname, current_subobj)) in enumerate(pending)
                    Threads.@spawn begin
                        Base.acquire(sem)
                        try
                            local_param = _local_gbc_param(parameter)
                            result_ref = Ref{Any}(nothing)
                            cuttime = @elapsed begin
                                remaining_time = _remaining_gbc_callback_time(solve_start_time, true_runtime)
                                _record_gbc_submodel_time_budget!(parameter.stats, remaining_time, true_runtime)
                                if remaining_time <= 0
                                    throw(TimeoutException("GBC callback reached the global time limit before separating subproblem $(subname)."))
                                end
                                feas, cut, bigMcut, pobj = genBenders_cut!(con, x_vals, local_param, remaining_time)
                                result_ref[] = (
                                    subname=subname,
                                    current_subobj=current_subobj,
                                    feas=feas,
                                    cut=cut,
                                    bigMcut=bigMcut,
                                    pobj=pobj,
                                    certificate=_connector_certificate(con),
                                    stats=local_param.stats,
                                )
                            end
                            base_result = result_ref[]
                            results[res_idx] = (; base_result..., cuttime=cuttime)
                        catch err
                            task_errors[res_idx] = (
                                subname=subname,
                                current_subobj=current_subobj,
                                err=err,
                                bt=stacktrace(catch_backtrace()),
                            )
                        finally
                            Base.release(sem)
                        end
                    end
                end

                failures = [failure for failure in task_errors if !isnothing(failure)]
                if !isempty(failures)
                    for failure in failures
                        @error _format_parallel_gbc_task_error(failure.subname, failure.err, failure.bt)
                    end
                    throw(first(failures).err)
                end
            else
                for (res_idx, (_, con, subname, current_subobj)) in enumerate(pending)
                    local_param = parameter.parallel_separation ? _local_gbc_param(parameter) : parameter
                    result_ref = Ref{Any}(nothing)
                    cuttime = @elapsed begin
                        remaining_time = _remaining_gbc_callback_time(solve_start_time, true_runtime)
                        _record_gbc_submodel_time_budget!(parameter.stats, remaining_time, true_runtime)
                        if remaining_time <= 0
                            throw(TimeoutException("GBC callback reached the global time limit before separating subproblem $(subname)."))
                        end
                        feas, cut, bigMcut, pobj = genBenders_cut!(con, x_vals, local_param, remaining_time)
                        result_ref[] = (
                            subname=subname,
                            current_subobj=current_subobj,
                            feas=feas,
                            cut=cut,
                            bigMcut=bigMcut,
                            pobj=pobj,
                            certificate=_connector_certificate(con),
                            stats=(parameter.parallel_separation ? local_param.stats : nothing),
                        )
                    end
                    base_result = result_ref[]
                    results[res_idx] = (; base_result..., cuttime=cuttime)
                end
            end

            for result in results
                if !isnothing(result.stats)
                    _merge_parallel_gbc_stats!(parameter.stats, result.stats)
                end
                _cache_subobj_value!(msol_subobj_mapping, msolkey, result.subname, result.current_subobj)
                per_solution_certificates = get!(msol_certificate_mapping, msolkey, Dict{String,Any}())
                per_solution_certificates[result.subname] = result.certificate
                if isfinite(result.certificate.error)
                    parameter.stats.data["ConnectorPricingMaxError"] = max(
                        Float64(get(parameter.stats.data, "ConnectorPricingMaxError", 0.0)),
                        Float64(result.certificate.error),
                    )
                else
                    parameter.stats.data["ConnectorPricingMaxError"] = Inf
                end
                if !isfinite(result.certificate.error) || result.certificate.error > 1e-9
                    parameter.stats.data["UsedInexactPricing"] = true
                    add_stat!(parameter.stats, "NInexactPricingCalls", 1)
                end
                add_stat!(parameter.stats, "SepaTimeCut", result.cuttime)

                if result.feas
                    cutfeas = @build_constraint(result.cut >= 1)
                    @debug "Adding feasibility cut $(cutfeas) to the master problem for sub $(result.subname)."
                    add_stat!(parameter.stats, "NFeasCuts", 1)
                    push!(lazy, cutfeas)
                else
                    if result.current_subobj + 1e-6 < result.pobj
                        cutopt = @build_constraint(result.cut <= subObj[result.subname])
                        @debug "Adding optimality cut $(cutopt) to the master problem for sub $(result.subname)."
                        add_stat!(parameter.stats, "NOptCuts", 1)
                        push!(lazy, cutopt)

                        # Only retained overestimating cuts restrict the master and
                        # therefore contribute to its global correction E.  Pricing
                        # calls that produce no cut must not enlarge that correction.
                        if parameter.connector_approximation == CONNECTOR_OVERESTIMATION
                            max_error_by_sub = get!(
                                parameter.stats.data,
                                "ConnectorPricingMaxErrorBySub",
                                Dict{String,Float64}(),
                            )
                            max_error_by_sub[result.subname] = max(
                                get(max_error_by_sub, result.subname, 0.0),
                                Float64(result.certificate.error),
                            )
                        end

                        if !isnothing(master.objL2) && !isnothing(result.bigMcut)
                            cutopt_blc = @build_constraint(master.objL2[result.subname] <= result.bigMcut)
                            @debug "Adding in addition to optimality cut also the BlC constraint $(cutopt_blc) to the master problem for sub $(result.subname)."
                            add_stat!(parameter.stats, "BlCLagCuts", 1)
                            push!(lazy_blc, cutopt_blc)
                        end
                    end
                end
            end

            # add lazy cuts to master model
            map(cu -> MOI.submit(master.model, MOI.LazyConstraint(cb_data), cu), lazy)
            map(cu -> MOI.submit(master.model, MOI.LazyConstraint(cb_data), cu), lazy_blc)

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
