
using JuMP
using BenchmarkTools
using CSV
using JSON
import MathOptInterface as MOI
using Base.Threads

const GBC_CONNECTOR_TIMEOUT_BUFFER = 1.0

function _request_gurobi_master_termination!(cb_data, reason::AbstractString)
    try
        Gurobi.GRBterminate(cb_data.model)
        @debug "Requested native Gurobi master termination: $(reason)"
    catch err
        @error "Could not request native Gurobi master termination after $(reason)" exception=(err, catch_backtrace())
    end
    return nothing
end

function _write_numerical_diagnostic(param, err::NumericalIssueException, bt)
    # Keep ordinary runs lightweight: write a structured diagnostic only when
    # numerical instability aborts the GBC solve.
    try
        diagnostic = Dict{String,Any}(
            "exception_type" => string(typeof(err)),
            "status" => err.status,
            "message" => err.message,
            "context" => err.context,
            "stacktrace" => string.(bt),
        )
        open(joinpath(param.output_folder_path, "numeric_termination.json"), "w") do io
            write(io, JSON.json(diagnostic, 2))
        end
    catch diagnostic_error
        @warn "Could not write numerical termination diagnostics: $(diagnostic_error)"
    end
end


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
    new_stat!(param.stats, "parallel_separation", param.parallel_separation)
    master_threads = resolve_nthreads!(param.stats, "threads_master", param.threads_master; context="the master MIP")
    sub_threads = resolve_nthreads!(param.stats, "threads_sub_con", param.threads_sub_con; context="the subproblem solvers")
    if param.parallel_separation
        parallel_workers = _resolve_parallel_workers!(param.stats, sub_threads)
        new_stat!(param.stats, "parallel_connector_workers_used", parallel_workers)
    else
        parallel_workers = 1
    end

    # Gurobi environments are not thread-safe.  Initialize the complete pool
    # before any callback can spawn separation workers; the connector models
    # are then assigned to these environments during construction below.
    if param.solver isa GurobiSolver
        ensure_worker_envs!(param.solver, parallel_workers)
    end

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

    # Use one absolute deadline for preprocessing, master optimization, and
    # all callback-side separation calls.  A callback may be entered many
    # times, so passing `param.runtime` directly to each callback would reset
    # the effective time limit on every master incumbent.
    deadline = time() + param.runtime

    # build the sub LPs for Benders subroutine
    clps, runtime_init = init_connectorLPs(subs, master.link_vars, subObj, param)
    new_stat!(param.stats, "runtime_preprocessingGBC", runtime_init)

    # Models supplied by callers may have been created with Gurobi's default
    # environment.  Rebind the master and follower models before optimization
    # so that the master has its own environment and each connector's complete
    # worker-side model set uses the connector's assigned worker environment.
    if param.solver isa GurobiSolver
        _assign_gurobi_worker_models!(master, subs, clps, param.solver, parallel_workers)
    end

    # Configure logging after any optimizer rebinding, since replacing a
    # JuMP optimizer clears optimizer-specific attributes such as LogFile.
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
    true_runtime = max(0.0, deadline - time())
    set_time_limit_sec(master.model, true_runtime)
    # Reserve a small margin for callback unwinding and model cleanup.  The
    # callback will request native termination if it is entered in this
    # buffer, and no ConnectorLP separation will be started.
    if true_runtime <= GBC_CONNECTOR_TIMEOUT_BUFFER
        param.stats.data["GBCStatus"] = "Timeout_Submodel"
        param.stats.data["Opt_status_override"] = "Timeout_Submodel"
    end
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
    callback_lock = ReentrantLock()
    worker_locks = [ReentrantLock() for _ in 1:parallel_workers]
    if true_runtime > 0
        @debug "Finished model construction. Now proceeding to optimization process with GBC. Remaining runtime is $true_runtime"
        set_attribute(master.model, MOI.LazyConstraintCallback(), cb -> begin
            # Gurobi may invoke callbacks from multiple master threads. The
            # callback owns shared cut caches and mutable ConnectorLP objects;
            # serialize callback entry while retaining parallel separation of
            # distinct connectors inside one callback.
            lock(callback_lock)
            try
                gbc_callback_function(cb, master, names, clps, subObj, msol_cuts_mapping, msol_cuts_mapping_blc, msol_subobj_mapping, param, deadline, worker_locks)
            finally
                unlock(callback_lock)
            end
        end)
    else
        @debug "We do not add any callbacks to GBCSolver as we run into a time out during the preprocessing."
        new_stat!(param.stats, "GBCStatus", "Timeout_Submodel")
    end

    try 
        optimize!(master.model)
    catch e
        if (e isa TimeoutException)
            @warn "We run into a timeout when solving a Submodel with GBCSolver. Note that this implies that the Subproblem run for the time passed by timelimit and did not terminate. "
            param.stats.data["GBCStatus"] = "Timeout_Submodel"
            param.stats.data["Opt_status_override"] = "Timeout_Submodel"
        elseif (e isa NumericalIssueException)
            @error "GBCSolver stopped due to a detected numerical issue: $(e.message)"
            # A status ending in _Numerics is already the final, qualified
            # status. Otherwise retain the solver's eventual base status and
            # let set_optimization_status_stats append the numerical marker.
            if endswith(e.status, "_Numerics")
                param.stats.data["GBCStatus"] = e.status
                param.stats.data["Opt_status_override"] = e.status
            else
                param.stats.data["GBCStatus"] = "Numerics"
                param.stats.data["Opt_status_override"] = "Numerics"
            end
            if get(param.stats.data, "GBCStatus", "") == "Terminate_Numerics"
                _write_numerical_diagnostic(param, e, catch_backtrace())
            end
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
                xsol = Dict(a => value(master.link_vars[a]) for a in master.A)
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
        if termination_status(master.model) == MOI.OPTIMAL || termination_status(master.model) == MOI.LOCALLY_SOLVED || termination_status(master.model) == MOI.TIME_LIMIT || termination_status(master.model) == MOI.INTERRUPTED
            if primal_status(master.model) == MOI.FEASIBLE_POINT
                try
                    mobj = objective_value(master.model)
                    xsol = Dict(a => value(master.link_vars[a]) for a in master.A)
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
        for (idx, s) in enumerate(subs)
            time_s = @elapsed begin
                worker_id = param.solver isa GurobiSolver ?
                    mod1(idx, max(1, get(param.stats.data, "parallel_connector_workers_used", 1))) : 1
                con = build_connectorLP(s, link_vars, subObjs[name(s)], param, timelimit_inner, worker_id)
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
function build_connectorLP(sub::SubSolver, link_vars_master::Dict, subObjvar, parameter::GBCparam, timelimit, worker_id::Integer=1)
    # Generate the ConnectorLP for the passed sub_problem
    @debug "Starting building of ConnectorLP for subpoblem $(name(sub))."

    # Construct LP (no objective or constraints)
    # Note that we use upper bounds to prevent unbounded solutions, see the ConnectorLP implimentation
    optimizer_factory = parameter.solver isa GurobiSolver ?
        (() -> get_worker_optimizer(parameter.solver, worker_id)) :
        (() -> get_next_optimizer(parameter.solver))
    myLP = Model(optimizer_factory)
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
        blc_generator = ConnectorLP_BlC(parameter, sub.A, link_vars_master, sub;
            optimizer_factory=optimizer_factory)
    end

    # build ConnectorLP obj
    return ConnectorLP(myLP, sub.A, link_vars_master, sub, lbm, blc_generator, Vector{ConSubsolCut}(), parameter.g_round_digit, Dict{Symbol,Any}())
end

"""Rebind worker-side models before the first threaded callback invocation."""
function _assign_gurobi_worker_models!(master, subs, clps, solver::GurobiSolver, nworkers::Integer)
    set_optimizer(master.model, () -> Gurobi.Optimizer(solver.env))
    for (idx, con) in enumerate(clps)
        worker_id = mod1(idx, nworkers)
        worker_factory = () -> get_worker_optimizer(solver, worker_id)
        # The ConnectorLP itself is constructed with this factory.  Rebinding
        # here also handles models constructed by an older/custom constructor.
        set_optimizer(con.lp, worker_factory)
        # MIP-based follower subsolvers own a JuMP model that must be rebound
        # to the worker environment. Algorithmic subsolvers such as AStarSolver
        # do not own a mip_model; their connector LP is still rebound above.
        if hasproperty(con.sub_solver, :mip_model)
            set_optimizer(con.sub_solver.mip_model, worker_factory)
        end
        # SubSolverBlCJuMP keeps a separate persistent oracle model created by
        # copy_model. It must use the same worker environment as the main
        # subsolver; otherwise oracle calls can recreate/finalize models in the
        # shared master environment from inside a threaded callback.
        if hasproperty(con.sub_solver, :oracle_model) &&
           !isnothing(con.sub_solver.oracle_model)
            set_optimizer(con.sub_solver.oracle_model, worker_factory)
            con.sub_solver.oracle_optimizer_bound[] = true
        end
        if !isnothing(con.blc_cut_generator)
            set_optimizer(con.blc_cut_generator.lp, worker_factory)
        end
    end
    return nothing
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

function _local_gbc_param(params::GBCparam; runtime=params.runtime)
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
        runtime,
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
    )
end

function _format_parallel_gbc_task_error(subname::AbstractString, err, bt; worker_id=nothing, phase=nothing, cleanup_errors=())
    io = IOBuffer()
    print(io, "Parallel GBC separator failed for subproblem ", subname, ". ")
    !isnothing(worker_id) && print(io, "worker=", worker_id, ". ")
    !isnothing(phase) && print(io, "phase=", phase, ". ")
    showerror(io, err, bt)
    for cleanup_error in cleanup_errors
        print(io, " Cleanup failure: ")
        showerror(io, cleanup_error.error, cleanup_error.backtrace)
    end
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

function gbc_callback_function(cb_data, master::Master, sub_names, clps, subObj, msol_cuts_mapping::Dict, msol_cuts_mapping_blc::Dict, msol_subobj_mapping::Dict, parameter::GBCparam, deadline::Real, worker_locks)
    # x are the linking variables and clps the connectors (one for each sub)
    # subObj are the obj. vars. in master (for each sub)

    # Do not start another separation after the global deadline.  Returning
    # normally lets Gurobi terminate through its own time-limit machinery.
    if deadline - time() <= GBC_CONNECTOR_TIMEOUT_BUFFER
        parameter.stats.data["GBCStatus"] = "Timeout_Submodel"
        parameter.stats.data["Opt_status_override"] = "Timeout_Submodel"
        _request_gurobi_master_termination!(cb_data, "the one-second connector timeout buffer was reached")
        return
    end

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
                # Timeout is cooperative. Tasks already inside separation are
                # allowed to drain; tasks waiting for a worker do not start new
                # Gurobi work once this flag is set.
                abort_separation = Threads.Atomic{Bool}(false)
                timeout_observed = Threads.Atomic{Bool}(false)
                sync_error = nothing
                try
                    @sync for (res_idx, (connector_idx, con, subname, current_subobj)) in enumerate(pending)
                    Threads.@spawn begin
                        worker_id = mod1(connector_idx, length(worker_locks))
                        semaphore_acquired = false
                        worker_locked = false
                        phase = "created"
                        primary_failure = nothing
                        cleanup_failures = NamedTuple{(:error, :backtrace),Tuple{Any,Any}}[]
                        try
                            phase = "acquire_semaphore"
                            Base.acquire(sem)
                            semaphore_acquired = true
                            if abort_separation[]
                                phase = "cancelled_before_lock"
                                task_errors[res_idx] = (
                                    subname=subname,
                                    current_subobj=current_subobj,
                                    err=TimeoutException("GBC separation cancelled after the global time limit was reached."),
                                    bt=nothing,
                                    worker_id=worker_id,
                                    phase=phase,
                                    cleanup_errors=cleanup_failures,
                                )
                                return
                            end
                            phase = "lock_worker"
                            lock(worker_locks[worker_id])
                            worker_locked = true
                            phase = "check_deadline"
                            remaining_time = deadline - time()
                            if remaining_time <= GBC_CONNECTOR_TIMEOUT_BUFFER
                                Threads.atomic_xchg!(abort_separation, true)
                                Threads.atomic_xchg!(timeout_observed, true)
                                phase = "deadline_reached"
                                task_errors[res_idx] = (
                                    subname=subname,
                                    current_subobj=current_subobj,
                                    err=TimeoutException("GBC global time limit reached before separating connector $(subname)."),
                                    bt=nothing,
                                    worker_id=worker_id,
                                    phase=phase,
                                    cleanup_errors=cleanup_failures,
                                )
                                return
                            end
                            phase = "build_local_parameters"
                            local_param = _local_gbc_param(parameter; runtime=remaining_time)
                            result_ref = Ref{Any}(nothing)
                            phase = "genBenders_cut"
                            cuttime = @elapsed begin
                                feas, cut, bigMcut, pobj = genBenders_cut!(con, x_vals, local_param, remaining_time)
                                result_ref[] = (
                                    subname=subname,
                                    current_subobj=current_subobj,
                                    feas=feas,
                                    cut=cut,
                                    bigMcut=bigMcut,
                                    pobj=pobj,
                                    stats=local_param.stats,
                                )
                            end
                            base_result = result_ref[]
                            results[res_idx] = (; base_result..., cuttime=cuttime)
                        catch err
                            primary_failure = (
                                err=err,
                                backtrace=stacktrace(catch_backtrace()),
                            )
                            if err isa TimeoutException
                                Threads.atomic_xchg!(abort_separation, true)
                                Threads.atomic_xchg!(timeout_observed, true)
                            end
                        finally
                            if worker_locked
                                phase = "unlock_worker"
                                try
                                    unlock(worker_locks[worker_id])
                                catch err
                                    push!(cleanup_failures, (error=err, backtrace=catch_backtrace()))
                                end
                            end
                            if semaphore_acquired
                                phase = "release_semaphore"
                                try
                                    Base.release(sem)
                                catch err
                                    push!(cleanup_failures, (error=err, backtrace=catch_backtrace()))
                                end
                            end
                            if !isnothing(primary_failure)
                                task_errors[res_idx] = (
                                    subname=subname,
                                    current_subobj=current_subobj,
                                    err=primary_failure.err,
                                    bt=primary_failure.backtrace,
                                    worker_id=worker_id,
                                    phase=phase,
                                    cleanup_errors=cleanup_failures,
                                )
                            elseif !isempty(cleanup_failures)
                                task_errors[res_idx] = (
                                    subname=subname,
                                    current_subobj=current_subobj,
                                    err=cleanup_failures[1].error,
                                    bt=cleanup_failures[1].backtrace,
                                    worker_id=worker_id,
                                    phase=phase,
                                    cleanup_errors=cleanup_failures,
                                )
                            end
                        end
                    end
                    end
                catch err
                    sync_error = (error=err, backtrace=catch_backtrace())
                    # A cleanup exception can mask the original
                    # TimeoutException.  If the synchronization barrier fails
                    # at the global deadline, preserve the timeout meaning for
                    # status reporting and retain the composite error above as
                    # a diagnostic.
                    if deadline - time() <= GBC_CONNECTOR_TIMEOUT_BUFFER
                        Threads.atomic_xchg!(timeout_observed, true)
                    end
                end

                if !isnothing(sync_error)
                    if timeout_observed[]
                        @error "Parallel GBC separation raised an exception after a subsolver timeout; classifying the run as Timeout_Submodel" exception=(sync_error.error, sync_error.backtrace)
                    else
                        throw(sync_error.error)
                    end
                end

                failures = [failure for failure in task_errors if !isnothing(failure)]
                timeout_failures = [failure for failure in failures if failure.err isa TimeoutException]
                other_failures = [failure for failure in failures if !(failure.err isa TimeoutException)]
                if !isempty(other_failures) && !timeout_observed[]
                    for failure in other_failures
                        @error _format_parallel_gbc_task_error(
                            failure.subname,
                            failure.err,
                            failure.bt;
                            worker_id=get(failure, :worker_id, nothing),
                            phase=get(failure, :phase, nothing),
                            cleanup_errors=get(failure, :cleanup_errors, ()),
                        )
                    end
                    throw(first(other_failures).err)
                elseif !isempty(timeout_failures) || abort_separation[] || timeout_observed[]
                    # Never throw a timeout from inside the Gurobi callback.
                    # Partial cuts are discarded and all spawned tasks have
                    # already drained because this code follows @sync.
                    parameter.stats.data["GBCStatus"] = "Timeout_Submodel"
                    parameter.stats.data["Opt_status_override"] = "Timeout_Submodel"
                    @debug "GBC separation reached its global time limit; connector tasks drained and no further cuts are submitted."
                    _request_gurobi_master_termination!(cb_data, "a connector subsolver reached its time limit")
                    return
                end
            else
                for (res_idx, (_, con, subname, current_subobj)) in enumerate(pending)
                    remaining_time = deadline - time()
                    if remaining_time <= GBC_CONNECTOR_TIMEOUT_BUFFER
                        parameter.stats.data["GBCStatus"] = "Timeout_Submodel"
                        parameter.stats.data["Opt_status_override"] = "Timeout_Submodel"
                        _request_gurobi_master_termination!(cb_data, "the one-second connector timeout buffer was reached in serial separation")
                        return
                    end
                    local_param = parameter.parallel_separation ? _local_gbc_param(parameter; runtime=remaining_time) : parameter
                    result_ref = Ref{Any}(nothing)
                    cuttime = @elapsed begin
                        feas, cut, bigMcut, pobj = genBenders_cut!(con, x_vals, local_param, remaining_time)
                        result_ref[] = (
                            subname=subname,
                            current_subobj=current_subobj,
                            feas=feas,
                            cut=cut,
                            bigMcut=bigMcut,
                            pobj=pobj,
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
