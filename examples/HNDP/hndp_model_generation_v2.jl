using JuBiC
using JuMP
using Graphs
using Gurobi
using Base.Threads

include("hndp_network_generation_v2.jl")
include("hndp_astar_wrapper_v2.jl")

const HNDP_BIGM_FIXED_NETWORK_PATH = :fixed_network_path
const HNDP_BIGM_N_MINUS_ONE = :n_minus_one_most_expensive
const HNDP_SUBPROBLEM_MIP = :mip
const HNDP_SUBPROBLEM_ASTAR = :astar
const HNDP_SUBPROBLEM_BLC_JUMP = :blc_jump
const HNDP_HYBRID_FALLBACK_SD = :sd

"""
    _fix_non_decision_arcs!(model, xvars, decision_arcs)

Fix all non-decision arcs to one. In the HNDP setting, only arcs in
`decision_arcs` are controlled by the leader; every other arc belongs to the
fixed network and is therefore always available.
"""
function _fix_non_decision_arcs!(model::JuMP.Model, xvars, decision_arcs)
    for arc in eachindex(xvars)
        if !(arc in decision_arcs)
            @constraint(model, xvars[arc] == 1, base_name = "fix_$(arc)")
        end
    end
end

"""
    build_hndp_blc_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH, subproblem_method=HNDP_SUBPROBLEM_MIP)

Build a JuBiC `BlCMaster` instance for the passed HNDP network. The follower
subproblems can be solved either as MIPs or via the generic JuBiC A* / labeling
framework.
"""
function build_hndp_blc_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
    subproblem_method::Symbol=HNDP_SUBPROBLEM_MIP,
)
    subproblem_method in (HNDP_SUBPROBLEM_MIP, HNDP_SUBPROBLEM_ASTAR) ||
        throw(ArgumentError("Unknown HNDP BlC subproblem method $(subproblem_method)."))

    all_arcs = _hndp_all_arcs(hndp)
    user_names = [string(user.uname) for user in hndp.users]

    hpr = Model(() -> get_next_optimizer(solver))
    @variable(hpr, x[all_arcs], Bin)
    _fix_non_decision_arcs!(hpr, x, hndp.edgeA)

    master_terms = Dict{String,Any}()
    sub_terms = Dict{String,Any}()
    for user in hndp.users
        leader_obj, follower_obj, _ = _add_user_flow_constraints!(
            hpr,
            user,
            hndp,
            all_arcs;
            xvars=x,
            binary_flow=true,
            add_strong_duality=false,
        )
        user_name = string(user.uname)
        opt_l1 = @variable(hpr, base_name = "optL1_$(user_name)")
        opt_l2 = @variable(hpr, base_name = "optL2_$(user_name)")
        @constraint(hpr, opt_l1 == leader_obj, base_name = "leader_obj_$(user_name)")
        @constraint(hpr, opt_l2 == follower_obj, base_name = "follower_obj_$(user_name)")
        master_terms[user_name] = opt_l1
        sub_terms[user_name] = opt_l2
    end

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(hpr, construction_cost_var)
    @constraint(hpr, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(hpr, Min, construction_cost_var + sum(values(master_terms)))

    big_m_by_user = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
    big_m_function(a, uname) = big_m_by_user[string(uname)]
    xdec = Dict(a => x[a] for a in hndp.edgeA)
    master = BlCMaster(hpr, hndp.edgeA, xdec, big_m_function, user_names, sub_terms)

    subs = Any[]
    for user in hndp.users
        push!(subs, _build_hndp_blc_subsolver(user, hndp, all_arcs, solver, subproblem_method))
    end

    return Instance(master, subs)
end

"""
    build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        include_objL2=false,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    )

Build a JuBiC `GBCSolver` instance for the passed HNDP network. The master uses
the same construction-cost objective as the other HNDP builders. When
`partial_decomposition=true`, the user flow constraints are added continuously
to the master so GBC starts from a tighter relaxation.

Supported follower subproblem options:
- `HNDP_SUBPROBLEM_MIP`: classic `SubSolverJuMP`
- `HNDP_SUBPROBLEM_BLC_JUMP`: `SubSolverBlCJuMP`, i.e. a bilevel-aware JuMP
  subsolver that enforces follower optimality inside the subproblem
- `HNDP_SUBPROBLEM_ASTAR`: `AStarSolver`
"""
function build_hndp_gbc_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    partial_decomposition::Bool=true,
    include_objL2::Bool=false,
    subproblem_method::Symbol=HNDP_SUBPROBLEM_MIP,
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
)
    subproblem_method in (HNDP_SUBPROBLEM_MIP, HNDP_SUBPROBLEM_BLC_JUMP, HNDP_SUBPROBLEM_ASTAR) ||
        throw(ArgumentError("Unknown HNDP GBC subproblem method $(subproblem_method)."))

    all_arcs = _hndp_all_arcs(hndp)
    user_names = [string(user.uname) for user in hndp.users]

    mm = Model(() -> get_next_optimizer(solver))
    @variable(mm, x[all_arcs], Bin)
    _fix_non_decision_arcs!(mm, x, hndp.edgeA)

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(mm, construction_cost_var)
    @constraint(mm, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(mm, Min, construction_cost_var)

    xdec_dict = Dict(a => x[a] for a in hndp.edgeA)
    objL2_map = Dict{String,Any}()

    if partial_decomposition
        partial = function (mmodel::JuMP.Model, mobj)
            for user in hndp.users
                leader_obj, follower_obj, _ = _add_user_flow_constraints!(
                    mmodel,
                    user,
                    hndp,
                    all_arcs;
                    xvars=x,
                    binary_flow=false,
                    add_strong_duality=false,
                )
                uname = string(user.uname)
                @constraint(mmodel, leader_obj == mobj[uname], base_name = "objLink_$(uname)")
                if include_objL2
                    opt_l2 = @variable(mmodel, base_name = "optL2_gbc_$(uname)")
                    @constraint(mmodel, opt_l2 == follower_obj, base_name = "objUser_$(uname)")
                    objL2_map[uname] = opt_l2
                end
            end
        end

        master = include_objL2 ?
            Master(mm, hndp.edgeA, xdec_dict, user_names, partial, objL2_map) :
            Master(mm, hndp.edgeA, xdec_dict, user_names, partial)
    else
        include_objL2 &&
            throw(ArgumentError("HNDP GBC can only expose objL2 information when partial_decomposition=true."))
        master = Master(mm, hndp.edgeA, xdec_dict, user_names)
    end

    big_m_by_user = subproblem_method == HNDP_SUBPROBLEM_BLC_JUMP ?
        _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode) :
        Dict{String,Float64}()

    subs = Any[]
    for user in hndp.users
        push!(subs, _build_hndp_gbc_subsolver(user, hndp, all_arcs, solver, subproblem_method, big_m_by_user))
    end

    return Instance(master, subs)
end

function _build_hndp_blc_subsolver(
    user::User,
    hndp::HNDPwC,
    all_arcs,
    solver::SolverWrapper,
    subproblem_method::Symbol,
)
    if subproblem_method == HNDP_SUBPROBLEM_MIP
        sub_model = Model(() -> get_next_optimizer(solver))
        leader_obj, follower_obj, flow = _add_user_flow_constraints!(
            sub_model,
            user,
            hndp,
            all_arcs;
            xvars=nothing,
            binary_flow=true,
            add_strong_duality=false,
        )
        @objective(sub_model, Min, follower_obj)
        set_silent(sub_model)
        y_vars = Dict(a => flow[a] for a in hndp.edgeA)
        return SubSolverJuMP(
            string(user.uname),
            sub_model,
            hndp.edgeA,
            y_vars,
            leader_obj,
            follower_obj,
        )
    elseif subproblem_method == HNDP_SUBPROBLEM_ASTAR
        return build_hndp_astar_user(user, hndp, hndp.edgeA)
    end

    throw(ArgumentError("Unknown HNDP BlC subproblem method $(subproblem_method)."))
end

function _build_hndp_gbc_subsolver(
    user::User,
    hndp::HNDPwC,
    all_arcs,
    solver::SolverWrapper,
    subproblem_method::Symbol,
    big_m_by_user::Dict{String,Float64},
)
    if subproblem_method == HNDP_SUBPROBLEM_MIP
        return _build_hndp_blc_subsolver(user, hndp, all_arcs, solver, HNDP_SUBPROBLEM_MIP)
    elseif subproblem_method == HNDP_SUBPROBLEM_ASTAR
        return _build_hndp_blc_subsolver(user, hndp, all_arcs, solver, HNDP_SUBPROBLEM_ASTAR)
    elseif subproblem_method == HNDP_SUBPROBLEM_BLC_JUMP
        sub_model = Model(() -> get_next_optimizer(solver))
        leader_obj, follower_obj, flow = _add_user_flow_constraints!(
            sub_model,
            user,
            hndp,
            all_arcs;
            xvars=nothing,
            binary_flow=true,
            add_strong_duality=false,
        )
        @objective(sub_model, Min, follower_obj)
        set_silent(sub_model)
        y_vars = Dict(a => flow[a] for a in hndp.edgeA)
        user_big_m = big_m_by_user[string(user.uname)]
        return SubSolverBlCJuMP(
            string(user.uname),
            sub_model,
            hndp.edgeA,
            y_vars,
            leader_obj,
            follower_obj,
            a -> user_big_m,
        )
    end

    throw(ArgumentError("Unknown HNDP GBC subproblem method $(subproblem_method)."))
end

"""
    build_hndp_hybrid_blc_instance(
        hndp,
        solver;
        enumeration_time_limit,
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )

Build a mixed BlC instance for HNDP. Users whose path sets are fully enumerated
before the shared wall-clock deadline are modeled exactly inside the HPR via the
compact path formulation. Users that miss the deadline remain genuine BlC
subproblems and are the only users listed in the `BlCMaster` bookkeeping.

This lets BlC focus its decomposition machinery on the hard users while keeping
the easy users exact and compact in the master.

Returns:
- `instance`: the generated JuBiC BlC instance
- `enum_runtime`: wall-clock time spent in the parallel precomputation stage
- `path_counts`: number of enumerated paths for each path-mode user
- `fallback_users`: names of users that were added as BlC subproblems
"""
function build_hndp_hybrid_blc_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    enumeration_time_limit::Float64,
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
    subproblem_method::Symbol=HNDP_SUBPROBLEM_MIP,
)
    subproblem_method in (HNDP_SUBPROBLEM_MIP, HNDP_SUBPROBLEM_ASTAR) ||
        throw(ArgumentError("Unknown HNDP BlC subproblem method $(subproblem_method)."))

    @info "Starting generation of the hybrid BlC HNDP reformulation."
    all_arcs = _hndp_all_arcs(hndp)
    deadline = time() + enumeration_time_limit

    tasks = Dict{String,Task}()
    for user in hndp.users
        local_user = user
        tasks[string(user.uname)] = Threads.@spawn _compute_user_path_data(local_user, hndp, all_arcs, solver, deadline)
    end

    user_results = Dict{String,Any}()
    for user in hndp.users
        user_results[string(user.uname)] = fetch(tasks[string(user.uname)])
    end

    user_names = [string(user.uname) for user in hndp.users]
    hpr = Model(() -> get_next_optimizer(solver))
    @variable(hpr, x[all_arcs], Bin)
    _fix_non_decision_arcs!(hpr, x, hndp.edgeA)

    master_terms = Dict{String,Any}()
    sub_terms = Dict{String,Any}()
    subs = Any[]
    path_counts = Dict{String,Int}()
    fallback_users = String[]

    for user in hndp.users
        uname = string(user.uname)
        result = user_results[uname]

        if !result.timed_out
            if result.used_fallback
                @warn "No feasible path was found in the fixed-arc network for user $(user.uname). Falling back to the n-1 edges heuristic bound $(result.bound) before path enumeration."
            end

            path_counts[uname] = length(result.paths)
            leader_obj, follower_obj = _add_user_path_constraints!(
                hpr,
                user,
                result.paths,
                x;
                add_strong_duality=true,
            )
            opt_l1 = @variable(hpr, base_name = "optL1_hybrid_blc_path_$(uname)")
            opt_l2 = @variable(hpr, base_name = "optL2_hybrid_blc_path_$(uname)")
            @constraint(hpr, opt_l1 == leader_obj, base_name = "leader_obj_hybrid_blc_path_$(uname)")
            @constraint(hpr, opt_l2 == follower_obj, base_name = "follower_obj_hybrid_blc_path_$(uname)")
            master_terms[uname] = opt_l1
            continue
        end

        push!(fallback_users, uname)
        @warn "Path enumeration reached the global deadline for user $(user.uname). Keeping this user as a BlC subproblem."

        leader_obj, follower_obj, _ = _add_user_flow_constraints!(
            hpr,
            user,
            hndp,
            all_arcs;
            xvars=x,
            binary_flow=true,
            add_strong_duality=false,
        )
        opt_l1 = @variable(hpr, base_name = "optL1_hybrid_blc_fallback_$(uname)")
        opt_l2 = @variable(hpr, base_name = "optL2_hybrid_blc_fallback_$(uname)")
        @constraint(hpr, opt_l1 == leader_obj, base_name = "leader_obj_hybrid_blc_fallback_$(uname)")
        @constraint(hpr, opt_l2 == follower_obj, base_name = "follower_obj_hybrid_blc_fallback_$(uname)")
        master_terms[uname] = opt_l1
        sub_terms[uname] = opt_l2

        push!(subs, _build_hndp_blc_subsolver(user, hndp, all_arcs, solver, subproblem_method))
    end

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(hpr, construction_cost_var)
    @constraint(hpr, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(hpr, Min, construction_cost_var + sum(values(master_terms)))

    big_m_by_user = Dict{String,Float64}()
    if !isempty(fallback_users)
        derived_big_m = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
        for uname in fallback_users
            big_m_by_user[uname] = derived_big_m[uname]
        end
    end
    big_m_function(a, uname) = get(big_m_by_user, string(uname), 0.0)
    xdec = Dict(a => x[a] for a in hndp.edgeA)
    master = BlCMaster(hpr, hndp.edgeA, xdec, big_m_function, fallback_users, sub_terms)

    enum_runtime = min(enumeration_time_limit, maximum((result.runtime for result in values(user_results)); init=0.0))
    @info "Finished generation of the hybrid BlC HNDP reformulation."
    return Instance(master, subs), enum_runtime, path_counts, fallback_users
end

"""
    build_hndp_sd_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH, indicator_constraints=false, bound_duals=true)

Build a strong-duality reformulation of the HNDP. This formulation is only
valid for instances without weight bounds, because then each follower problem is
an LP shortest-path model.
"""
function build_hndp_sd_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
    indicator_constraints::Bool=false,
    bound_duals::Bool=true,
)
    any(!isnothing(user.weighlimit) for user in hndp.users) &&
        throw(ArgumentError("The strong-duality HNDP formulation is currently only supported for instances without weight bounds."))

    all_arcs = _hndp_all_arcs(hndp)
    mip = Model(() -> get_next_optimizer(solver))
    @variable(mip, x[all_arcs], Bin)
    _fix_non_decision_arcs!(mip, x, hndp.edgeA)

    big_m_by_user = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
    big_m_function(a, uname) = big_m_by_user[string(uname)]

    leader_terms = Dict{String,Any}()
    for user in hndp.users
        leader_obj, follower_obj, _ = _add_user_flow_constraints!(
            mip,
            user,
            hndp,
            all_arcs;
            xvars=x,
            binary_flow=false,
            add_strong_duality=true,
            big_m_function=big_m_function,
            indicator_constraints=indicator_constraints,
            bound_duals=bound_duals,
        )
        user_name = string(user.uname)
        opt_l2 = @variable(mip, base_name = "optL2_$(user_name)")
        @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_$(user_name)")
        leader_terms[user_name] = leader_obj
    end

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(mip, construction_cost_var)
    @constraint(mip, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(mip, Min, construction_cost_var + sum(values(leader_terms)))

    return Instance(MIPMaster(mip), nothing)
end

"""
    build_hndp_sd_auto_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH)

Build an HNDP strong-duality reformulation using the experimental generic JuBiC
wrapper around `SubSolverJuMP` LP followers.

This builder is currently experimental. It is useful for prototyping and
cross-checking formulations, but it is not yet considered part of the stable
HNDP model-generation pipeline.
"""
function build_hndp_sd_auto_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
)
    @info "The automated HNDP strong-duality reformulation is currently experimental."
    any(!isnothing(user.weighlimit) for user in hndp.users) &&
        throw(ArgumentError("The automated strong-duality HNDP formulation is currently only supported for instances without weight bounds."))

    all_arcs = _hndp_all_arcs(hndp)
    master_model = Model()
    @variable(master_model, x[all_arcs], Bin)
    _fix_non_decision_arcs!(master_model, x, hndp.edgeA)
    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @objective(master_model, Min, constructioncost)

    subs = Any[]
    for user in hndp.users
        sub_model = Model()
        leader_obj, follower_obj, flow = _add_user_flow_constraints!(
            sub_model,
            user,
            hndp,
            all_arcs;
            xvars=nothing,
            binary_flow=false,
            add_strong_duality=false,
        )
        @objective(sub_model, Min, follower_obj)
        set_silent(sub_model)
        y_vars = Dict(a => flow[a] for a in hndp.edgeA)
        push!(
            subs,
            SubSolverJuMP(
                string(user.uname),
                sub_model,
                hndp.edgeA,
                y_vars,
                leader_obj,
                follower_obj,
            ),
        )
    end

    big_m_values = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
    big_m_function(a, uname) = big_m_values[string(uname)]
    xdec = Dict(a => x[a] for a in hndp.edgeA)
    return build_strong_duality_mip_instance(master_model, xdec, subs, big_m_function)
end

"""
    build_hndp_path_instance(hndp, solver; enumeration_time_limit, parallelize=false)

Wrapper around the sequential and parallel path-model builders. Set
`parallelize=true` to use the parallel precomputation variant.
"""
function build_hndp_path_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    enumeration_time_limit::Float64,
    parallelize::Bool=false,
)
    if parallelize
        return build_hndp_path_instance_parallel(hndp, solver; enumeration_time_limit=enumeration_time_limit)
    end
    return build_hndp_path_instance_sequential(hndp, solver; enumeration_time_limit=enumeration_time_limit)
end

"""
    build_hndp_path_instance_sequential(hndp, solver; enumeration_time_limit)

Build the path-based HNDP reformulation. For each user, JuBiC first computes a
path-length bound from the fixed-arc network. Then it enumerates all feasible
paths in the full network whose follower cost is below that bound and builds a
path-based strong-duality model from this LP formulation.

If no feasible path exists in the fixed-arc network, JuBiC falls back to the
`n - 1` heuristic bound based on the user costs and logs a warning.

The `enumeration_time_limit` is a global budget for all fixed-arc path-bound
computations and all path-enumeration work inside this function. If the budget
is exhausted, JuBiC returns the model built so far and reports
`enum_runtime == enumeration_time_limit` by convention.

Returns:
- `instance`: the generated JuBiC MIP instance
- `enum_runtime`: the cumulative runtime spent on path-bound computation and
  path enumeration
- `path_counts`: number of enumerated paths for each user
"""
function build_hndp_path_instance_sequential(
    hndp::HNDPwC,
    solver::SolverWrapper;
    enumeration_time_limit::Float64,
)
    @info "Starting generation of the path-based HNDP reformulation."
    all_arcs = _hndp_all_arcs(hndp)

    mip = Model(() -> get_next_optimizer(solver))
    @variable(mip, x[all_arcs], Bin)
    _fix_non_decision_arcs!(mip, x, hndp.edgeA)

    total_enum_runtime = 0.0
    path_counts = Dict{String,Int}()
    leader_terms = Dict{String,Any}()
    deadline = time() + enumeration_time_limit
    for user in hndp.users
        if time() >= deadline
            @warn "The global path-enumeration time budget was exhausted before processing user $(user.uname). Returning the partial path-based model built so far."
            total_enum_runtime = enumeration_time_limit
            break
        end

        bound, bound_runtime, used_fallback, timed_out = _fixed_arc_path_bound(user, hndp, all_arcs, solver, deadline)
        total_enum_runtime += bound_runtime
        total_enum_runtime = min(total_enum_runtime, enumeration_time_limit)
        if timed_out
            @warn "The global path-enumeration time budget was exhausted while computing the fixed-arc path bound for user $(user.uname). Returning the partial path-based model built so far."
            total_enum_runtime = enumeration_time_limit
            break
        end
        if used_fallback
            @warn "No feasible path was found in the fixed-arc network for user $(user.uname). Falling back to the n-1 edges heuristic bound $(bound)."
        end

        paths, enum_runtime, timed_out = _enumerate_user_paths_v2(user, hndp, bound, deadline)
        total_enum_runtime += enum_runtime
        total_enum_runtime = min(total_enum_runtime, enumeration_time_limit)
        path_counts[string(user.uname)] = length(paths)

        if isempty(paths)
            @warn "No path was enumerated for user $(user.uname) before the global time budget was exhausted. Returning the partial path-based model built so far."
            total_enum_runtime = enumeration_time_limit
            break
        end

        if timed_out
            @warn "The global path-enumeration time budget was exhausted while enumerating paths for user $(user.uname). The partial set of paths found so far will be used, and no further users will be processed."
            total_enum_runtime = enumeration_time_limit
        end

        leader_obj, follower_obj = _add_user_path_constraints!(
            mip,
            user,
            paths,
            x;
            add_strong_duality=true,
        )
        user_name = string(user.uname)
        opt_l2 = @variable(mip, base_name = "optL2_path_$(user_name)")
        @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_path_$(user_name)")
        leader_terms[user_name] = leader_obj

        timed_out && break
    end

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(mip, construction_cost_var)
    @constraint(mip, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(mip, Min, construction_cost_var + sum(values(leader_terms)))

    @info "Finished generation of the path-based HNDP reformulation."
    return Instance(MIPMaster(mip), nothing), total_enum_runtime, path_counts
end

"""
    build_hndp_path_instance_parallel(hndp, solver; enumeration_time_limit)

Parallel version of the path-model builder. User-specific path-bound
computations and path enumerations are run in parallel, while the final JuMP
model is built serially afterwards. The time budget is interpreted as wall-clock
time for the precomputation stage only.

If any worker times out, the function returns an empty MIP model together with
`enum_runtime == enumeration_time_limit` and an empty path-count dictionary.
"""
function build_hndp_path_instance_parallel(
    hndp::HNDPwC,
    solver::SolverWrapper;
    enumeration_time_limit::Float64,
)
    @info "Starting generation of the parallel path-based HNDP reformulation."
    all_arcs = _hndp_all_arcs(hndp)
    deadline = time() + enumeration_time_limit

    tasks = Dict{String,Task}()
    for user in hndp.users
        local_user = user
        tasks[string(user.uname)] = Threads.@spawn _compute_user_path_data(local_user, hndp, all_arcs, solver, deadline)
    end

    user_results = Dict{String,Any}()
    any_timeout = false
    for user in hndp.users
        result = fetch(tasks[string(user.uname)])
        user_results[string(user.uname)] = result
        any_timeout |= result.timed_out
    end

    if any_timeout
        @warn "Parallel path precomputation hit the global wall-clock time limit. Returning an empty path-based model by convention."
        return _empty_path_instance(solver), enumeration_time_limit, Dict{String,Int}()
    end

    mip = Model(() -> get_next_optimizer(solver))
    @variable(mip, x[all_arcs], Bin)
    _fix_non_decision_arcs!(mip, x, hndp.edgeA)

    leader_terms = Dict{String,Any}()
    path_counts = Dict{String,Int}()
    for user in hndp.users
        result = user_results[string(user.uname)]
        if result.used_fallback
            @warn "No feasible path was found in the fixed-arc network for user $(user.uname). Falling back to the n-1 edges heuristic bound $(result.bound)."
        end

        path_counts[string(user.uname)] = length(result.paths)
        leader_obj, follower_obj = _add_user_path_constraints!(
            mip,
            user,
            result.paths,
            x;
            add_strong_duality=true,
        )
        user_name = string(user.uname)
        opt_l2 = @variable(mip, base_name = "optL2_path_$(user_name)")
        @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_path_$(user_name)")
        leader_terms[user_name] = leader_obj
    end

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(mip, construction_cost_var)
    @constraint(mip, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(mip, Min, construction_cost_var + sum(values(leader_terms)))

    enum_runtime = min(enumeration_time_limit, maximum((result.runtime for result in values(user_results)); init=0.0))
    @info "Finished generation of the parallel path-based HNDP reformulation."
    return Instance(MIPMaster(mip), nothing), enum_runtime, path_counts
end

"""
    build_hndp_hybrid_instance(
        hndp,
        solver;
        enumeration_time_limit,
        fallback_mode=HNDP_HYBRID_FALLBACK_SD,
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
        indicator_constraints=false,
        bound_duals=true,
    )

Build the parallel hybrid HNDP reformulation. Every user gets the same global
wall-clock deadline for the path precomputation stage. Users whose feasible
paths are fully enumerated before the deadline are modeled with the compact
path formulation; users that time out fall back to an alternative formulation.

Returns:
- `instance`: the generated JuBiC MIP instance
- `enum_runtime`: wall-clock time spent in the parallel precomputation stage
- `path_counts`: number of enumerated paths for each path-mode user
- `fallback_users`: names of users that were modeled with the fallback
"""
function build_hndp_hybrid_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    enumeration_time_limit::Float64,
    fallback_mode::Symbol=HNDP_HYBRID_FALLBACK_SD,
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
    indicator_constraints::Bool=false,
    bound_duals::Bool=true,
)
    _validate_hybrid_fallback_mode(fallback_mode)
    @info "Starting generation of the parallel hybrid HNDP reformulation."

    all_arcs = _hndp_all_arcs(hndp)
    deadline = time() + enumeration_time_limit

    tasks = Dict{String,Task}()
    for user in hndp.users
        local_user = user
        tasks[string(user.uname)] = Threads.@spawn _compute_user_path_data(local_user, hndp, all_arcs, solver, deadline)
    end

    user_results = Dict{String,Any}()
    for user in hndp.users
        user_results[string(user.uname)] = fetch(tasks[string(user.uname)])
    end

    mip = Model(() -> get_next_optimizer(solver))
    @variable(mip, x[all_arcs], Bin)
    _fix_non_decision_arcs!(mip, x, hndp.edgeA)

    leader_terms = Dict{String,Any}()
    path_counts = Dict{String,Int}()
    fallback_users = String[]
    big_m_values = nothing

    for user in hndp.users
        uname = string(user.uname)
        result = user_results[uname]

        if !result.timed_out
            if result.used_fallback
                @warn "No feasible path was found in the fixed-arc network for user $(user.uname). Falling back to the n-1 edges heuristic bound $(result.bound) before path enumeration."
            end

            path_counts[uname] = length(result.paths)
            leader_obj, follower_obj = _add_user_path_constraints!(
                mip,
                user,
                result.paths,
                x;
                add_strong_duality=true,
            )
            opt_l2 = @variable(mip, base_name = "optL2_hybrid_path_$(uname)")
            @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_hybrid_path_$(uname)")
            leader_terms[uname] = leader_obj
            continue
        end

        push!(fallback_users, uname)
        @warn "Path enumeration reached the global deadline for user $(user.uname). Using fallback formulation $(fallback_mode) for this user."

        if isnothing(big_m_values)
            big_m_values = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
        end
        big_m_function(a, fallback_uname) = big_m_values[string(fallback_uname)]

        leader_obj, follower_obj = _add_hybrid_fallback_user!(
            mip,
            user,
            hndp,
            all_arcs,
            x,
            solver;
            fallback_mode=fallback_mode,
            big_m_function=big_m_function,
            indicator_constraints=indicator_constraints,
            bound_duals=bound_duals,
        )
        opt_l2 = @variable(mip, base_name = "optL2_hybrid_fallback_$(uname)")
        @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_hybrid_fallback_$(uname)")
        leader_terms[uname] = leader_obj
    end

    constructioncost = sum((hndp.edge_price[a] * x[a] for a in hndp.edgeA); init=0.0)
    @variable(mip, construction_cost_var)
    @constraint(mip, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(mip, Min, construction_cost_var + sum(values(leader_terms)))

    enum_runtime = min(enumeration_time_limit, maximum((result.runtime for result in values(user_results)); init=0.0))
    @info "Finished generation of the parallel hybrid HNDP reformulation."
    return Instance(MIPMaster(mip), nothing), enum_runtime, path_counts, fallback_users
end

function _hndp_all_arcs(hndp::HNDPwC)
    return [(src(e), dst(e)) for e in edges(hndp.mygraph)]
end

function _empty_path_instance(solver::SolverWrapper)
    mip = Model(() -> get_next_optimizer(solver))
    @objective(mip, Min, 0)
    return Instance(MIPMaster(mip), nothing)
end

function _validate_hybrid_fallback_mode(fallback_mode::Symbol)
    fallback_mode == HNDP_HYBRID_FALLBACK_SD ||
        throw(ArgumentError("Unsupported HNDP hybrid fallback mode $(fallback_mode). Currently implemented: $(HNDP_HYBRID_FALLBACK_SD)."))
end

"""
    _add_hybrid_fallback_user!(...)

Dispatch helper for hybrid fallback formulations. The public hybrid builder uses
this helper so additional fallback formulations can be added later without
rewriting the model-assembly logic.
"""
function _add_hybrid_fallback_user!(
    mip::JuMP.Model,
    user::User,
    hndp::HNDPwC,
    all_arcs,
    xvars,
    solver::SolverWrapper;
    fallback_mode::Symbol,
    big_m_function,
    indicator_constraints::Bool,
    bound_duals::Bool,
)
    fallback_mode == HNDP_HYBRID_FALLBACK_SD || _validate_hybrid_fallback_mode(fallback_mode)

    !isnothing(user.weighlimit) &&
        throw(ArgumentError("The HNDP hybrid fallback mode :sd is currently only supported for users without weight bounds. User $(user.uname) would require the fallback."))

    leader_obj, follower_obj, _ = _add_user_flow_constraints!(
        mip,
        user,
        hndp,
        all_arcs;
        xvars=xvars,
        binary_flow=false,
        add_strong_duality=true,
        big_m_function=big_m_function,
        indicator_constraints=indicator_constraints,
        bound_duals=bound_duals,
    )
    return leader_obj, follower_obj
end

"""
    _add_user_flow_constraints!(...)

Add the arc-based follower formulation for one HNDP user. Depending on the
flags, this can be either the binary shortest-path / constrained-shortest-path
subproblem used in decomposition models or its continuous strong-duality
reformulation used in compact MIP models.
"""
function _add_user_flow_constraints!(
    model::JuMP.Model,
    user::User,
    hndp::HNDPwC,
    all_arcs;
    xvars,
    binary_flow::Bool,
    add_strong_duality::Bool,
    big_m_function=nothing,
    indicator_constraints::Bool=false,
    bound_duals::Bool=true,
)
    if binary_flow
        flow = @variable(model, [a in all_arcs], Bin, base_name = "f_$(user.uname)")
    else
        flow = @variable(model, [a in all_arcs], lower_bound = 0, upper_bound = 1, base_name = "f_$(user.uname)")
    end

    if !isnothing(user.weighlimit)
        @constraint(
            model,
            sum(user.mweight[a...] * flow[a] for a in all_arcs) <= user.weighlimit,
            base_name = "weight_$(user.uname)",
        )
    end

    for node in vertices(hndp.mygraph)
        rhs = node == user.origin ? 1 : (node == user.destination ? -1 : 0)
        outgoing = sum((flow[(node, neigh)] for neigh in outneighbors(hndp.mygraph, node)); init=0.0)
        incoming = sum((flow[(neigh, node)] for neigh in inneighbors(hndp.mygraph, node)); init=0.0)
        @constraint(model, outgoing - incoming == rhs, base_name = "flow_$(user.uname)_$(node)")
    end

    if !isnothing(xvars)
        for a in all_arcs
            @constraint(model, flow[a] <= xvars[a], base_name = "link_$(user.uname)_$(a)")
        end
    end

    leader_obj = @expression(model, sum(user.mrisk[a...] * flow[a] for a in all_arcs))
    follower_obj = @expression(model, sum(user.mcost[a...] * flow[a] for a in all_arcs))

    if add_strong_duality
        isnothing(big_m_function) && throw(ArgumentError("A big-M function is required when adding strong-duality constraints."))
        isnothing(xvars) && throw(ArgumentError("Strong-duality constraints require first-level x-variables."))

        potentials = @variable(model, [node in vertices(hndp.mygraph)], lower_bound = 0, base_name = "pi_$(user.uname)")
        if bound_duals
            max_big_m = maximum(big_m_function(a, string(user.uname)) for a in all_arcs)
            for node in vertices(hndp.mygraph)
                set_upper_bound(potentials[node], max_big_m)
            end
        end

        for a in all_arcs
            if indicator_constraints
                @constraint(
                    model,
                    xvars[a] --> {potentials[a[1]] - potentials[a[2]] <= user.mcost[a...]},
                    base_name = "dual_ind_$(user.uname)_$(a)",
                )
            else
                @constraint(
                    model,
                    potentials[a[1]] - potentials[a[2]] <= user.mcost[a...] + big_m_function(a, string(user.uname)) * (1 - xvars[a]),
                    base_name = "dual_$(user.uname)_$(a)",
                )
            end
        end

        @constraint(
            model,
            follower_obj == potentials[user.origin] - potentials[user.destination],
            base_name = "strong_duality_$(user.uname)",
        )
    end

    return leader_obj, follower_obj, flow
end

"""
    _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)

Compute one big-M value per user for the arc-based strong-duality and BlC
models. The current implementation intentionally uses MIP solves for the
fixed-network path bound, because this keeps the logic uniform for constrained
shortest paths with optional weight limits.
"""
function _derive_user_big_m_values(hndp::HNDPwC, all_arcs, solver::SolverWrapper, big_m_mode::Symbol)
    big_m_mode in (HNDP_BIGM_FIXED_NETWORK_PATH, HNDP_BIGM_N_MINUS_ONE) ||
        throw(ArgumentError("Unknown HNDP big-M mode $(big_m_mode)."))

    values = Dict{String,Float64}()
    for user in hndp.users
        uname = string(user.uname)
        if big_m_mode == HNDP_BIGM_FIXED_NETWORK_PATH
            path_cost = _fixed_network_path_cost(user, hndp, all_arcs, solver)
            isnothing(path_cost) &&
                throw(ArgumentError("No feasible path exists in the fixed network for user $(uname), so the fixed-network big-M mode is not available."))
            values[uname] = path_cost + 1.0
        else
            costs = sort([user.mcost[a...] for a in all_arcs], rev=true)
            n_take = min(length(costs), max(nv(hndp.mygraph) - 1, 1))
            values[uname] = sum(costs[1:n_take])
        end
    end
    return values
end

"""
    _fixed_network_path_cost(user, hndp, all_arcs, solver)

Solve the follower problem on the fixed network only, i.e. with all decision
arcs removed. We only use the resulting objective value as a valid path-cost
bound; because costs are assumed to be nonnegative, potential zero-cost cycles
do not invalidate the bound.
"""
function _fixed_network_path_cost(user::User, hndp::HNDPwC, all_arcs, solver::SolverWrapper)
    model = Model(() -> get_next_optimizer(solver))
    @variable(model, flow[a in all_arcs], Bin)

    for a in hndp.edgeA
        @constraint(model, flow[a] == 0, base_name = "forbid_$(user.uname)_$(a)")
    end

    if !isnothing(user.weighlimit)
        @constraint(
            model,
            sum(user.mweight[a...] * flow[a] for a in all_arcs) <= user.weighlimit,
            base_name = "weight_$(user.uname)",
        )
    end

    for node in vertices(hndp.mygraph)
        rhs = node == user.origin ? 1 : (node == user.destination ? -1 : 0)
        outgoing = sum((flow[(node, neigh)] for neigh in outneighbors(hndp.mygraph, node)); init=0.0)
        incoming = sum((flow[(neigh, node)] for neigh in inneighbors(hndp.mygraph, node)); init=0.0)
        @constraint(model, outgoing - incoming == rhs, base_name = "flow_$(user.uname)_$(node)")
    end

    @objective(model, Min, sum(user.mcost[a...] * flow[a] for a in all_arcs))
    set_silent(model)
    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        return objective_value(model)
    elseif status == MOI.INFEASIBLE
        return nothing
    end
    throw(ArgumentError("Unexpected status $(status) while computing fixed-network big-M for user $(user.uname)."))
end

function _n_minus_one_user_bound(user::User, hndp::HNDPwC, all_arcs)
    costs = sort([user.mcost[a...] for a in all_arcs], rev=true)
    n_take = min(length(costs), max(nv(hndp.mygraph) - 1, 1))
    return Float64(sum(costs[1:n_take]))
end

"""
    _fixed_arc_path_bound(user, hndp, all_arcs, solver, deadline)

Compute the path-length threshold used by the path reformulation. JuBiC first
tries to solve the follower problem on the fixed-arc network. If no such path
exists, the function falls back to the `n - 1` heuristic bound and reports that
fallback to the caller. The solve respects the remaining global deadline and
reports whether that budget was exhausted.
"""
function _fixed_arc_path_bound(user::User, hndp::HNDPwC, all_arcs, solver::SolverWrapper, deadline::Float64)
    start_time = time()
    fixed_arcs = Set(a for a in all_arcs if !(a in hndp.edgeA))
    remaining_time = deadline - start_time
    if remaining_time <= 0
        return _n_minus_one_user_bound(user, hndp, all_arcs), 0.0, true, true
    end
    if isempty(fixed_arcs)
        return _n_minus_one_user_bound(user, hndp, all_arcs), time() - start_time, true, false
    end

    model = Model(() -> get_next_optimizer(solver))
    set_time_limit_sec(model, remaining_time)
    set_attribute(model, MOI.NumberOfThreads(), 1)
    @variable(model, flow[a in all_arcs], Bin)

    for a in all_arcs
        if !(a in fixed_arcs)
            @constraint(model, flow[a] == 0, base_name = "forbid_comp_$(user.uname)_$(a)")
        end
    end

    if !isnothing(user.weighlimit)
        @constraint(
            model,
            sum(user.mweight[a...] * flow[a] for a in all_arcs) <= user.weighlimit,
            base_name = "weight_comp_$(user.uname)",
        )
    end

    for node in vertices(hndp.mygraph)
        rhs = node == user.origin ? 1 : (node == user.destination ? -1 : 0)
        outgoing = sum((flow[(node, neigh)] for neigh in outneighbors(hndp.mygraph, node)); init=0.0)
        incoming = sum((flow[(neigh, node)] for neigh in inneighbors(hndp.mygraph, node)); init=0.0)
        @constraint(model, outgoing - incoming == rhs, base_name = "flow_comp_$(user.uname)_$(node)")
    end

    @objective(model, Min, sum(user.mcost[a...] * flow[a] for a in all_arcs))
    set_silent(model)
    optimize!(model)

    status = termination_status(model)
    runtime = time() - start_time
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        return objective_value(model), runtime, false, false
    elseif status == MOI.INFEASIBLE
        return _n_minus_one_user_bound(user, hndp, all_arcs), runtime, true, false
    elseif status == MOI.TIME_LIMIT
        return _n_minus_one_user_bound(user, hndp, all_arcs), runtime, true, true
    end
    throw(ArgumentError("Unexpected status $(status) while computing the fixed-arc path bound for user $(user.uname)."))
end

struct HNDPPath
    arcs::Vector{Tuple{Int,Int}}
    cost::Float64
    risk::Float64
    weight::Float64
end

"""
    _enumerate_user_paths_v2(user, hndp, length_bound, deadline)

Enumerate all feasible simple paths for one user whose cost is at most
`length_bound`. The DFS keeps a `visited` node set, so paths with cycles are
never generated. This is important because the path reformulation should work on
simple follower paths even if the graph itself contains cycles. The search
stops once the shared global deadline is reached and then returns the paths
found so far together with a timeout flag.
"""
function _enumerate_user_paths_v2(user::User, hndp::HNDPwC, length_bound::Float64, deadline::Float64)
    start_time = time()
    mincostpairs = floyd_warshall_shortest_paths(hndp.mygraph, user.mcost)
    paths = HNDPPath[]
    visited = Set{Int}([user.origin])
    timed_out = false

    function dfs(node::Int, arcs::Vector{Tuple{Int,Int}}, cost::Float64, risk::Float64, weight::Float64)
        if time() >= deadline
            timed_out = true
            return
        end

        if node == user.destination
            push!(paths, HNDPPath(copy(arcs), cost, risk, weight))
            return
        end

        for next_node in outneighbors(hndp.mygraph, node)
            # Restrict the enumeration to simple paths by preventing node revisits.
            next_node in visited && continue
            edge = (node, next_node)
            edge_cost = user.mcost[edge...]
            edge_weight = isnothing(user.mweight) ? 0.0 : user.mweight[edge...]
            new_cost = cost + edge_cost
            new_weight = weight + edge_weight

            if !isnothing(user.weighlimit) && new_weight > user.weighlimit
                continue
            end
            if new_cost > length_bound
                continue
            end

            remaining_cost = mincostpairs.dists[next_node, user.destination]
            if isfinite(remaining_cost) && new_cost + remaining_cost > length_bound
                continue
            end

            push!(arcs, edge)
            push!(visited, next_node)
            dfs(
                next_node,
                arcs,
                new_cost,
                risk + user.mrisk[edge...],
                new_weight,
            )
            pop!(arcs)
            delete!(visited, next_node)
            timed_out && return
        end
    end

    dfs(user.origin, Tuple{Int,Int}[], 0.0, 0.0, 0.0)
    if isempty(paths) && !timed_out
        throw(ArgumentError("No feasible path was enumerated for user $(user.uname) within the path bound $(length_bound)."))
    end
    return paths, min(time() - start_time, max(0.0, deadline - start_time)), timed_out
end

function _path_precompute_solver(solver::SolverWrapper)
    return solver
end

function _path_precompute_solver(solver::GurobiSolver)
    # Gurobi environments are not thread-safe, so each worker gets its own
    # environment for the precomputation phase.
    return GurobiSolver()
end

function _compute_user_path_data(user::User, hndp::HNDPwC, all_arcs, solver::SolverWrapper, deadline::Float64)
    local_solver = _path_precompute_solver(solver)
    start_time = time()
    bound, bound_runtime, used_fallback, timed_out = _fixed_arc_path_bound(user, hndp, all_arcs, local_solver, deadline)
    if timed_out
        return (
            bound=bound,
            used_fallback=true,
            paths=HNDPPath[],
            timed_out=true,
            runtime=min(time() - start_time, max(0.0, deadline - start_time)),
        )
    end

    paths, enum_runtime, timed_out = _enumerate_user_paths_v2(user, hndp, bound, deadline)
    total_runtime = min(bound_runtime + enum_runtime, max(0.0, deadline - start_time))
    return (
        bound=bound,
        used_fallback=used_fallback,
        paths=paths,
        timed_out=timed_out,
        runtime=total_runtime,
    )
end

"""
    _add_user_path_constraints!(mip, user, paths, xvars; add_strong_duality)

Add the path-based follower formulation for one user. Path variables form a
convex combination of enumerated feasible paths, while the linking constraints
ensure that a path can only be selected if each of its arcs is available in the
leader solution.
"""
function _add_user_path_constraints!(
    mip::JuMP.Model,
    user::User,
    paths::Vector{HNDPPath},
    xvars;
    add_strong_duality::Bool,
)
    λ = @variable(mip, [p in eachindex(paths)], lower_bound = 0, base_name = "λ_$(user.uname)")
    @constraint(mip, sum(λ) == 1, base_name = "cover_$(user.uname)")

    for (p_idx, path) in enumerate(paths)
        for a in path.arcs
            @constraint(mip, λ[p_idx] <= xvars[a], base_name = "path_link_$(user.uname)_$(p_idx)_$(a)")
        end
    end

    leader_obj = @expression(mip, sum(paths[p].risk * λ[p] for p in eachindex(paths)))
    follower_obj = @expression(mip, sum(paths[p].cost * λ[p] for p in eachindex(paths)))

    if add_strong_duality
        max_path_cost = maximum(path.cost for path in paths)
        π = @variable(mip, base_name = "pi_path_$(user.uname)")
        for (p_idx, path) in enumerate(paths)
            @constraint(
                mip,
                π <= path.cost + sum((max_path_cost - path.cost) * (1 - xvars[a]) for a in path.arcs),
                base_name = "path_dual_$(user.uname)_$(p_idx)",
            )
        end
        @constraint(mip, follower_obj == π, base_name = "path_sd_$(user.uname)")
    end

    return leader_obj, follower_obj
end
