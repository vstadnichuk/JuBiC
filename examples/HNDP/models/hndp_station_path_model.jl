# Station-based path formulation for sharing-station location problems.
#
# This file complements the legacy path-based builders in
# `hndp_model_generation.jl` and leaves them unchanged. It implements the
# formulation of the paper "Bilevel Feasible Paths for the Station Location
# Problem":
#
#   * one binary opening variable per station (decision group), indexed by the
#     anchor arc `first(group)`; no variables for fixed arcs,
#   * station-wise linking  sum_{p : s in S(p)} f_p <= x_s  per user and station,
#   * optimality constraints  theta <= c_q + M_q * sum_{s in S(q)} (1 - x_s)
#     with M_q = max_p c_p - c_q over the enumerated path set,
#   * path enumeration with dominance on the set of required stations S(p)
#     instead of the set of traversed decision arcs.
#
# Requires the definitions of `hndp_model_generation.jl` (HNDPPath,
# `_try_eliminate_fixed_path_user!`, `_fixed_arc_path_bound`, ...).

"""
    HNDPStationData

Map between decision arcs and stations. Every decision group forms one station.
Decision arcs that are not part of any group form a singleton station that does
not consume the opening budget, mirroring the legacy budget semantics.

Fields:
- `arc_station`: station index of each decision arc,
- `stations`: decision arcs of each station; the first arc is the anchor arc
  that indexes the station's opening variable,
- `budgeted`: whether the station counts towards `hndp.decision_budget`.
"""
struct HNDPStationData
    arc_station::Dict{Tuple{Int,Int},Int}
    stations::Vector{Vector{Tuple{Int,Int}}}
    budgeted::Vector{Bool}
end

"""Anchor arc whose opening variable represents station `s`."""
_hndp_station_anchor(data::HNDPStationData, s::Int) = first(data.stations[s])

"""
    hndp_station_data(hndp)

Build the arc-to-station map of an `HNDPwC` instance. Groups must be disjoint,
which holds for all sharing-station instances (one entry and one return arc per
station).
"""
function hndp_station_data(hndp::HNDPwC)
    decision_set = Set(hndp.edgeA)
    arc_station = Dict{Tuple{Int,Int},Int}()
    stations = Vector{Vector{Tuple{Int,Int}}}()
    budgeted = Bool[]
    if !isnothing(hndp.decision_groups)
        for group in hndp.decision_groups
            isempty(group) && continue
            push!(stations, collect(group))
            push!(budgeted, true)
            station = length(stations)
            for arc in group
                arc in decision_set || throw(ArgumentError(
                    "Grouped transfer arc $(arc) is missing from edgeA."))
                haskey(arc_station, arc) && throw(ArgumentError(
                    "Decision arc $(arc) belongs to more than one station; the station-based formulation requires disjoint groups."))
                arc_station[arc] = station
            end
        end
    end
    for arc in hndp.edgeA
        haskey(arc_station, arc) && continue
        push!(stations, [arc])
        push!(budgeted, false)
        arc_station[arc] = length(stations)
    end
    return HNDPStationData(arc_station, stations, budgeted)
end

"""Set of stations required by a path, i.e., S(p) = {s : p contains an arc of s}."""
function hndp_path_stations(path::HNDPPath, data::HNDPStationData)
    required = BitSet()
    for arc in path.arcs
        station = get(data.arc_station, arc, 0)
        station > 0 && push!(required, station)
    end
    return required
end

"""Label of a partial path in the station-dominance enumeration."""
struct HNDPStationLabel
    stations::BitSet
    cost::Float64
    risk::Float64
    weight::Float64
end

"""
    _dominates_station_label(a, b; use_weight)

Label `a` dominates label `b` (both ending at the same node) if `a` requires a
subset of the stations of `b` and is either strictly cheaper or equally
expensive with a leader contribution (`risk`) that is not worse. With an
active resource limit, `a` must also not use more of the resource.
"""
function _dominates_station_label(a::HNDPStationLabel, b::HNDPStationLabel; use_weight::Bool)
    cost_better = a.cost < b.cost
    same_cost_better_risk = a.cost == b.cost && a.risk <= b.risk
    (cost_better || same_cost_better_risk) || return false
    use_weight && !(a.weight <= b.weight) && return false
    issubset(a.stations, b.stations) || return false
    return true
end

"""
    _cleanup_station_dominated_hndp_paths(paths, data)

Remove every complete path for which a strictly cheaper complete path exists
that requires a subset of its stations. Equal-cost alternatives are retained,
and the order of the surviving paths is preserved.
"""
function _cleanup_station_dominated_hndp_paths(paths::Vector{HNDPPath}, data::HNDPStationData)
    length(paths) <= 1 && return paths
    requirements = [hndp_path_stations(path, data) for path in paths]
    order = sortperm(eachindex(paths); by=index -> paths[index].cost)
    frontier = Int[]
    retained = trues(length(paths))
    for index in order
        if any(other -> paths[other].cost < paths[index].cost &&
                        issubset(requirements[other], requirements[index]), frontier)
            retained[index] = false
        else
            push!(frontier, index)
        end
    end
    all(retained) && return paths
    return paths[retained]
end

"""
    _enumerate_user_paths_station_dominance(user, hndp, length_bound, deadline; kwargs...)

Enumerate the simple paths of `user` with cost at most `length_bound`, analogous
to `_enumerate_user_paths_v2`, but with dominance on required stations:

* a partial path is discarded if a stored label at the same node requires a
  subset of its stations and is strictly cheaper, or equally expensive with a
  leader contribution that is not worse (`use_station_dominance=true`),
* a partial path is discarded if it requires more budgeted stations than the
  opening budget (`check_decision_budget=true`),
* after a completed enumeration, complete paths with a strictly cheaper
  complete path requiring a subset of their stations are removed
  (`use_station_dominance=true`).

With `use_station_dominance=false` and `check_decision_budget=false`, all simple
paths within the cost bound are returned. Keyword `station_data` can pass a
precomputed `HNDPStationData` to avoid rebuilding it for every user.

Returns `(paths, runtime, timed_out)`.
"""
function _enumerate_user_paths_station_dominance(
    user::User,
    hndp::HNDPwC,
    length_bound::Float64,
    deadline::Float64;
    use_station_dominance::Bool=true,
    check_decision_budget::Bool=true,
    station_data::Union{Nothing,HNDPStationData}=nothing,
    heuristic_distances=nothing,
    search_observer=nothing,
)
    hndp = _hndp_user_instance(user, hndp)
    data = isnothing(station_data) ? hndp_station_data(hndp) : station_data
    start_time = time()
    effective_bound = length_bound + _hndp_path_bound_epsilon(length_bound)
    remaining = isnothing(heuristic_distances) ?
        _hndp_reverse_dijkstra_distances(hndp.mygraph, user.mcost, user.destination) :
        heuristic_distances
    budget = hndp.decision_budget
    budget_limit = (check_decision_budget && !isnothing(budget)) ? floor(Int, budget) : typemax(Int)
    if check_decision_budget && !isnothing(budget)
        (budget isa Real && isfinite(budget) && budget >= 0) ||
            throw(ArgumentError("The station-opening budget must be finite and nonnegative."))
    end
    use_weight = !isnothing(user.weighlimit)

    paths = HNDPPath[]
    visited = Set{Int}([user.origin])
    # Refs avoid boxed captured variables in the recursive closure.
    timed_out = Ref(false)
    budget_pruned = Ref(false)
    expanded_states = Ref(0)
    last_progress = Ref(start_time)
    labels_by_node = Dict{Int,Vector{HNDPStationLabel}}()

    function prune_or_register!(node::Int, label::HNDPStationLabel)
        use_station_dominance || return false
        bucket = get!(labels_by_node, node, HNDPStationLabel[])
        for existing in bucket
            _dominates_station_label(existing, label; use_weight=use_weight) && return true
        end
        filter!(existing -> !_dominates_station_label(label, existing; use_weight=use_weight), bucket)
        push!(bucket, label)
        return false
    end

    function dfs(node::Int, arcs::Vector{Tuple{Int,Int}}, cost::Float64, risk::Float64,
                 weight::Float64, required::BitSet, budget_used::Int)
        expanded_states[] += 1
        if !isnothing(search_observer) && expanded_states[] % 10000 == 0 && time() - last_progress[] >= 5.0
            search_observer((paths_found=length(paths), expanded_states=expanded_states[], elapsed_seconds=time() - start_time))
            last_progress[] = time()
        end
        if time() >= deadline
            timed_out[] = true
            return
        end
        if node == user.destination
            push!(paths, HNDPPath(copy(arcs), cost, risk, weight))
            return
        end
        for next_node in outneighbors(hndp.mygraph, node)
            next_node in visited && continue
            edge = (node, next_node)
            new_cost = cost + user.mcost[edge...]
            new_weight = weight + (isnothing(user.mweight) ? 0.0 : user.mweight[edge...])
            use_weight && new_weight > user.weighlimit && continue
            new_cost > effective_bound && continue
            to_go = remaining[next_node]
            isfinite(to_go) && new_cost + to_go > effective_bound && continue

            next_required = required
            next_budget_used = budget_used
            station = get(data.arc_station, edge, 0)
            if station > 0 && !(station in required)
                if data.budgeted[station]
                    next_budget_used += 1
                    if next_budget_used > budget_limit
                        budget_pruned[] = true
                        continue
                    end
                end
                next_required = copy(required)
                push!(next_required, station)
            end

            new_risk = risk + user.mrisk[edge...]
            label = HNDPStationLabel(next_required, new_cost, new_risk, new_weight)
            prune_or_register!(next_node, label) && continue

            push!(arcs, edge)
            push!(visited, next_node)
            dfs(next_node, arcs, new_cost, new_risk, new_weight, next_required, next_budget_used)
            pop!(arcs)
            delete!(visited, next_node)
            timed_out[] && return
        end
    end

    dfs(user.origin, Tuple{Int,Int}[], 0.0, 0.0, 0.0, BitSet(), 0)
    # Label dominance does not retract complete paths found earlier; the final
    # cleanup removes the remaining station-dominated complete paths.
    result = (use_station_dominance && !timed_out[]) ?
        _cleanup_station_dominated_hndp_paths(paths, data) : paths
    if isempty(result) && !timed_out[] && !budget_pruned[]
        throw(ArgumentError("No feasible path was enumerated for user $(user.uname) within path bound $(length_bound) (effective bound $(effective_bound)); this indicates an inconsistent network/bound or a path-enumeration bug."))
    end
    return result, min(time() - start_time, max(0.0, deadline - start_time)), timed_out[]
end

"""
    _fixed_network_bound_station_model(user, hndp, solver, deadline; method)

Fixed-network path-cost bound. Uses the precomputed fixed-network distances of
`HNDPProfileCosts` users directly and falls back to `_fixed_arc_path_bound`
otherwise. Returns `(bound, runtime, used_fallback, timed_out)`.
"""
function _fixed_network_bound_station_model(user::User, hndp::HNDPwC, solver::SolverWrapper,
                                            deadline::Float64; method::Symbol=:labeling)
    if method == :labeling && user.mcost isa HNDPProfileCosts && isnothing(user.weighlimit)
        start_time = time()
        time() >= deadline && return Inf, 0.0, true, true
        value = Float64(user.mcost.fixed_distances[user.origin, user.destination])
        isfinite(value) && return value, time() - start_time, false, false
    end
    user_hndp = _hndp_user_instance(user, hndp)
    return _fixed_arc_path_bound(user, user_hndp, _hndp_all_arcs(user_hndp), solver, deadline; method=method)
end

"""Key of the completion bounds of a user: user graph, cost matrix, and destination."""
_station_completion_key(user::User, hndp::HNDPwC) =
    (objectid(_hndp_user_instance(user, hndp).mygraph), objectid(user.mcost), user.destination)

"""
    _station_completion_bounds(hndp)

Reverse-Dijkstra completion bounds for every distinct (user graph, cost matrix,
destination) combination, computed once and in parallel. Users that share these
three objects (e.g., all users of one class with the same destination) share
the bound vector instead of recomputing it in every enumeration.
"""
function _station_completion_bounds(hndp::HNDPwC, users=hndp.users)
    representative = Dict{Tuple{UInt,UInt,Int},User}()
    for user in users
        get!(representative, _station_completion_key(user, hndp), user)
    end
    keys_ = collect(keys(representative))
    values_ = Vector{Any}(undef, length(keys_))
    Threads.@threads :dynamic for i in eachindex(keys_)
        user = representative[keys_[i]]
        graph = _hndp_user_instance(user, hndp).mygraph
        # Store the objects themselves so that lookups can verify identity.
        values_[i] = (graph, user.mcost, user.destination,
                      _hndp_reverse_dijkstra_distances(graph, user.mcost, user.destination))
    end
    return Dict(zip(keys_, values_))
end

"""Completion bounds of `user` from the cache, or `nothing` if not cached (identity-checked)."""
function _station_cached_completion_bounds(cache, user::User, hndp::HNDPwC)
    isnothing(cache) && return nothing
    entry = get(cache, _station_completion_key(user, hndp), nothing)
    isnothing(entry) && return nothing
    graph, costs, destination, distances = entry
    (graph === _hndp_user_instance(user, hndp).mygraph && costs === user.mcost &&
        destination == user.destination) || return nothing
    return distances
end

"""Per-user precomputation (bound and station-dominance enumeration)."""
function _compute_user_station_path_data(
    user::User,
    hndp::HNDPwC,
    solver::SolverWrapper,
    deadline::Float64,
    station_data::HNDPStationData;
    use_station_dominance::Bool=true,
    fixed_path_bound_method::Symbol=:labeling,
    search_observer=nothing,
    completion_bounds=nothing,
    path_enumerator=nothing,
)
    start_time = time()
    bound, bound_runtime, used_fallback, timed_out = _fixed_network_bound_station_model(
        user, hndp, solver, deadline; method=fixed_path_bound_method)
    if timed_out
        return (bound=bound, used_fallback=true, paths=HNDPPath[], timed_out=true,
                bound_runtime=bound_runtime, enumeration_runtime=0.0,
                runtime=min(time() - start_time, max(0.0, deadline - start_time)))
    end
    paths, enum_runtime, timed_out = !isnothing(path_enumerator) ?
        path_enumerator(user, bound, deadline) : _enumerate_user_paths_station_dominance(
        user, hndp, bound, deadline;
        use_station_dominance=use_station_dominance,
        check_decision_budget=true,
        station_data=station_data,
        heuristic_distances=_station_cached_completion_bounds(completion_bounds, user, hndp),
        search_observer=search_observer,
    )
    return (bound=bound, used_fallback=used_fallback, paths=paths, timed_out=timed_out,
            bound_runtime=bound_runtime, enumeration_runtime=enum_runtime,
            runtime=min(bound_runtime + enum_runtime, max(0.0, deadline - start_time)))
end

"""
    _add_user_station_path_constraints!(mip, user, paths, xstation, station_data)

Add the path-based lower-level formulation of one user:

* `sum(λ) == 1`,
* one linking constraint per station s used by any of the paths:
  `sum_{p : s in S(p)} λ_p <= x_s`,
* `θ == sum_p c_p λ_p` and, for every path q,
  `θ <= c_q + M_q * sum_{s in S(q)} (1 - x_s)` with `M_q = max_p c_p - c_q`.

Path variables are continuous: for binary station decisions, the remaining
constraints only admit convex combinations of minimum-cost available paths,
so an optimal vertex selects a single path.
Returns `(leader_obj, follower_obj)`.
"""
function _add_user_station_path_constraints!(
    mip::JuMP.Model,
    user::User,
    paths::Vector{HNDPPath},
    xstation,
    station_data::HNDPStationData,
)
    uname = string(user.uname)
    λ = @variable(mip, [p in eachindex(paths)], lower_bound = 0, base_name = "λ_$(uname)")
    path_variables = get!(mip.ext, :hndp_path_variables, Dict{String,Vector{VariableRef}}())
    path_variables[uname] = collect(λ)
    @constraint(mip, sum(λ) == 1, base_name = "cover_$(uname)")

    required = [hndp_path_stations(path, station_data) for path in paths]
    users_of_station = Dict{Int,Vector{Int}}()
    for (p, stations) in enumerate(required), s in stations
        push!(get!(users_of_station, s, Int[]), p)
    end
    for s in sort!(collect(keys(users_of_station)))
        @constraint(mip, sum(λ[p] for p in users_of_station[s]) <= xstation[s],
                    base_name = "station_link_$(uname)_$(s)")
    end

    leader_obj = @expression(mip, sum(paths[p].risk * λ[p] for p in eachindex(paths)))
    follower_obj = @expression(mip, sum(paths[p].cost * λ[p] for p in eachindex(paths)))
    max_path_cost = maximum(path.cost for path in paths)
    θ = @variable(mip, base_name = "pi_path_$(uname)")
    for (p, path) in enumerate(paths)
        big_m = max_path_cost - path.cost
        if isempty(required[p])
            @constraint(mip, θ <= path.cost, base_name = "path_opt_$(uname)_$(p)")
        else
            @constraint(mip, θ <= path.cost + big_m * sum(1 - xstation[s] for s in required[p]),
                        base_name = "path_opt_$(uname)_$(p)")
        end
    end
    @constraint(mip, follower_obj == θ, base_name = "path_sd_$(uname)")
    return leader_obj, follower_obj
end

"""
    build_hndp_station_path_instance(hndp, solver; enumeration_time_limit, kwargs...)

Station-based counterpart of `build_hndp_path_instance_parallel`. User-specific
bounds and path enumerations run in parallel; the JuMP model is assembled
serially afterwards. Differences to the legacy builder:

* opening variables exist only for stations (indexed by their anchor arc
  `first(group)`, so `mip[:x][first(group)]` remains valid); fixed arcs carry no
  variables and no constraints,
* station-wise linking constraints,
* path enumeration with station dominance (`use_station_dominance=true`).

Users whose only path uses fixed arcs are removed from the MIP and restored in
the solution statistics (`eliminate_single_path_users=true`), exactly as in the
legacy builder. Returns `(instance, enum_runtime, path_counts)`. If any user
times out, an empty instance and an empty count dictionary are returned.

Optional keywords for very large networks:

* `completion_bound_batch=k`: compute the shared completion bounds for at most
  `k` (graph, cost, destination) keys at a time (bounded memory),
* `path_enumerator=f`: replace the labeling by `f(user, bound, deadline)`
  returning `(paths, runtime, timed_out)`, e.g., a problem-specific enumerator
  with the same dominance-reduced result.
"""
function build_hndp_station_path_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    enumeration_time_limit::Float64,
    use_station_dominance::Bool=true,
    fixed_path_bound_method::Symbol=:labeling,
    eliminate_single_path_users::Bool=true,
    enumeration_observer=nothing,
    completion_bound_batch::Union{Nothing,Int}=nothing,
    path_enumerator=nothing,
)
    @info "Starting generation of the station-based path reformulation."
    station_data = hndp_station_data(hndp)
    deadline = time() + enumeration_time_limit
    precompute_start = time()
    worker_count = max(1, min(Threads.nthreads(), length(hndp.users)))
    worker_solvers = [_path_precompute_solver(solver) for _ in 1:worker_count]
    user_results = Vector{Any}(undef, length(hndp.users))
    function enumerate_users!(indices, completion_bounds)
        next_user = Threads.Atomic{Int}(1)
        @sync for worker in 1:worker_count
            worker_solver = worker_solvers[worker]
            Threads.@spawn begin
                while true
                    position = Threads.atomic_add!(next_user, 1)
                    position > length(indices) && break
                    index = indices[position]
                    user = hndp.users[index]
                    isnothing(enumeration_observer) || enumeration_observer(:started, index, user, nothing)
                    user_results[index] = _compute_user_station_path_data(
                        user, hndp, worker_solver, deadline, station_data;
                        use_station_dominance=use_station_dominance,
                        fixed_path_bound_method=fixed_path_bound_method,
                        completion_bounds=completion_bounds,
                        path_enumerator=path_enumerator,
                        search_observer=isnothing(enumeration_observer) ? nothing :
                            data -> enumeration_observer(:progress, index, user, data),
                    )
                    isnothing(enumeration_observer) || enumeration_observer(:finished, index, user, user_results[index])
                end
            end
        end
    end
    if !isnothing(path_enumerator)
        # A problem-specific enumerator `(user, bound, deadline) -> (paths, runtime,
        # timed_out)` replaces the labeling; it needs no completion bounds.
        enumerate_users!(eachindex(hndp.users), nothing)
    elseif isnothing(completion_bound_batch)
        enumerate_users!(eachindex(hndp.users), _station_completion_bounds(hndp))
    else
        # Bounded memory: completion bounds are computed for at most
        # `completion_bound_batch` (graph, cost, destination) keys at a time, the
        # users of these keys are enumerated, and the bounds are released.
        completion_bound_batch > 0 || throw(ArgumentError("completion_bound_batch must be positive"))
        by_key = Dict{Tuple{UInt,UInt,Int},Vector{Int}}()
        for (index, user) in enumerate(hndp.users)
            push!(get!(by_key, _station_completion_key(user, hndp), Int[]), index)
        end
        for chunk in Iterators.partition(sort!(collect(keys(by_key))), completion_bound_batch)
            indices = reduce(vcat, (by_key[key] for key in chunk))
            bounds = _station_completion_bounds(hndp, (hndp.users[first(by_key[key])] for key in chunk))
            enumerate_users!(indices, bounds)
            bounds = nothing
        end
    end
    enum_runtime = min(enumeration_time_limit, time() - precompute_start)
    if any(result.timed_out for result in user_results)
        @warn "Station path precomputation hit the global wall-clock time limit. Returning an empty path-based model by convention."
        return _empty_path_instance(solver), enumeration_time_limit, Dict{String,Int}()
    end

    isnothing(enumeration_observer) || enumeration_observer(:assembly_started, 0, nothing, nothing)
    mip = Model(() -> get_next_optimizer(solver))
    mip.ext[:hndp_single_path_reduction_enabled] = eliminate_single_path_users
    mip.ext[:hndp_station_data] = station_data
    anchors = [_hndp_station_anchor(station_data, s) for s in eachindex(station_data.stations)]
    @variable(mip, x[anchors], Bin)
    xstation = [x[anchor] for anchor in anchors]
    budget = hndp.decision_budget
    if !isnothing(budget)
        budget_terms = [xstation[s] for s in eachindex(xstation) if station_data.budgeted[s]]
        isempty(budget_terms) || @constraint(mip, sum(budget_terms) <= budget, base_name = "decision_budget")
    end

    # The leader objective is accumulated in place. Summing one expression per
    # user with `sum` copies the growing expression at every step (quadratic in
    # the number of users) and dominated the assembly time on large instances.
    leader_objective = AffExpr(0.0)
    path_counts = Dict{String,Int}()
    for (index, user) in enumerate(hndp.users)
        result = user_results[index]
        if result.used_fallback
            @warn "No feasible path was found in the fixed-arc network for user $(user.uname). Falling back to the n-1 edges heuristic bound $(result.bound)."
        end
        uname = string(user.uname)
        path_counts[uname] = length(result.paths)
        if eliminate_single_path_users && _try_eliminate_fixed_path_user!(
            mip, user, result.paths, hndp.edgeA; complete=!result.timed_out)
            add_to_expression!(leader_objective, only(result.paths).risk)
            isnothing(enumeration_observer) || enumeration_observer(:assembled, index, user, nothing)
            user_results[index] = nothing
            continue
        end
        if isempty(result.paths)
            # No admissible path (e.g., pruned by the budget): the MIP is infeasible.
            @constraint(mip, AffExpr(0.0) == 1.0, base_name = "cover_$(uname)")
            path_variables = get!(mip.ext, :hndp_path_variables, Dict{String,Vector{VariableRef}}())
            path_variables[uname] = VariableRef[]
            user_results[index] = nothing
            continue
        end
        leader_obj, follower_obj = _add_user_station_path_constraints!(
            mip, user, result.paths, xstation, station_data)
        opt_l2 = @variable(mip, base_name = "optL2_path_$(uname)")
        @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_path_$(uname)")
        add_to_expression!(leader_objective, leader_obj)
        isnothing(enumeration_observer) || enumeration_observer(:assembled, index, user, nothing)
        user_results[index] = nothing
    end

    construction = sum((sum(hndp.edge_price[a] for a in station_data.stations[s]) * xstation[s]
                        for s in eachindex(xstation)); init=0.0)
    @variable(mip, construction_cost_var)
    @constraint(mip, construction_cost_var == construction, base_name = "construction_cost")
    add_to_expression!(leader_objective, construction_cost_var)
    @objective(mip, Min, leader_objective)
    @info "Finished generation of the station-based path reformulation."
    return Instance(MIPMaster(mip), nothing), enum_runtime, path_counts
end
