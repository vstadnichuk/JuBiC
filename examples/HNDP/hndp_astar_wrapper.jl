using JuBiC
using Graphs

import JuBiC: neighbours_w
import JuBiC: heuristic_w
import JuBiC: cost_w
import JuBiC: isgoal_w
import JuBiC: hashfn_w
import JuBiC: start_state
import JuBiC: end_state
import JuBiC: used_resource
import JuBiC: risk
import JuBiC: cost_subs
import JuBiC: prepare
import JuBiC: reconstract_path
import JuBiC: dominate_w
import JuBiC: validate_nonnegative_arc_costs

const HNDP_ASTAR_NUMERIC_COST_TOLERANCE = 1e-4

function _hndp_clamp_tiny_negative_cost(cost::Real)
    value = Float64(cost)
    if value < 0.0 && value >= -HNDP_ASTAR_NUMERIC_COST_TOLERANCE
        return 0.0
    end
    return value
end

"""
    HNDPAStarLabel

Label used by the generic JuBiC `AStarSolver` wrapper for HNDP follower
problems. The predecessor pointer allows path reconstruction and optional
cycle checks.
"""
struct HNDPAStarLabel
    node::Int
    arc::Union{Tuple{Int,Int},Nothing}
    weight::Float64
    predecessor::Union{HNDPAStarLabel,Nothing}
end

"""
    HNDPAStarStructure

Problem data required by the generic JuBiC A* / labeling routines for one HNDP
user.
"""
mutable struct HNDPAStarStructure
    graph::DiGraph
    origin::Int
    destination::Int
    max_weight::Float64
    check_cycles::Bool
    mrisk
    mcost
    mweight
    shortest_cost
    shortest_risk
    shortest_weight
    shortest_adaptive
end

function Base.:(==)(a::HNDPAStarLabel, b::HNDPAStarLabel)
    return a.node == b.node
end

function Base.hash(label::HNDPAStarLabel)
    return Base.hash(label.node)
end

function neighbours_w(state::HNDPAStarLabel, xsol, structure::HNDPAStarStructure, params::SolverParam)
    neighbours = HNDPAStarLabel[]
    for next_node in outneighbors(structure.graph, state.node)
        arc = (state.node, next_node)
        if !_hndp_edge_allowed(xsol, arc)
            continue
        end

        new_label = HNDPAStarLabel(
            next_node,
            arc,
            state.weight + structure.mweight[arc...],
            state,
        )
        if _hndp_label_feasible(new_label, structure)
            push!(neighbours, new_label)
        end
    end
    return neighbours
end

function heuristic_w(state::HNDPAStarLabel, goal::HNDPAStarLabel, xsol, structure::HNDPAStarStructure, cs::CostStructure, params::SolverParam)
    adaptive = structure.shortest_adaptive
    # The normal weighted/ConnectorLP path uses a node-to-current-sink vector;
    # the compatibility fallback for negative adaptive costs is an all-pairs
    # matrix.
    return ndims(adaptive) == 1 ? adaptive[state.node] : adaptive[state.node, goal.node]
end

function cost_w(current::HNDPAStarLabel, neighbour::HNDPAStarLabel, structure::HNDPAStarStructure, cs::CostStructure, params::SolverParam)
    if cs.cost_state == MASTER_LEVEL
        return structure.mrisk[current.node, neighbour.node]
    elseif cs.cost_state == SUB_PROBLEM_LEVEL
        return structure.mcost[current.node, neighbour.node]
    end

    @assert cs.cost_state == CONNECTOR_BASED
    reduced_cost =
        structure.mcost[current.node, neighbour.node] * cs.gval +
        structure.mrisk[current.node, neighbour.node]
    return _hndp_clamp_tiny_negative_cost(reduced_cost + get(cs.kvals, (current.node, neighbour.node), 0.0))
end

function isgoal_w(state::HNDPAStarLabel, goal::HNDPAStarLabel, params::SolverParam)
    return state == goal
end

function hashfn_w(state::HNDPAStarLabel)
    return hash(state)
end

function start_state(sol::AStarSolver)
    return HNDPAStarLabel(sol.structure.origin, nothing, 0.0, nothing)
end

function end_state(sol::AStarSolver)
    return HNDPAStarLabel(sol.structure.destination, nothing, -1.0, nothing)
end

function used_resource(state::HNDPAStarLabel, A, structure::HNDPAStarStructure)
    if !isnothing(state.arc) && state.arc in A
        return [state.arc]
    end
    return Tuple{Int,Int}[]
end

function risk(current::HNDPAStarLabel, neighbour::HNDPAStarLabel, structure::HNDPAStarStructure)
    return structure.mrisk[current.node, neighbour.node]
end

function cost_subs(current::HNDPAStarLabel, neighbour::HNDPAStarLabel, structure::HNDPAStarStructure)
    return structure.mcost[current.node, neighbour.node]
end

function prepare(structure::HNDPAStarStructure, cs::CostStructure, params::SolverParam)
    if cs.cost_state == MASTER_LEVEL
        structure.shortest_adaptive = structure.shortest_risk
    elseif cs.cost_state == SUB_PROBLEM_LEVEL
        structure.shortest_adaptive = structure.shortest_cost
    else
        @assert cs.cost_state == CONNECTOR_BASED
        structure.shortest_adaptive = _hndp_calculate_shortest_matrix(structure, cs)
    end
end

function reconstract_path(goal::HNDPAStarLabel)
    path = HNDPAStarLabel[]
    current = goal
    while !isnothing(current)
        pushfirst!(path, current)
        current = current.predecessor
    end
    return path
end

function dominate_w(::HNDPAStarLabel, a, b, structure::HNDPAStarStructure, cs::CostStructure, params::SolverParam)
    if a.me.node != b.me.node
        return false
    end

    return a.me.weight >= b.me.weight && a.cost >= b.cost
end

"""
    build_hndp_astar_user(user, hndp, decision_arcs)

Build an `AStarSolver` wrapper for one HNDP follower. The resulting subsolver is
compatible with `BlCSolver` via the standard `solve_sub_for_x` interface.
"""
function build_hndp_astar_user(user::User, hndp::HNDPwC, decision_arcs)
    has_weight_limit = !isnothing(user.weighlimit)
    n = nv(hndp.mygraph)
    weight_matrix = has_weight_limit ? Float64.(user.mweight) : zeros(Float64, n, n)
    # These matrices are only used as admissible heuristics.  For an
    # unweighted user, leave them empty and let the labeling routine run as
    # ordinary Dijkstra (zero heuristic).  In particular, do not run an
    # unnecessary all-pairs shortest-path computation for every user.
    # Keep the static master/subproblem heuristic representation unchanged;
    # the single-sink optimization below is specifically for the adaptive
    # ConnectorLP objective.  For an unweighted user these are never read.
    cost_heuristic = floyd_warshall_shortest_paths(hndp.mygraph, Float64.(user.mcost)).dists
    risk_heuristic = floyd_warshall_shortest_paths(hndp.mygraph, Float64.(user.mrisk)).dists
    # Unlike the objective heuristic, resource feasibility needs the minimum
    # weight between every intermediate node and the destination. Keep this
    # one all-pairs matrix; it is not recomputed during connector separation.
    weight_heuristic = has_weight_limit ? floyd_warshall_shortest_paths(hndp.mygraph, weight_matrix).dists : zeros(Float64, n, n)
    max_weight = has_weight_limit ? Float64(user.weighlimit) : Inf

    structure = HNDPAStarStructure(
        hndp.mygraph,
        user.origin,
        user.destination,
        max_weight,
        true,
        user.mrisk,
        user.mcost,
        weight_matrix,
        cost_heuristic,
        risk_heuristic,
        weight_heuristic,
        cost_heuristic,
    )

    capacities = Dict(a => 1 for a in decision_arcs)
    return AStarSolver(string(user.uname), decision_arcs, structure, capacities, Inf, false)
end

function _hndp_calculate_shortest_matrix(structure::HNDPAStarStructure, cs::CostStructure)
    @assert cs.cost_state == CONNECTOR_BASED
    adaptive_costs = cs.gval .* structure.mcost + structure.mrisk
    for a in keys(cs.kvals)
        adaptive_costs[a...] += cs.kvals[a]
    end

    # Keep the heuristic consistent with cost_w: tiny negative values caused
    # by ConnectorLP bound tolerances are treated as zero.
    for i in 1:size(adaptive_costs, 1), j in 1:size(adaptive_costs, 2)
        adaptive_costs[i, j] = _hndp_clamp_tiny_negative_cost(adaptive_costs[i, j])
    end

    # Dijkstra is valid only after all active transition costs have been
    # screened as nonnegative.  Keep the legacy Floyd-Warshall fallback for a
    # negative adaptive objective: this is not an A* mode, but it avoids
    # turning a caller-side validation/fallback path into an inadmissible
    # zero-heuristic search.
    for edge in edges(structure.graph)
        if adaptive_costs[src(edge), dst(edge)] < 0.0
            return floyd_warshall_shortest_paths(structure.graph, adaptive_costs).dists
        end
    end
    # Unweighted connector searches retain the established all-pairs fallback
    # because connector costs may be signed even when the static arc data are
    # nonnegative. The resource-constrained case uses the single-sink vector.
    isinf(structure.max_weight) &&
        return floyd_warshall_shortest_paths(structure.graph, adaptive_costs).dists
    return _hndp_distances_to_sink(structure.graph, adaptive_costs, structure.destination)
end

function _hndp_matrix_nonnegative_on_graph(graph::DiGraph, matrix::AbstractMatrix)
    return all(matrix[src(edge), dst(edge)] >= 0.0 for edge in edges(graph))
end

"""Return shortest distances from every node to `sink` in a directed graph."""
function _hndp_distances_to_sink(graph::DiGraph, costs::AbstractMatrix, sink::Int)
    n = nv(graph)
    distances = fill(Inf, n)
    settled = falses(n)
    distances[sink] = 0.0

    # This is the reverse-graph equivalent of a multi-source/single-sink
    # Dijkstra call.  The small O(n²) implementation avoids an extra graph
    # allocation and is adequate for the HNDP graph sizes.
    for _ in 1:n
        current = 0
        current_distance = Inf
        for node in 1:n
            if !settled[node] && distances[node] < current_distance
                current = node
                current_distance = distances[node]
            end
        end
        current == 0 && break
        settled[current] = true

        for predecessor in inneighbors(graph, current)
            settled[predecessor] && continue
            candidate = current_distance + costs[predecessor, current]
            if candidate < distances[predecessor]
                distances[predecessor] = candidate
            end
        end
    end
    return distances
end

function validate_nonnegative_arc_costs(sol::AStarSolver, xmapping, cs::CostStructure, params::SolverParam)
    structure = sol.structure
    graph = structure.graph

    for u in vertices(graph)
        for v in outneighbors(graph, u)
            arc = (u, v)
            if !_hndp_edge_allowed(xmapping, arc)
                continue
            end

            arc_cost = if cs.cost_state == MASTER_LEVEL
                structure.mrisk[arc...]
            elseif cs.cost_state == SUB_PROBLEM_LEVEL
                structure.mcost[arc...]
            else
                @assert cs.cost_state == CONNECTOR_BASED
                structure.mcost[arc...] * cs.gval + structure.mrisk[arc...] + get(cs.kvals, arc, 0.0)
            end

            if arc_cost < 0 && arc_cost >= -HNDP_ASTAR_NUMERIC_COST_TOLERANCE
                @warn "The A*-based HNDP subsolver encountered a numerically insignificant negative connector arc cost; treating it as zero. Subproblem=$(sol.name), arc=$(arc), raw_cost=$(arc_cost), tolerance=$(HNDP_ASTAR_NUMERIC_COST_TOLERANCE)."
                continue
            elseif arc_cost < 0
                throw(ArgumentError(
                    "The A*-based HNDP subsolver does not support negative arc costs for the active objective. " *
                    "Found arc $(arc) with cost $(arc_cost) in cost state $(cs.cost_state) for subproblem $(sol.name).",
                ))
            end
        end
    end
    return nothing
end

function _hndp_edge_allowed(xsol, arc::Tuple{Int,Int})
    return get(xsol, arc, 1) > 0.5
end

function _hndp_has_node(label::HNDPAStarLabel, node::Int)
    current = label
    while !isnothing(current)
        if current.node == node
            return true
        end
        current = current.predecessor
    end
    return false
end

function _hndp_label_feasible(label::HNDPAStarLabel, structure::HNDPAStarStructure)
    min_remaining_weight = structure.shortest_weight[label.node, structure.destination]
    if label.weight + min_remaining_weight > structure.max_weight
        return false
    end

    if structure.check_cycles && !isnothing(label.predecessor) && _hndp_has_node(label.predecessor, label.node)
        return false
    end

    return true
end
