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
    return structure.shortest_adaptive[state.node, goal.node]
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
    return reduced_cost + get(cs.kvals, (current.node, neighbour.node), 0.0)
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
    cost_apsp = floyd_warshall_shortest_paths(hndp.mygraph, user.mcost)
    risk_apsp = floyd_warshall_shortest_paths(hndp.mygraph, user.mrisk)

    weight_matrix = isnothing(user.weighlimit) ? zeros(Float64, nv(hndp.mygraph), nv(hndp.mygraph)) : Float64.(user.mweight)
    weight_apsp = floyd_warshall_shortest_paths(hndp.mygraph, weight_matrix)
    max_weight = isnothing(user.weighlimit) ? Inf : Float64(user.weighlimit)

    structure = HNDPAStarStructure(
        hndp.mygraph,
        user.origin,
        user.destination,
        max_weight,
        true,
        user.mrisk,
        user.mcost,
        weight_matrix,
        cost_apsp.dists,
        risk_apsp.dists,
        weight_apsp.dists,
        cost_apsp.dists,
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

    apsp = floyd_warshall_shortest_paths(structure.graph, adaptive_costs)
    return apsp.dists
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

            if arc_cost < 0
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
