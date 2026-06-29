# HNDP A* Subsolver Tutorial

The HNDP example uses JuBiC's generic [`AStarSolver`](../../subsolvers/astar.md)
wrapper to implement a graph-specific follower oracle. This page explains the
pattern so the HNDP code can serve as a template for custom labeling
subsolvers.

The implementation is in `examples/HNDP/hndp_astar_wrapper.jl`.

## When This Pattern Applies

Use an `AStarSolver`-style subsolver when the follower problem is more naturally
solved by a custom labeling algorithm than by repeatedly building and solving a
JuMP model.

For HNDP, each follower is a shortest-path or constrained-shortest-path problem.
The A* wrapper is therefore a good fit:

- states are partial paths,
- transitions are graph arcs,
- resources track path weight,
- the objective changes depending on the JuBiC oracle call,
- dominance can prune labels at the same node.

## Step 1: Define the Label

The HNDP label stores the current node, the arc used to reach it, accumulated
weight, and a predecessor pointer for path reconstruction:

```julia
struct HNDPAStarLabel
    node::Int
    arc::Union{Tuple{Int,Int},Nothing}
    weight::Float64
    predecessor::Union{HNDPAStarLabel,Nothing}
end
```

For another application, this is the place to store the state information that
the labeling algorithm needs for feasibility, dominance, and solution recovery.

## Step 2: Store Problem Data

HNDP stores graph data and precomputed shortest-path information in
`HNDPAStarStructure`:

```julia
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
```

The important point is that the JuBiC wrapper remains generic. Application data
is stored in the `structure` object and accessed by the methods below.

## Step 3: Implement Start and Goal States

HNDP defines the start and end labels from the user origin and destination:

```julia
function start_state(sol::AStarSolver)
    return HNDPAStarLabel(sol.structure.origin, nothing, 0.0, nothing)
end

function end_state(sol::AStarSolver)
    return HNDPAStarLabel(sol.structure.destination, nothing, -1.0, nothing)
end
```

The goal test compares labels by node:

```julia
function isgoal_w(state::HNDPAStarLabel, goal::HNDPAStarLabel, params::SolverParam)
    return state == goal
end
```

## Step 4: Generate Feasible Neighbor Labels

`neighbours_w` expands a partial path over outgoing graph arcs. HNDP first
checks whether the arc is available under the current first-level solution and
then checks the weight bound and cycle rule:

```julia
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
```

This is the main application-specific part of a labeling subsolver. For a
different problem, replace the expansion and feasibility test by the relevant
state-transition logic.

## Step 5: Provide Costs for Different Oracle Calls

JuBiC calls subsolvers in different cost states. HNDP distinguishes:

- `MASTER_LEVEL`: optimize the first-level contribution `risk`.
- `SUB_PROBLEM_LEVEL`: optimize the true follower cost.
- `CONNECTOR_BASED`: optimize the connector-LP reduced cost used by GBC cuts.

```julia
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
```

The same distinction is used in `prepare(...)` to select or recompute the
heuristic matrix before a solve.

## Step 6: Add Dominance

HNDP labels dominate each other only at the same node. A label is dominated when
it has at least as much weight and at least as much cost as another label:

```julia
function dominate_w(::HNDPAStarLabel, a, b, structure::HNDPAStarStructure, cs::CostStructure, params::SolverParam)
    if a.me.node != b.me.node
        return false
    end

    return a.me.weight >= b.me.weight && a.cost >= b.cost
end
```

Dominance is optional but important for performance. It must be conservative:
discarding a label is only valid when the discarded label cannot lead to a
better feasible solution.

## Step 7: Validate Algorithm Assumptions

The HNDP A* implementation assumes nonnegative active arc costs. This is checked
once before the relevant solve:

```julia
function validate_nonnegative_arc_costs(sol::AStarSolver, xmapping, cs::CostStructure, params::SolverParam)
    # throws ArgumentError if an available arc has negative active cost
end
```

This guard is intentionally explicit. Negative reduced costs can invalidate the
A* heuristic and dominance logic.

## Step 8: Build the JuBiC Subsolver

Finally, HNDP packages the structure into an `AStarSolver`:

```julia
function build_hndp_astar_user(user::User, hndp::HNDPwC, decision_arcs)
    cost_apsp = floyd_warshall_shortest_paths(hndp.mygraph, user.mcost)
    risk_apsp = floyd_warshall_shortest_paths(hndp.mygraph, user.mrisk)

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
```

The returned object can be used by HNDP model builders as
`HNDP_SUBPROBLEM_ASTAR`.

## Implementation Checklist

To implement a custom labeling subsolver through JuBiC, provide:

- a label type,
- a structure type containing application data,
- `start_state` and `end_state`,
- `neighbours_w`,
- `cost_w`,
- `heuristic_w`,
- `isgoal_w`,
- `hashfn_w`,
- `reconstract_path`,
- optional resource helpers such as `used_resource`,
- optional dominance through `dominate_w`,
- optional validation such as `validate_nonnegative_arc_costs`.

The HNDP implementation is deliberately small enough to be used as a reference
for this pattern.
