# Instances for variants of the Hazmat Network Design Problem (HNDP)
using Graphs


struct HNDPwC
    mygraph::Any
    users::Any  # list of User objects
    edgeA::Any  # subset of arcs that forms the resources. Arcs are stored as tuples (NOT as Edge objects)!
    edge_price::Dict  # for each arc the cost of constructing it, i.e., including it into the network

    minweights::Any  # a matrix with the shortest path according to weight for each node pair. If no weight parameter is used, this is 'nothing'. (ASSUMPTIONS: All users have same weight matrix)
    # Optional groups of decision arcs that must be enabled together.  This is
    # used by multimodal instances where one station controls two transfers.
    decision_groups::Any
    decision_budget::Any
end

# Keep the historical five-argument constructor compatible with existing
# instance generators.
HNDPwC(mygraph, users, edgeA, edge_price, minweights) =
    HNDPwC(mygraph, users, edgeA, edge_price, minweights, nothing, nothing)
HNDPwC(mygraph, users, edgeA, edge_price, minweights, decision_groups) =
    HNDPwC(mygraph, users, edgeA, edge_price, minweights, decision_groups, nothing)

struct User
    uname::Any
    origin::Any
    destination::Any

    mrisk::Any  # risk matrix (should support [i, j] access)
    mcost::Any  # cost matrix (should support [i, j] access)
    mweight::Any  # weight matrix (should support [i, j] access)

    weighlimit::Any  # maximum weight a path may have. If no weights are given, it is 'nothing'
end


"""
    hndp_sharing_station_decisions(groups; construction_cost=1.0)

Include both transfer directions in the controllable arc set. Charge each
station once, on its entry arc (the first member of the group). Costs may be
one common nonnegative number or one number per station, in group order.
Coupling equalities and the station budget are added by the model builders.
"""
function hndp_sharing_station_decisions(groups; construction_cost=1.0)
    costs = construction_cost isa Real ? fill(construction_cost, length(groups)) : collect(construction_cost)
    length(costs) == length(groups) || throw(ArgumentError("Provide one construction cost per station."))
    decisions = Tuple{Int,Int}[]
    prices = Dict{Tuple{Int,Int},Float64}()
    for (group, cost) in zip(groups, costs)
        length(group) == 2 || throw(ArgumentError("A sharing station needs entry and return transfer arcs."))
        entry, exit = group
        entry == reverse(exit) || throw(ArgumentError("Station entry and return must be reverse transfer arcs."))
        cost isa Real && isfinite(cost) && cost >= 0 || throw(ArgumentError("Station construction costs must be finite and nonnegative."))
        any(haskey(prices, a) for a in group) && throw(ArgumentError("Station groups must not share transfer arcs."))
        append!(decisions, group)
        prices[entry] = Float64(cost)
        prices[exit] = 0.0
    end
    return decisions, prices
end

# Small in-memory solver regression fixtures; no external topology files.
"""
    build_toy_HNDPwC()

A simple toy instance for the HNDPwC containing 4 nodes, 5 arcs, and 3 paths from 1 to 4. Path 1-2-3-4 is cost optimal but high risk, Path 1-3-4 is infeasible due to weight restrictions, and Path 1-2-4 should be the solution.
"""
function build_toy_HNDPwC()
    # TODO: Small toy example for HNDPwC
    G = DiGraph(4)
    add_edge!(G, 1, 2)
    add_edge!(G, 2, 4)
    add_edge!(G, 1, 3)
    add_edge!(G, 3, 4)
    add_edge!(G, 2, 3)

    # set cost, risk, and weight matrix
    rrisk = [0 20 10 0; 0 0 110 40; 0 0 0 30; 0 0 0 0]
    rcost = [0 3 6 0; 0 0 1 11; 0 0 0 2; 0 0 0 0]
    rweight = [0 3000 9000 0; 0 0 1000 4000; 0 0 0 2000; 0 0 0 0]

    # compute minimal path legths according to rweight matrix
    minweight = floyd_warshall_shortest_paths(G, rweight)

    # generate user
    userslist = []
    origin = 1
    destination = 4
    bound = 10000
    push!(userslist, User("U0", origin, destination, rrisk, rcost, rweight, bound))

    rarcs = [(src(e), dst(e)) for e in edges(G)]
    return HNDPwC(G, userslist, rarcs, Dict((src(e), dst(e)) => 0 for e in edges(G)), minweight)
end

"""
    build_toy_HNDPneg()

A simple toy instance for the HNDP with neg. first level objective cycle. 
"""
function build_toy_HNDPneg()
    # TODO: Small toy example for HNDPwC
    G = DiGraph(5)
    add_edge!(G, 1, 2)
    add_edge!(G, 2, 3)
    add_edge!(G, 3, 4)
    add_edge!(G, 4, 2)
    add_edge!(G, 2, 5)
    add_edge!(G, 1, 5)

    # set cost, risk, and weight matrix
    rrisk = [0 -10 0 0 -1; 0 0 -3 0 -10; 0 0 0 -3 0; 0 -3 0 0 0;0 0 0 0 0]
    rcost = [0 1 0 0 7;0 0 1 0 1;0 0 0 1 0;0 1 0 0 0;0 0 0 0 0]
    rweight = [0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0]

    # compute minimal path legths according to rweight matrix
    minweight = nothing

    # generate user
    userslist = []
    origin = 1
    destination = 5
    bound = nothing
    push!(userslist, User("U0", origin, destination, rrisk, rcost, rweight, bound))

    rarcs = [(1,2), (2,3), (3,4), (4,2), (2,5)]  # all arcs except (1,5) are decision arcs
    return HNDPwC(G, userslist, rarcs, Dict((src(e), dst(e)) => 0 for e in edges(G)), minweight)
end

# A shared user-specific topology keeps unavailable modes out of every model.
# Preserve the historical User layout, including serialized OSM instances.
struct HNDPProfileCosts{M,G} <: AbstractMatrix{Float64}
    values::M
    graph::G
    fixed_distances::Matrix{Float64}
end
Base.size(m::HNDPProfileCosts) = size(m.values)
Base.getindex(m::HNDPProfileCosts, i::Int, j::Int) = Float64(m.values[i,j])

struct HNDPWeightedRisk{M} <: AbstractMatrix{Float64}
    values::M
    demand::Float64
end
Base.size(m::HNDPWeightedRisk) = size(m.values)
Base.getindex(m::HNDPWeightedRisk, i::Int, j::Int) = m.demand * m.values[i,j]

function _hndp_user_instance(user::User, hndp::HNDPwC)
    user.mcost isa HNDPProfileCosts || return hndp
    graph = user.mcost.graph
    graph === hndp.mygraph && return hndp
    return HNDPwC(graph, hndp.users, hndp.edgeA, hndp.edge_price,
                 hndp.minweights, hndp.decision_groups, hndp.decision_budget)
end
