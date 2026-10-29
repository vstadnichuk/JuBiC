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
