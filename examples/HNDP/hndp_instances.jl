# Instances for variants of the Hazmat Network Design Problem (HNDP)
using JuBiC


struct HNDPwC
    mygraph::Any
    users::Any  # list of User objects
    edgeA::Any  # subset of arcs that forms the resources. Arcs are stored as tuples (NOT as Edge objects)!
    edge_price::Dict  # for each arc the cost of constructing it, i.e., including it into the network

    minweights::Any  # a matrix with the shortest path according to weight for each node pair. If no weight parameter is used, this is 'nothing'. (ASSUMPTIONS: All users have same weight matrix)
end

struct User
    uname::Any
    origin::Any
    destination::Any

    mrisk::Any  # risk matrix (should support [i, j] access)
    mcost::Any  # cost matrix (should support [i, j] access)
    mweight::Any  # weight matrix (should support [i, j] access)

    weighlimit::Any  # maximum weight a path may have. If no weights are given, it is 'nothing'
end


########## Create Instances of HNDPwC ##########
using Graphs
using CSV
using Random
using GraphPlot

function parse_Sioux_Falls_csv(file_path)
    # Function to parse custom CSV-like input (private function)

    # Read file lines
    lines = readlines(file_path)

    # Find index of "<END OF METADATA>"
    start_idx = findfirst(x -> occursin("Init node", x), lines) + 1

    # Extract relevant data starting after metadata line
    edge_data = lines[start_idx:end]

    # Remove any leading '~' or whitespace characters and split by tab/space
    processed_data = [split(strip(line), r"\s+") for line in edge_data if !isempty(line)]

    return processed_data
end

function create_graph_from_edges_SiouxFalls(edge_list)
    # Function to create directed graph from parsed data (private function)
    g = DiGraph(24)

    for edge in edge_list
        if length(edge) <= 1
            continue  # scip empty lines 
        end
        init_node = parse(Int, edge[1])
        term_node = parse(Int, edge[2])

        add_edge!(g, init_node, term_node)
    end

    return g
end

"""
    build_random_SiouxFalls(nusers, alpha;seed=25, max_cost=100, maxrisk=100, maxweight=100)

Generate a random instance of the 'HNDPwC' with the Sioux Falls network. 
Note that the arc parameters are not symmetric, e.g., cost(i, j) != cost(j,i).
- 'nusers': Number of users (with random edge parameters)
- 'alpha': The alpha in [0;1] parameter for the weight bound on the paths. If 'alpha=0', the weight bound is equal to that of shortest path (according to weight function). 
If 'alpha=1', you can assume the bound to be non-existent.

Additional parameters that can be set are: 
- 'seed': The random seed used to generate arc parameters.
- 'max_cost': cost are randomly generated within 1:'max_cost'
- 'maxrisk': risk are randomly generated within 1:'max_cost'
- 'maxweight': weights are randomly generated within 1:'max_cost'
- 'constructioncost': The cost for including an arc into the network are generated randomly within 0:'constructioncost'
- 'two_stage=false': If true, use same value for risk and cost, ie.e, boths level fully cooperate
"""
function build_random_SiouxFalls(nusers, alpha; seed=25, max_cost=100, maxrisk=100, maxweight=100, constructioncost=0, two_stage=false)
    @assert alpha >= 0 && alpha <= 1
    @assert max_cost >= 0 && maxrisk >= 0 && maxweight >= 0

    # read in graph
    file_path = "examples/data/SF_DNDP_10_base.txt"
    edge_list = parse_Sioux_Falls_csv(file_path)
    graph = create_graph_from_edges_SiouxFalls(edge_list)

    # generate random cost, risk, and weight data
    Random.seed!(seed)
    rcost = rand(1:max_cost, 24, 24)
    rrisk = rand(1:maxrisk, 24, 24)
    rweight = rand(1:maxweight, 24, 24)

    # small hack to model two-stage setting
    if two_stage
        rrisk = rcost
    end

    # set entries in matrices to 0 if no arc for them exists
    for (e1, e2) in Iterators.product(vertices(graph), vertices(graph))
        if !has_edge(graph, e1, e2)
            rcost[e1, e2] = 0
            rrisk[e1, e2] = 0
            rweight[e1, e2] = 0
        end
    end

    # compute minimal path legths according to rweight matrix
    minweight = floyd_warshall_shortest_paths(graph, rweight)

    # create users
    Random.seed!(seed)  # reset seed just to make sure
    userslist = []
    for user = 1:nusers
        origin = rand(1:24)
        # a small hack to generate o != d 
        destination = origin
        while destination == origin
            destination = rand(1:24)
        end

        # compute the maximal weight a path may have for this user
        minweight_user = minweight.dists[origin, destination]
        bound = minweight_user + alpha * (24 * maxweight - minweight_user)

        # create this user
        push!(
            userslist,
            User("U$(user)", origin, destination, rrisk, rcost, rweight, bound),
        )
    end

    # generate construction cost vector
    xconstruction = Dict((src(e), dst(e)) => rand(0:constructioncost) for e in edges(graph))

    # finish 
    rarcs = [(src(e), dst(e)) for e in edges(graph)]
    return HNDPwC(graph, userslist, rarcs, xconstruction, minweight)
end

"""
    build_random_layer_SiouxFalls(nusers, alpha; seed=25, max_cost=100, maxrisk=100, maxweight=100, withweight=false, beta=0.8)

Generate a random instance of the HNDP with the Sioux Falls network consisting of two layers and additional base layer.
The users start and end their trips on a base layer, which does have only connection arcs to other layers.
The transfer arcs connect base and non-base layers. Only the second layer can be modified by the first level.
Note that the arc parameters are not symmetric, e.g., cost(i, j) != cost(j,i).
- 'nusers': Number of users (with random edge parameters)
- 'alpha': The alpha in [0;1] parameter for the weight bound on the paths. If 'alpha=0', the weight bound is equal to that of shortest path (according to weight function). 
If 'alpha=1', you can assume the bound to be non-existent.

Additional parameters that can be set are: 
- 'seed': The random seed used to generate arc parameters.
- 'max_cost': cost are randomly generated within 1:'max_cost'. Cost of connection arcs are 0.
- 'maxrisk': risk are randomly generated within 1:'max_cost'. Note that the risk is set to be negative (profit for first level). 
Further rules: Risk is 0 for all first level arcs, connection arcs to and from first level, and connection arcs down to second level. 
- 'maxweight': weights are randomly generated within 1:'max_cost'. Weights are 0 for connection arcs
- 'constructioncost': The cost for including an arc into the network are generated randomly within 0:'constructioncost'.
- 'withweight': If false, no weight parameters for arcs are generated.
- 'beta': Rescale cost of decision arcs by this factor (value is rounded to integer)
"""
function build_random_layer_SiouxFalls(nusers, alpha; seed=25, max_cost=100, maxrisk=100, maxweight=100, constructioncost=0, withweight=false, beta=0.8)
    @assert alpha >= 0 && alpha <= 1
    @assert max_cost >= 0 && maxrisk >= 0 && maxweight >= 0

    # read in sioux falls graph
    file_path = "examples/data/SF_DNDP_10_base.txt"
    @info "Reading in Sioux Falls graph structure from file $file_path"
    edge_list = parse_Sioux_Falls_csv(file_path)
    sfgraph = create_graph_from_edges_SiouxFalls(edge_list)

    # build sioux falls graphs with base and 2 layer
    @info "Generated single-layer Sioux falls graph strcuture. Starting generating multi-layer construction"
    ## partition nodes
    graph = DiGraph(3*nv(sfgraph))
    nlbase = 1:nv(sfgraph)
    nl1 = (nv(sfgraph)+1):(2*nv(sfgraph))
    nl2 = (2*nv(sfgraph)+1):(3*nv(sfgraph))

    # connection arcs
    for n in nlbase
        # layer 1 connection
        add_edge!(graph, n, n + nv(sfgraph))
        add_edge!(graph, n + nv(sfgraph), n)
        # layer 2 connection
        add_edge!(graph, n, n + 2*nv(sfgraph))
        add_edge!(graph, n + 2*nv(sfgraph), n)
    end

    # arcs for each layer  
    for l in [1, 2]
        for e in edges(sfgraph)
            s = src(e) + l*nv(sfgraph)
            d = dst(e) + l*nv(sfgraph)
            add_edge!(graph, s, d)
        end
    end

    # define set of connect arcs
    rarcs = [(src(e), dst(e)) for e in edges(graph) if src(e) in nl2 && dst(e) in nl2]
    @debug "The decision arcs are $rarcs"


    # generate random cost, risk, and weight data
    @info "Multi-layer graph finished. Now generating parmaeter matrices for the arcs"
    Random.seed!(seed)
    rcost = rand(1:max_cost, nv(graph), nv(graph))
    rrisk = -rand(1:maxrisk, nv(graph), nv(graph))
    rweight = rand(1:maxweight, nv(graph), nv(graph))

    # adjust arc cost to rules 
    for (e1, e2) in Iterators.product(vertices(graph), vertices(graph))
        # set entries in matrices to 0 if no arc for them exists
        if !has_edge(graph, e1, e2)
            rcost[e1, e2] = 0
            rrisk[e1, e2] = 0
            rweight[e1, e2] = 0
        end

        # no profit if no vertex from second layer involved
        if e1 <= 2*nv(sfgraph) && e1 <= 2*nv(sfgraph) 
            rrisk[e1, e2] = 0
        end

        # down arcs from second layer have all parameters 0
        if e1 in nl2 && e1 in nlbase
            rcost[e1, e2] = 0
            rrisk[e1, e2] = 0
            rweight[e1, e2] = 0 
        # up arcs have cost and weigh 0 (free for users)
        elseif (e1 in nl1 || e1 in nl2) && (e2 in nlbase)
            rcost[e1, e2] = 0
            rweight[e1, e2] = 0 
        # no profit for up arcs to layer 1 already ensured
        end

        # rescale decision arcs
        if (e1, e2) in rarcs
            old_val = rcost[e1, e2]
            rcost[e1, e2] = round(beta*rcost[e1, e2])
            @debug "When generating instance scaled down cost of arc $((e1, e2)) from $old_val down to $(rcost[e1, e2])"
        end
    end

    # compute minimal path legths according to rweight matrix if weigh parameters used
    minweight = floyd_warshall_shortest_paths(graph, rweight)
    

    # create users
    @info "Graph generation finished. Now building up $nusers users. Random seed is $seed"
    Random.seed!(seed)  # reset seed just to make sure
    userslist = []
    for user = 1:nusers
        # generate od-pairs only on base layer
        origin = rand(nlbase)
        # a small hack to generate o != d 
        destination = origin
        while destination == origin
            destination = rand(nlbase)
        end

        @debug "Next user number $user has origin=$origin and $destination=$destination"

        # compute the maximal weight a path may have for this user
        bound = -1
        if withweight
            minweight_user = minweight.dists[origin, destination]
            bound = minweight_user + alpha * (24 * maxweight - minweight_user)
        else 
            bound = nothing
        end
        
        @debug "Generated new user with weigh bound $bound. Note that parameter withweight=$withweight"

        # create this user
        push!(
            userslist,
            User("U$(user)", origin, destination, rrisk, rcost, rweight, bound),
        )
    end
    @info "Finished generating users. Therefore, generation of multi-layer Sioux Falls network with $nusers users, alpha=$alpha $seed as random seed complete."


    # if no weight parameter requested, remove them
    if !withweight
        @debug "Removed weight matrices as parameter withweight=$withweight"
        minweight = nothing
        rrisk = nothing
    end

    # generate construction cost vector. As fixed network induces constant cost, we can assign random cost to it 
    xconstruction = Dict((src(e), dst(e)) => rand(0:constructioncost) for e in edges(graph))

    # finish
    return HNDPwC(graph, userslist, rarcs, xconstruction, minweight)
end


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