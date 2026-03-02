using JuBiC, JuMP                       # For modeling and solving the bilevel problem
using Gurobi                            # Or any other supported solver
using Graphs, SimpleWeightedGraphs      # For Dijkstra's algorithm


# Helper function to extract the shortest path from the Dijkstra result
function _extract_path(distances, source, target)
    # If source and target are the same, the path is empty
    if source == target
        return []
    end

    path = []
    current = target
    last = nothing

    # Backtrack from target to source using parent information
    while current != source
        # Add the edge (current, last) to the path if we have a previous node
        if last !== nothing
            push!(path, (current, last))
        end

        # Move to the parent node
        last = current
        current = distances.parents[current]
    end

    # Add the final edge from source to its child
    push!(path, (source, last))
    return path
end


function _define_upper_level(arcs, opening_costs, sub_names)
    # Create the master model
    master_model = Model(Gurobi.Optimizer)

    # Define binary variables of first level problem
    @variable(master_model, x[arcs], Bin)

    # Define the objective function of the first level problem without the second level variables
    @objective(master_model, Min, sum(opening_costs[a] * x[a] for a in arcs))

    # Create a dictionary for the linking variables
    xdict = Dict(a => x[a] for a in arcs)

    # Create the Master object
    master = Master(master_model, arcs, xdict, sub_names)

    return master
end

struct OwnSubSolver <: JuBiC.SubSolver
    name::String                                        # Name of the sub-problem
    A::Vector{Tuple{Int,Int}}                           # The set of arcs in the graph

    n_nodes::Int                                        # Number of nodes in the graph
    source::Int                                         # The source node
    target::Int                                         # The target node

    transport_costs::Dict{Tuple{Int,Int},Float64}       # The transportation costs for each arc
    transport_time::Dict{Tuple{Int,Int},Float64}        # The transportation time for each arc
end

function JuBiC.name(sub_solver::OwnSubSolver)
    return sub_solver.name
end

function JuBiC.check(sub_solver::OwnSubSolver, params::SolverParam)
    # Check that source and target nodes are valid
    if sub_solver.source > sub_solver.n_nodes || sub_solver.source < 1
        error("The source node $(sub_solver.source) is not in the node set")
    end
    if sub_solver.target > sub_solver.n_nodes || sub_solver.target < 1
        error("The target node $(sub_solver.target) is not in the node set")
    end

    # Check that all arcs have transport costs and transport times defined
    for a in sub_solver.A
        if !haskey(sub_solver.transport_costs, a)
            error("Transport cost for arc $(a) is not defined")
        end
        if !haskey(sub_solver.transport_time, a)
            error("Transport time for arc $(a) is not defined")
        end
    end
end

function JuBiC.capacity_linking(sub_solver::OwnSubSolver, a, params::SolverParam)
    return 1
end

function JuBiC.compute_lower_bound_master_contribution(sub_solver::OwnSubSolver, params::SolverParam, time_limit)
    # Create weighted graph, where the arcs weights are the transport costs
    graph = SimpleWeightedDiGraph(sub_solver.n_nodes)
    for (a, cost) in sub_solver.transport_costs
        u, v = a
        add_edge!(graph, u, v, cost)
    end

    # Compute the lower bound as the shortest path from source to target
    distances = dijkstra_shortest_paths(graph, sub_solver.source)
    lbm = distances.dists[sub_solver.target]

    # Return the computed lower bound
    return lbm
end

function JuBiC.separation!(sub_solver::OwnSubSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)
    # Create weighted graph
    graph = SimpleWeightedDiGraph(sub_solver.n_nodes)
    for a in sub_solver.A
        u, v = a
        cost = sub_solver.transport_time[a] * gvals + sub_solver.transport_costs[a] + kvals[a]
        add_edge!(graph, u, v, cost)
    end

    # Compute the shortest path from source to target
    distances = dijkstra_shortest_paths(graph, sub_solver.source)

    # Get the objective value of the shortest path
    obj_value = distances.dists[sub_solver.target]

    # Extract the shortest path
    shortest_path = _extract_path(distances, sub_solver.source, sub_solver.target)

    # Check if the Benders cut is violated
    is_violated = !(sval < obj_value + 1e-6)

    # Compute risk and cost of the shortest path
    risk = 0.0
    cost = 0.0
    for a in shortest_path
        risk += sub_solver.transport_costs[a]
        cost += sub_solver.transport_time[a]
    end

    return JuBiC.SubSolution(is_violated, risk, cost, shortest_path)
end

function JuBiC.solve_sub_for_x(sub_solver::OwnSubSolver, xvals, params::SolverParam, time_limit)
    # Build the graph including only opened arcs
    graph = SimpleWeightedDiGraph(sub_solver.n_nodes)
    for (a, cost) in sub_solver.transport_time
        if xvals[a] > 0.5
            u, v = a
            add_edge!(graph, u, v, cost)
        end
    end

    # Compute the shortest path from source to target
    distances = dijkstra_shortest_paths(graph, sub_solver.source)
    distance = distances.dists[sub_solver.target]

    if distance == Inf
        # No feasible solution exists (problem is infeasible)
        return false, 0, Dict()
    end

    # Extract the shortest path
    shortest_path = _extract_path(distances, sub_solver.source, sub_solver.target)

    # Create resource mapping: 1.0 if arc is used, 0.0 otherwise
    rmap = Dict(a => (a in shortest_path ? 1.0 : 0.0) for a in sub_solver.A)

    # Return the feasible solution
    return true, distance, rmap
end


function define_instance(arcs, opening_costs, n_nodes, sources, targets, transport_costs, transport_time, sub_names)
    @assert length(sources) == length(targets) == length(sub_names) "The number of sources, targets and subsolver names must be equal"

    # Initialize an empty list to hold the sub-solvers
    sub_solvers = []

    # Create a sub-solver for each source-target pair
    for i in eachindex(sources)
        # Create an OwnSubSolver instance with the provided data
        sub_solver = OwnSubSolver(
            sub_names[i],           # Name of the sub-solver
            arcs,                   # Arcs in the graph
            n_nodes,                # Number of nodes in the graph
            sources[i],             # Source node for this sub-problem
            targets[i],             # Target node for this sub-problem
            transport_costs[i],     # Transportation costs specific to this sub-problem
            transport_time[i]       # Transportation times specific to this sub-problem
        )
        push!(sub_solvers, sub_solver)
    end

    # Define the upper-level (master) problem
    master = _define_upper_level(arcs, opening_costs, sub_names)

    # Create the instance structure combining master and sub-solvers
    instance = Instance(master, sub_solvers)

    return instance
end

n_nodes = 5
arcs = [(1, 2), (1, 3), (2, 4), (3, 4), (2, 5), (3, 5)]

opening_costs = Dict(
    (1, 2) => 10,
    (1, 3) => 15,
    (2, 4) => 5,
    (3, 4) => 8,
    (2, 5) => 12,
    (3, 5) => 9
)

sources = [1, 1]
targets = [4, 5]

transport_costs = [
    Dict((1, 2) => 2, (1, 3) => 3, (2, 4) => 6, (3, 4) => 4, (2, 5) => 7, (3, 5) => 7),
    Dict((1, 2) => 4, (1, 3) => 5, (2, 4) => 8, (3, 4) => 5, (2, 5) => 5, (3, 5) => 9)
]

transport_time = [
    Dict((1, 2) => 5, (1, 3) => 10, (2, 4) => 6, (3, 4) => 10, (2, 5) => 1, (3, 5) => 7),
    Dict((1, 2) => 8, (1, 3) => 12, (2, 4) => 7, (3, 4) => 12, (2, 5) => 1, (3, 5) => 11)
]

sub_names = ["OwnSubSolver1", "OwnSubSolver2"]

instance = define_instance(arcs, opening_costs, n_nodes, sources, targets, transport_costs, transport_time, sub_names)


# Set up the solver parameters
params = GBCparam(GurobiSolver(Gurobi.Env()), true, "./output", "lp", PARETO_OPTIMALITY_ONLY)

mkpath("./output")  # Create the output directory

result = solve_instance!(instance, params)  # Solve the instance
