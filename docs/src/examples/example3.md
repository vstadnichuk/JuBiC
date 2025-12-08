# Example 3: Custom Pricing Solver

This document presents a step-by-step guide on implementing a custom pricing solver within the `JuBiC` framework. 

---

## Problem Definition

Imagine a logistics company that is managing transportation networks which connects various distribution centers to retail stores. While they manage the possible available routes for transportation, individual delivery drivers can choose between these routes based on their own preferences and constraints to transport certain resources. The company aims to optimize its operations by minimizing the overall costs associated with opening specific routes and utilizing the opened routes for transportation of different types of goods.

We model the problem as a two-level optimization problem:
- **First-Level Problem**: Decide which routes should be opened for transporting goods.
- **Second-Level Problems**: For each type of resource, determine the optimal transportation plan given the opened routes from the first-level decision.

The problem is represented on a directed graph $\mathcal{G} = (\mathcal{N}, \mathcal{A})$, where:
- Nodes $\mathcal{N}$ denote distinct regions (e.g., cities or delivery hubs), and
- Arcs $\mathcal{A} \subseteq \mathcal{N} \times \mathcal{N}$ denote the set of directed routes between these regions.

In addition, we have several resources $\mathcal{U}$ to be transported, each from a source node $s$ to a target node $t$. The balance vector $b_v^u$ indicates the supply and demand of resource $u$ at node $v$, which is $1$ at the source node, $-1$ at the target node, and $0$ elsewhere.
Each route $a \in \mathcal{A}$ has an associated opening costs $c_a$, as well as costs $c_a^u$ for utilization and transportation time $\tau_a^u$ with each resource type $u \in \mathcal{U}$.

The first-level decision-maker needs to decide which routes should be opened for transporting goods, indicated by the variable $x_a$ for every $a \in \mathcal{A}$. Conversely, the second-level variable $y_a^u$ indicates whether route $a$ is used for transportation with resource type $u$. The objective of the first-level problem is to minimize the total costs and the problem can be formulated as follows:

```math
\begin{aligned}
\min \quad & \sum_{a \in \mathcal{A}} c_{a} x_{a} + \sum_{u \in \mathcal{U}} \sum_{a \in \mathcal{A}} c_{a}^u y_{a}^u \\
s.t. \quad & x_{a} \in \{0,1\}   \quad & \forall a \in \mathcal{A}
\end{aligned}
```

In addition, for each resource type $u \in \mathcal{U}$, there is a corresponding second-level problem, which aims to find the minimal transportation time for resource $u$ given the opened routes from the first-level decision. That is, the $u$-th second-level problem can be formulated as follows:

```math
\begin{aligned}
\min \quad & \sum_{u \in \mathcal{U}} \sum_{a \in \mathcal{A}} \tau_{a}^u y_{a}^u \\
s.t. \quad & \sum_{\mathclap{\substack{a \in \delta^{+}(v)}}} y_{a}^u - \sum_{\mathclap{\substack{a \in \delta^{-}(v)}}} y_{a}^u = b_{v}^u   \quad & \forall v \in \mathcal{N} \\
    \quad & y_{a}^u \leq x_a     \quad & \forall a \in \mathcal{A} \\
    \quad & y_{a}^u \in \{0,1\}  \quad & \forall a \in \mathcal{A}
\end{aligned}
```

Note that the second-level problems consists of shortest path calculations with additional constraints ensuring that transportation can only occur on open routes. For these second-level problems, we define a custom sub-solver that uses Dijkstra's algorithm to find optimal solutions.

---

## Step-by-Step Implementation

### Step 1: Load Required Packages

Load the necessary packages including Dijkstra's algorithm via the packages `Graphs` and `SimpleWeightedGraphs`.

```julia
using JuBiC, JuMP                       # For modeling and solving the bilevel problem
using Gurobi                            # Or any other supported solver
using Graphs, SimpleWeightedGraphs      # For Dijkstra's algorithm
```

### Step 2: Define Helper Function to Extract Path

While the provided packages can be used for calculating the costs of a shortest path, they do not provide a solution for the path. Therefore, we implement the helper function `_extract_path()` to extract the path from the Dijkstra result.

```julia
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
```

### Step 3: Define the First-Level Problem

Create the master model representing the first-level problem, including binary variables for route openings and the objective function.

```julia
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
```

### Step 4: Define the Custom SubSolver

Define the custom sub-solver inheriting from the `JuBiC.SubSolver`. It should contain the variable `A` as the set of all identifiers of linking variables, and has to provide the following methods:

- `JuBiC.name(sub_solver::OwnSubSolver)`: Return the name of the sub-problem
- `JuBiC.check(sub_solver::OwnSubSolver, params::SolverParam)`: Check if the sub-solver type is correctly defined
- `JuBiC.capacity_linking(sub_solver::OwnSubSolver, a, params::SolverParam)`: Return the capacity of each linking variable
- `JuBiC.compute_lower_bound_master_contribution(sub_solver::OwnSubSolver, params::SolverParam, time_limit)`: Compute a lower bound on the contribution of the sub-solver to the master problem
- `JuBiC.separation!(sub_solver::OwnSubSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)`: Perform the separation for the given solution
- `JuBiC.solve_sub_for_x(sub_solver::OwnSubSolver, xvals, params::SolverParam, time_limit)`: Solve the second-level problem for the given first-level solution

#### Define the Custom SubSolver Struct (`OwnSubSolver`):

This struct has to contain all necessary data for solving the second-level problem. In this case, we need to store the name, the set of arcs, number of nodes, source and target nodes, transportation costs, and transportation times.

```julia
struct OwnSubSolver <: JuBiC.SubSolver
    name::String                                        # Name of the sub-problem
    A::Vector{Tuple{Int,Int}}                           # The set of arcs in the graph

    n_nodes::Int                                        # Number of nodes in the graph
    source::Int                                         # The source node
    target::Int                                         # The target node

    transport_costs::Dict{Tuple{Int,Int},Float64}       # The transportation costs for each arc
    transport_time::Dict{Tuple{Int,Int},Float64}        # The transportation time for each arc
end
```

#### Implement the Required `name()` Method:

As the name is already stored in the struct, we can simply return it.

```julia
function JuBiC.name(sub_solver::OwnSubSolver)
    return sub_solver.name
end
```

#### Implement the Required `check()` Method:

The `check()` method verifies that the sub-solver is correctly defined. That is, it should ensure that all the necessary data is provided and valid. For this example, we check that the source and target nodes are within the valid range and that all arcs have defined transportation costs and times.

```julia
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
```

#### Implement the Required `capacity_linking()` Method:

For this problem, the capacity $C_a^u$ of each linking variable is `1` and can be returned directly.

```julia
function JuBiC.capacity_linking(sub_solver::OwnSubSolver, a, params::SolverParam)
    return 1
end
```

#### Implement the Required `compute_lower_bound_master_contribution()` Method:

This method computes a lower bound on the contribution of the sub-solver to the master problem. In this case, we can compute the shortest path from source to target considering only transportation costs $c_a^u$ using Dijkstra's algorithm.

```julia
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
```

#### Implement the Required `separation!()` Method:

```julia
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
```

#### Implement the Required `solve_sub_for_x()` Method:

Given a first-level solution `xvals`, this method solves the second-level problem by constructing the graph with only the opened arcs and computing the shortest path from source to target. If a feasible path exists, it returns the total transportation time and the resource mapping. Otherwise, it indicates infeasibility.

```julia
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
```

### Step 5: Define an Instance

With the custom sub-solver defined, we can now create a function `define_instance()` that sets up a complete instance of the bilevel problem, including both the first-level and second-level problems.

```julia
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
```

### Step 6: Run a Concrete Example

We have now defined all necessary components to set up and solve a concrete instance, which can be solved using the `solve_instance!()` function from JuBiC.

```julia
# Define the problem data for a concrete example with two resources

# Define the graph structure
n_nodes = 5
arcs = [(1, 2), (1, 3), (2, 4), (3, 4), (2, 5), (3, 5)]

# Define the opening costs for each arc
opening_costs = Dict(
    (1, 2) => 10,
    (1, 3) => 15,
    (2, 4) => 5,
    (3, 4) => 8,
    (2, 5) => 12,
    (3, 5) => 9
)

# Define sources and targets for each resource type
sources = [1, 1]
targets = [4, 5]

# Define transportation costs and times for each resource type and arc
transport_costs = [
    Dict((1, 2) => 2, (1, 3) => 3, (2, 4) => 6, (3, 4) => 4, (2, 5) => 7, (3, 5) => 7),
    Dict((1, 2) => 4, (1, 3) => 5, (2, 4) => 8, (3, 4) => 5, (2, 5) => 5, (3, 5) => 9)
]
transport_time = [
    Dict((1, 2) => 5, (1, 3) => 10, (2, 4) => 6, (3, 4) => 10, (2, 5) => 1, (3, 5) => 7),
    Dict((1, 2) => 8, (1, 3) => 12, (2, 4) => 7, (3, 4) => 12, (2, 5) => 1, (3, 5) => 11)
]

# Define the names of the sub-solvers
sub_names = ["OwnSubSolver1", "OwnSubSolver2"]

# Define the instance
instance = define_instance(arcs, opening_costs, n_nodes, sources, targets, transport_costs, transport_time, sub_names)


# Set up the solver parameters
params = GBCparam(GurobiSolver(Gurobi.Env()), true, "./output", "lp", PARETO_OPTIMALITY_ONLY)

# Create output directory if it doesn't exist
mkpath("./output")

# Solve the bilevel optimization problem
result = solve_instance!(instance, params)
```