# Implementation of SubSolvers

The `JuBiC` framework enables users to implement various sub-solvers tailored to different second-level problems in bilevel optimization. This documentation outlines how to create a custom sub-solver and describes the two provided sub-solver interfaces: `SubSolverJuMP` and `AStarSolver`.


## Custom SubSolver

A custom sub-solver inherits from the `SubSolver` struct in JuBiC. An example of creating a custom sub-solver can be found in [the third example](./examples/example3.md), where a custom sub-solver is implemented to solve the shortest path problem. The own sub-solver should contain the variable `A`, which holds all identifiers of linking variables, and must provide the following methods:

```@docs
JuBiC.name
JuBiC.check
JuBiC.capacity_linking
JuBiC.compute_lower_bound_master_contribution
JuBiC.separation!
JuBiC.solve_sub_for_x
```


## JuMP-based SubSolver (`SubSolverJuMP`)

This sub-solver interface is designed for solving second-level problems using JuMP models and already provides all the required methods for a sub-solver. The type consists of the following fields:

- `name`: The name (unique identifier) of this sub-problem
- `mip_model`: The JuMP model representing the sub-problem
- `A`: The set of shared resources
- `link_varsC`: A dictionary mapping each shared resource to its corresponding first-level linking variable
- `y_vars`: A dictionary mapping each shared resource to its associated second-level variable
- `link_constraints_capacities`: The capacity parameters in interdiction linking constraints
- `r_objterm`: The objective function term of the first-level problem that depends on the second-level variables
- `c_objterm`: The objective function of the second-level problem
- `extra_cuts`: A user-defined function that adds additional cuts to the master problem based on the solution of the sub-problem, and called after the `mip_model` is solved. It has the signature `timelimit -> (resolve, time)`, where `resolve` indicates whether to re-solve, and `time` is the total time spent in this function.


## A* Search SubSolver (`AStarSolver`)

The A* search sub-solver interface is designed for solving second-level problems modeled as shortest path problems within a state-space graph using an A* search algorithm. This interface requires implementing several wrapper functions that define the problem-specific behavior. The type consists of the following fields:

- `name`: The name (unique identifier) of this sub-problem
- `A`: The set of shared resources
- `structure`: The structure of the state-space graph
- `link_constraints_capacities`: The capacity parameters in interdiction linking constraints
- `max_cost`: The maximum cost allowed; path exceeding this are considered infeasible
- `check_cycles`: Boolean flag indicating whether to check for cycles in the path

In addition to these fields, the A* search sub-solver requires the implementation of several wrapper functions:

```@docs
JuBiC.neighbours_w
JuBiC.heuristic_w
JuBiC.cost_w
JuBiC.isgoal_w
JuBiC.hashfn_w
JuBiC.start_state
JuBiC.end_state
JuBiC.used_resource
JuBiC.risk
JuBiC.cost_subs
JuBiC.prepare
JuBiC.reconstract_path
JuBiC.dominate_w
```
