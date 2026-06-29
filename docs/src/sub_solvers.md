# SubSolver Interface

This page documents the technical subsolver interface used by JuBiC. The concrete built-in subsolver wrappers are documented on their own pages:

- [`SubSolverJuMP`](subsolvers/subsolver_jump.md)
- [`SubSolverBlCJuMP`](subsolvers/subsolver_blc_jump.md)
- [`SubSolverMiBS`](subsolvers/subsolver_mibs.md)
- [`AStarSolver`](subsolvers/astar.md)

## Abstract Interface

All subsolvers inherit from `SubSolver`.

JuBiC expects a subsolver object to expose the shared resource set as field `A`
and to implement the functions below. Built-in wrappers implement these methods
for JuMP models, MiBS models, and A*-style labeling oracles.

## Required Interface Functions

Custom subsolvers should implement:

- `name(sub_solver)`: return the unique follower name. This name is used to match subsolvers with `sub_names` stored in the master wrapper and to label statistics and diagnostic output.
- `check(sub_solver, params)`: validate the subsolver before optimization starts. This should catch structural modeling errors early, for example missing linking variables, inconsistent resource sets, or unsupported objective/constraint types.
- `capacity_linking(sub_solver, a, params)`: return the capacity coefficient `C_a` of the linking constraint for resource `a`, i.e. the coefficient in `y_a <= C_a x_a`. Most current JuBiC examples use unit capacities.
- `compute_lower_bound_master_contribution(sub_solver, params, time_limit)`: compute a valid lower bound on the follower contribution to the first-level objective. In `GBC`, this is a lower bound on the value function approximated by the generated cuts. Connector models use this value to build safe cuts and numerical bounds.
- `solve_sub_for_x(sub_solver, xvals, params, time_limit)`: solve the follower problem for fixed first-level linking values `xvals`. It returns whether a feasible follower solution exists, the follower objective value, the first-level contribution of the returned follower solution, and the follower linking-variable solution.
- `separation!(sub_solver, sval, gvals, kvals, params, time_limit)`: solve the GBC connector separation problem. The inputs define the current connector objective and resource prices; the method returns a `SubSolution` describing whether a violated connector constraint was found and the corresponding objective values/resource set.
- `separation_BlC!(sub_solver, sval, kvals, params, time_limit)`: solve the BlC/BlCLag connector separation problem. The returned solution must be bilevel-feasible for the follower problem; otherwise the generated BlC coefficients are not meaningful.
- `supports_bilevel_subproblem_solver(sub_solver)`: return `true` only if the subsolver implements the bilevel separation functionality required by `separation_BlC!`. The default is `false`.

Subsolvers that support thread control can also implement `set_nthreads` and
`set_singlethread`. Solvers without internal parallelism may leave these as
no-ops.

## Exceptions Used in Subsolver Execution

The abstract subsolver layer defines custom exception types that solver drivers
catch and translate into run statuses:

- `TimeoutException`: thrown when a subsolver, connector LP, or helper oracle reaches the provided `time_limit`. This does not mean the model is infeasible; it means JuBiC could not complete the required oracle call within the remaining time budget. Solver drivers usually translate it into a timeout or terminate status.
- `NumericalIssueException`: thrown when connector-cut generation detects a numerical inconsistency that makes the generated cut unsafe, for example a coefficient that should be nonnegative but is materially negative.
- `MibSFailureException`: thrown by MiBS-based wrappers when the external MiBS call fails or does not return the required optimal solution for the requested oracle.

Other ordinary Julia exceptions, such as `ArgumentError`, are used for
invalid input or unsupported model structures and usually indicate that the
instance or wrapper implementation needs to be corrected.

## Subsolver Capability Requirement of `BlCLag`

`BlCLag` requires a subsolver for which:

- `supports_bilevel_subproblem_solver(sub) == true`

because it calls `separation_BlC!` and the associated bilevel connector logic.
