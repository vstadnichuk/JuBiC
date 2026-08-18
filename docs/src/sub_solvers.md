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
- `compute_lower_bound_master_contribution(sub_solver, params, time_limit)`: compute the proven optimal minimum follower contribution to the first-level objective. In `GBC`, this value is used as a global lower bound on the value function in every connector model, so a heuristic incumbent is not sufficient. The built-in `SubSolverJuMP` explicitly disables its separation gap for this initialization solve.
- `solve_sub_for_x(sub_solver, xvals, params, time_limit)`: solve the follower problem for fixed first-level linking values `xvals`. It returns whether a feasible follower solution exists, the follower objective value, the first-level contribution of the returned follower solution, and the follower linking-variable solution.
- `solve_sub_for_x_optimistic(sub_solver, xvals, params, time_limit)`: optional final-incumbent evaluator. It must first minimize the follower objective and then minimize the first-level contribution over all follower-optimal solutions. GBC detects support through method applicability; no capability flag is required.
- `separation!(sub_solver, sval, gvals, kvals, params, time_limit)`: solve the GBC connector separation problem. The method returns a `SubSolution` containing a feasible pricing solution, its objective, and a certified lower bound for the minimization pricing problem.
- `separation_BlC!(sub_solver, sval, kvals, params, time_limit)`: solve the BlC/BlCLag connector separation problem. The returned solution must be bilevel-feasible for the follower problem; otherwise the generated BlC coefficients are not meaningful.
- `supports_bilevel_subproblem_solver(sub_solver)`: return `true` only if the subsolver implements the bilevel separation functionality required by `separation_BlC!`. The default is `false`.
- `separation_exact!`: GBC uses this method for feasibility-cut pricing. Its
  default delegates to `separation!`; a custom subsolver whose ordinary
  separation can stop early must overload it and return equal certified bounds.

There is deliberately no `is_heuristic` capability flag. ConnectorLP infers
whether a pricing call is exact exclusively from the returned bounds. Likewise,
`solve_sub_for_x` must solve the fixed-``x`` follower problem exactly, but it may
return any follower-optimal solution. The optional optimistic method is called
only after GBC terminates. If it is unavailable, JuBiC uses the ordinary cached
response as a valid upper bound and warns that the interval width cannot be
bounded solely by the master gap and connector-pricing errors.

## Certified pricing result

`SubSolution.obj_compare` is the objective ``U`` of the returned feasible
pricing solution. Since connector pricing is a minimization problem, it is an
upper bound. `SubSolution.obj_bound` is a certified lower bound ``L``:

```math
L \le p^* \le U.
```

The five-argument `SubSolution` constructor remains available for exact custom
subsolvers and sets `obj_bound == obj_compare`. A pricing method without a finite
certificate must use `-Inf` for `obj_bound`; JuBiC then reports an infinite
approximation error. `SubSolverJuMP` reads the bound from
`JuMP.objective_bound`. Numerical preprocessing that changes the pricing
objective or feasible region is disabled during certified inexact calls.

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
