# `SubSolverJuMP`

`SubSolverJuMP` lets a follower problem be passed to JuBiC as a JuMP model. It
is intended for follower problems that can be solved directly as a single-level
JuMP model for different first-level decisions.

When the wrapper is constructed, JuBiC adds internal binary copies of the
first-level linking variables to the JuMP model and adds linking constraints of
the form `y[a] <= x_copy[a]`. Solver routines then change the objective or fix
these internal copies depending on whether the subsolver is used for fixed-`x`
evaluation, GBC separation, or BlC separation.

## Constructor Parameters

The main constructor is:

```julia
SubSolverJuMP(
    name,
    mip_model,
    A,
    y_vars,
    r_objterm,
    c_objterm,
    extra_cuts;
    heuristic=false,
    mip_gap=nothing,
)
```

The shorter form without `extra_cuts` uses a default callback that does not add
additional cuts.

The parameters are:

- `name`: unique name of the follower problem. It must match the corresponding entry in the master's `sub_names`.
- `mip_model`: JuMP model of the follower problem. The model should be a minimization model. JuBiC modifies this model during the solve by changing objectives, adding temporary fixing constraints, and adding internal linking-copy variables.
- `A`: shared resource set used for the linking variables. The keys of `y_vars` must match this set.
- `y_vars`: follower-side variables appearing in the linking constraints. For every `a in A`, JuBiC creates an internal binary copy `x_copy[a]` and adds `y_vars[a] <= x_copy[a]`.
- `r_objterm`: expression for the follower solution's contribution to the first-level objective. This is the quantity returned to the master as the follower contribution.
- `c_objterm`: expression for the follower objective. This is the objective used when solving the follower problem for fixed first-level values.
- `extra_cuts`: optional callback used for branch-and-check style strengthening inside the subsolver MIP. It receives a time limit and returns `(need_resolve, time_spent)`. If `need_resolve` is `true`, JuBiC resolves the follower MIP after the callback added cuts.
- `heuristic`: permit ordinary `GBC` connector pricing to stop with a feasible
  but suboptimal solution. This is an implementation setting of this wrapper,
  not a subsolver capability advertised to the master.
- `mip_gap`: optional relative MIP gap used only for heuristic connector
  pricing. JuBiC reads the actual feasible objective and
  `JuMP.objective_bound`, so the cut certificate is based on the realized
  absolute pricing interval rather than on this requested relative tolerance.

Deprecated constructors with explicit linking-copy variables still exist for
backward compatibility, but JuBiC now ignores the passed copies and creates its
own internal copies with unit linking capacities.

## Certified Inexact Pricing

For a minimization pricing problem, `SubSolverJuMP` returns an incumbent ``U``
and solver bound ``L``. `ConnectorLP` adds a violated row when ``U<s``;
otherwise it records ``\delta=\max(0,s-L)``. See the [GBC approximation-mode
documentation](../solvers/gbc.md#Certified-Inexact-Connector-Pricing) for the
safe-underestimation and tight-overestimation policies and their final
objective intervals.

Connector numerical preprocessing is bypassed during certified inexact
pricing because fixing variables or temporarily replacing the objective would
make the solver's bound invalid for the original pricing problem. Pareto and
BlC coefficient-refinement stages are skipped whenever the returned pricing
bounds show a nonzero connector error. If an infeasible first-level point is
encountered, JuBiC calls the same subsolver in exact mode to generate a proven
feasibility cut. Fixed-``x`` evaluation always uses the ordinary exact
`solve_sub_for_x` method and accepts any follower-optimal response.

For the final GBC incumbent, `SubSolverJuMP` also implements
`solve_sub_for_x_optimistic`. It first solves the fixed-``x`` model to proven
optimality for `c_objterm`, adds an equality fixing that optimal follower value,
and then minimizes `r_objterm`. Thus the returned response minimizes the
first-level contribution over the complete follower-optimal set. This
lexicographic two-solve formulation does not use an epsilon or weighted penalty
term, avoiding problem-dependent scaling and numerical ambiguity. Both stages
force zero relative and absolute MIP gaps and restore the original settings.

The connector initialization value returned by
`compute_lower_bound_master_contribution` must be a valid global lower bound.
JuBiC therefore forces `MIPGap = MIPGapAbs = 0` for this solve and requires an optimal
termination status, even when the `SubSolverJuMP` instance enables heuristic
pricing or its underlying model was configured with a nonzero MIP gap. The
previous solver setting is restored afterwards. The `heuristic` and `mip_gap`
options affect only ordinary connector-pricing calls.

## Numerical Preprocessing in `GBC`

`SubSolverJuMP` supports the `GBC` numerical preprocessing option for extreme
connector coefficients. The behavior is documented in
[Numerics and Status Codes](../numerics_and_status.md).

## Related Pages

- [`GBC`](../solvers/gbc.md)
- [`BlC`](../solvers/blc.md)
- [SubSolver Interface](../sub_solvers.md)
