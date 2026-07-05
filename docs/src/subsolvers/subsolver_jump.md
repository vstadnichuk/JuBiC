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
    extra_cuts,
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

Deprecated constructors with explicit linking-copy variables still exist for
backward compatibility, but JuBiC now ignores the passed copies and creates its
own internal copies with unit linking capacities.

## Numerical Preprocessing in `GBC`

`SubSolverJuMP` supports the `GBC` numerical preprocessing option for extreme
connector coefficients. The behavior is documented in
[Numerics and Status Codes](../numerics_and_status.md).

## Related Pages

- [`GBC`](../solvers/gbc.md)
- [`BlC`](../solvers/blc.md)
- [SubSolver Interface](../sub_solvers.md)
