# `SubSolverMiBS`

`SubSolverMiBS` uses the same bilevel-separation idea as
[`SubSolverBlCJuMP`](subsolver_blc_jump.md), but delegates the resulting
subproblem solves to [`MiBS`](https://github.com/coin-or/MibS).

## Constructor Parameters

The main constructor is:

```julia
SubSolverMiBS(
    name,
    bi_model,
    A,
    y_vars,
    r_objterm,
    c_objterm,
)
```

The parameters are:

- `name`: unique name of the follower problem. It must match the corresponding entry in the master's `sub_names`.
- `bi_model`: `BilevelJuMP.BilevelModel` representing the follower-side bilevel subproblem solved by MiBS.
- `A`: shared resource set used for the linking variables. The keys of `y_vars` must match this set.
- `y_vars`: lower-level variables appearing in the linking constraints. For every `a in A`, JuBiC creates an upper-level binary copy `x_copy[a]` and adds the lower-level linking constraint `y_vars[a] <= x_copy[a]`.
- `r_objterm`: expression for the follower solution's contribution to the first-level objective.
- `c_objterm`: expression for the original follower objective. This objective remains in the lower level when MiBS solves the subproblem.

`SubSolverMiBS` declares support for bilevel subproblem separation:

```julia
supports_bilevel_subproblem_solver(::SubSolverMiBS) = true
```

This is why it can be used with `BlCLag`.

## Related Pages

- [`BlCLag`](../solvers/blc_lag.md)
- [`Direct MiBS Wrapper`](../solvers/mibs.md)
- [SubSolver Interface](../sub_solvers.md)
