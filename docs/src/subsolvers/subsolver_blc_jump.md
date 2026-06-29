# `SubSolverBlCJuMP`

`SubSolverBlCJuMP` lets a follower problem be passed to JuBiC as a JuMP model
while still supporting bilevel-aware separation routines.

The motivation is that routines such as `separation_BlC!` change the objective
inside the subsolver. Even under this changed objective, JuBiC must return a
solution that is optimal for the original follower objective. `SubSolverBlCJuMP`
handles this by using a multi-tree implementation of BlC inside the subsolver:
the working MIP is solved with the current separation objective, then an oracle
solves the original follower problem for the same linking-variable pattern. If
the candidate is not optimal for the original follower objective, a persistent
BlC cut is added and the working MIP is resolved.

The theoretical idea is the same BlC concept described on the
[`BlC`](../solvers/blc.md) solver page, but applied inside a subsolver.

## Constructor Parameters

The main constructor is:

```julia
SubSolverBlCJuMP(
    name,
    mip_model,
    A,
    y_vars,
    r_objterm,
    c_objterm,
    big_m,
)
```

There is also a constructor that accepts a custom `oracle_solve` callback:

```julia
SubSolverBlCJuMP(
    name,
    mip_model,
    A,
    y_vars,
    r_objterm,
    c_objterm,
    oracle_solve,
    big_m,
)
```

The parameters are:

- `name`: unique name of the follower problem. It must match the corresponding entry in the master's `sub_names`.
- `mip_model`: JuMP model used as the working MIP. JuBiC modifies this model during the solve by changing objectives, adding temporary fixing constraints, adding internal linking-copy variables, and adding persistent BlC cuts.
- `A`: shared resource set used for the linking variables. The keys of `y_vars` must match this set.
- `y_vars`: follower-side variables appearing in the linking constraints. For every `a in A`, JuBiC creates an internal binary copy `x_copy[a]` and adds `y_vars[a] <= x_copy[a]`.
- `r_objterm`: expression for the follower solution's contribution to the first-level objective.
- `c_objterm`: expression for the original follower objective. This is the objective used by the oracle to verify lower-level optimality.
- `oracle_solve`: optional callback that solves the original follower problem for fixed linking-variable values. If omitted, JuBiC builds a default oracle by copying `mip_model` and minimizing `c_objterm`.
- `big_m`: big-M callback used to generate persistent BlC cuts inside the subsolver. This is the same type of coefficient used in the `BlC` solver: `big_m(a)` must return a valid coefficient for resource `a`.

`SubSolverBlCJuMP` declares support for bilevel subproblem separation:

```julia
supports_bilevel_subproblem_solver(::SubSolverBlCJuMP) = true
```

This is why it can be used with `BlCLag`.

## Related Pages

- [`BlC`](../solvers/blc.md)
- [`BlCLag`](../solvers/blc_lag.md)
- [SubSolver Interface](../sub_solvers.md)
