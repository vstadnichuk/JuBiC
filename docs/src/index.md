# JuBiC Documentation

## Problem Class

JuBiC addresses mixed-integer bilevel optimization problems with binary linking variables. A generic problem in the supported class can be written as

```math
\begin{aligned}
\min_{x,\, y^1,\dots,y^K}\quad & F(x) + \sum_{k \in \mathcal{K}} G_k(x, y^k) \\
\text{s.t.}\quad & x \in X, \\
& y^k \in \arg\min_{\bar y \in Y_k(x)} f_k(x, \bar y) \qquad \forall k \in \mathcal{K},
\end{aligned}
```

where:

- `x` denotes the first-level decision variables,
- `\mathcal{K}` is the set of second-level problems,
- `Y_k(x)` is the feasible region of follower `k`, parameterized by the first-level linking decisions,
- `F` and `G_k` define the first-level objective,
- and `f_k` is the follower objective.

JuBiC builds on top of [JuMP](https://jump.dev/JuMP.jl/stable/), so the first-level model `X` and the follower models `Y_k(x)` are ordinary JuMP models. The important structural restriction for the algorithms currently implemented in JuBiC is that the linking variables connecting both levels are binary. In a standard form, this is represented by

```math
x_a \in \{0,1\} \qquad \forall a \in A,
```

and, for each follower `k`,

```math
y^k_a \in \{0,1\}, \qquad y^k_a \le x_a \qquad \forall a \in A.
```

This coupling means that the first level activates or blocks second-level binary decisions through the shared resource set `A`.

The following are the settings that JuBiC has been tested on so far:

- a minimization problem at both levels,
- linear first-level and second-level constraints,
- linear objective terms,
- and binary linking variables that activate or restrict second-level decisions.

JuBiC provides several solver interfaces for this problem class. The package is
organized around a small number of common concepts:

- **model wrappers**, which describe how the first-level model and follower problems are passed to JuBiC,
- **solver methods**, which define the algorithmic route used to solve an `Instance`,
- **subsolvers**, which implement follower-side oracle calls for decomposition methods,
- and **run statistics**, which report solver status, objective values, runtimes, and numerical fallback information.

The documentation follows this structure. Start with a small complete example,
then read the wrapper overview, and then move to the solver or subsolver page
that matches the algorithm you want to use.

## Installation

The current development version can be installed directly from GitHub:

```julia-repl
julia> import Pkg

julia> Pkg.add(url="https://github.com/vstadnichuk/JuBiC")
```

## Verify Installation

To run the package tests:

```julia-repl
julia> Pkg.test("JuBiC")
```

The current test suite uses `Gurobi`, so running the tests requires a working `Gurobi` installation and a valid license.

## Documentation Map

- [Getting Started](getting_started.md)
  A small bilevel model written in JuBiC syntax and solved with representative solver routes.
- [Using JuBiC](model_objects.md)
  The main conceptual overview: `Instance`, master wrappers, subsolver wrappers, solver families, and experiment output.
- [Solver Methods](solvers/gbc.md)
  Algorithm-specific pages for `GBC`, `BlC`, `BlCLag`, the compact MIP wrapper, and the direct `MiBS` wrapper.
- [SubSolvers](sub_solvers.md)
  The follower-oracle interface and the currently implemented subsolver wrappers.
- [Numerics and Status Codes](numerics_and_status.md)
  Numerical safeguards, fallback logic, and the status values written to `RunStats`.
- [Examples](examples/hndp/motivation.md)
  Application examples built on top of the JuBiC package, starting with HNDP.
- [Core Solver API](solver_api.md)
  A generated technical reference for exported JuBiC objects.

## Citation

If you use JuBiC in research, cite:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online. [Link](https://optimization-online.org/?p=28877)
