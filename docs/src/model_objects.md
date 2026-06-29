# Using JuBiC

The central idea is that JuBiC does not ask for one universal bilevel model representation for every algorithm. Instead, it provides a small set of wrapper objects that match different solver ideas.

## Why JuBiC Uses Explicit Wrapper Objects

Several JuBiC solvers are decomposition-based. In those methods, the first level and the follower problems play different algorithmic roles:

- the first level is handled by a master problem,
- the followers are handled through one or more subproblems,
- and the solver exchanges information between them through cuts.

That is why JuBiC introduces explicit wrapper objects.

This split is especially useful when:

- the same first-level decisions affect many follower problems,
- and we want to exploit algorithmically the structure of follower problems.

Not every solver in JuBiC uses decomposition, but this master/subproblem view is the main organizational idea of the package.

## Technical Structure: `Instance`, `Master`, and `SubSolver`

At the technical level, JuBiC organizes a solve around three object layers:

- `Instance`
- one master wrapper
- zero or more subsolver wrappers

The basic pattern is:

- `Instance(master, [sub1, sub2, ...])` for decomposition-based solvers
- `Instance(master, nothing)` for direct solver routes without JuBiC subproblems

The master wrapper stores the first-level representation used by the selected solver. Subsolver wrappers store the follower-side models or follower-side algorithms.

The role of `Instance` is to collect these pieces into one object that can be passed to JuBiC's common solve interface, `solve_instance!`.

## Solver Families and Their Modeling View

At a high level, the included solvers can be grouped as follows.

### Decomposition Solvers

- [`GBC`](solvers/gbc.md)
  - Uses logic-based Benders cuts in the natural space of the first-level variables.
  - Modeled with `Master` plus follower subsolvers.

- [`BlC`](solvers/blc.md)
  - Uses a high-point relaxation and follower subproblems.
  - The user provides the big-M values used in the Benders-like cuts.
  - Modeled with `BlCMaster` plus follower subsolvers.

- [`BlCLag`](solvers/blc_lag.md)
  - Uses the `BlC` structure, but tries to generate big-M information from the given problem structure.
  - Modeled with `BlCLagMaster` plus bilevel-capable follower subsolvers.

### Single-Model Solver Routes

- [Compact MIP wrapper](solvers/mip.md)
  - This is JuBiC's generic wrapper for solving compact MIP models directly.
  - Modeled with `MIPMaster`.

- [Direct `MiBS`](solvers/mibs.md)
  - JuBiC provides a direct wrapper for [`MiBS`](https://github.com/coin-or/MibS).
  - Built around [`BilevelJuMP`](https://joaquimg.github.io/BilevelJuMP.jl/stable/) and modeled with `MibSMaster`.

So the wrapper you choose is tightly connected to the algorithm you intend to run.

## Master Wrappers

JuBiC provides several master-side wrapper types. Their fields and intended use are documented on the solver pages.

Decomposition solver wrappers:

- [`Master` for `GBC`](solvers/gbc.md)
- [`BlCMaster` for `BlC`](solvers/blc.md)
- [`BlCLagMaster` for `BlCLag`](solvers/blc_lag.md)

Single-model solver wrappers:

- [`MIPMaster` for the compact MIP wrapper](solvers/mip.md)
- [`MibSMaster` for the direct `MiBS` wrapper](solvers/mibs.md)

## Subsolver Wrappers

All JuBiC follower-side wrappers inherit from the abstract type `SubSolver`.

The required interface is documented in [SubSolver Interface](sub_solvers.md). Built-in subsolvers are grouped by whether the follower oracle itself needs to solve a bilevel problem.

Non-bilevel subsolvers:

- [`SubSolverJuMP`](subsolvers/subsolver_jump.md)
  - Generic JuMP-based follower wrapper for single-level follower solves.
- [`AStarSolver`](subsolvers/astar.md)
  - Labeling / A*-style follower oracle for shortest-path type problems.

Bilevel-capable subsolvers:

- [`SubSolverBlCJuMP`](subsolvers/subsolver_blc_jump.md)
  - Generic JuMP-based follower wrapper where the follower problem is modeled as JuMP model.
- [`SubSolverMiBS`](subsolvers/subsolver_mibs.md)
  - File-based follower wrapper that delegates bilevel subproblem solves to `MiBS`.

## MIP Solver

Most JuBiC routines rely on an available MIP solver. The direct `MiBS` wrapper is currently the main exception, because the full bilevel solve is delegated to `MiBS`.

MIP solvers are represented through [`SolverWrapper`](solver_api.md) subtypes. Currently, `GurobiSolver` is the only MIP solver wrapper that has been extensively tested. The implementation is generic, so other MIP solvers can be used by implementing a custom `SolverWrapper` object.

For callback-based routines, JuBiC expects the selected MIP solver to support [JuMP's solver-independent callback interface](https://jump.dev/JuMP.jl/stable/manual/callbacks/).

## Statistics and Experiment Output

`solve_instance!` returns a `RunStats` object. It stores the main status, objective value, runtime, and solver-specific statistics collected during the solve.

JuBiC also provides benchmark pipeline utilities for automated experiment runs. These pipelines execute generated instances through the same solver interface and write collected statistics to CSV files. The HNDP documentation shows this workflow in more detail in [HNDP Benchmark Pipeline](examples/hndp/benchmarks.md).

## Typical Decomposition Build Pattern

A generic decomposition-based workflow looks like this:

```julia
master = Master(master_model, A, xdict, ["Sub0"])
sub = SubSolverJuMP("Sub0", follower_model, A, y, first_level_term, follower_obj, timelimit -> (false, 0))
instance = Instance(master, [sub])
params = GBCparam(GurobiSolver(), false, outdir, "lp", PARETO_OPTIMALITY_ONLY, 60.0)
stats = solve_instance!(instance, params)
```
