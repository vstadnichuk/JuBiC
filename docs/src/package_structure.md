# Package Structure

The root module is defined in `src/JuBiC.jl`. It includes and exports the core object types, solver drivers, wrappers, and utility functions.

## Source Layout

The internal package is organized into the following groups:

## Solver Parameter and Wrapper Layer

- `src/solvers/solver_parameters.jl`
- `src/solvers/solver_wrapper.jl`
- `src/solvers/solver_statistics.jl`

This layer provides:

- solver parameter structs,
- wrapped external optimizer environments,
- and result/statistics storage.

## Model Objects

- `src/model_structs/master.jl`
- `src/model_structs/blc_master.jl`
- `src/model_structs/blc_lag_master.jl`
- `src/model_structs/mip_master.jl`
- `src/model_structs/mibs_master.jl`
- `src/model_structs/instance.jl`
- `src/model_structs/sub_solver.jl`

This layer defines the objects that are composed into a JuBiC `Instance`.

## Connector LPs

- `src/model_structs/connector_lp.jl`
- `src/model_structs/connector_lp_blc.jl`

These files implement the iterative separation models used by `GBC` and `BlCLag`.

## Solver Drivers

- `src/solvers/gbc_solver.jl`
- `src/solvers/blc_solver.jl`
- `src/solvers/blclag_solver.jl`
- `src/solvers/mip_solver.jl`
- `src/solvers/mibs_solver.jl`
- `src/solvers/solvers.jl`

This layer contains the algorithm entry points and the dispatch from `solve_instance!`.

## Built-in Subsolver Implementations

- `src/solvers/sub_solvers/sub_solver_mip.jl`
- `src/solvers/sub_solvers/sub_solver_blc_mip.jl`
- `src/solvers/sub_solvers/sub_solver_mibs.jl`
- `src/solvers/sub_solvers/a_star_search.jl`
- `src/solvers/sub_solvers/labeling.jl`

## MiBS Parsing and File Generation

- files under `src/mibs_parser`

These files implement the transformation to and from the `MiBS` file interface.

## Utilities

- `src/Auxilliaries.jl`
- `src/utils/logging_utils.jl`
- `src/utils/gbc_instance_io.jl`
- `src/solvers/batch_solve.jl`

## Dependency Roles

The current implementation uses:

- `JuMP` for model construction,
- `MathOptInterface` for status and callback integration,
- `Gurobi` through `SolverWrapper`,
- and `BilevelJuMP` for the direct `MiBS` master wrapper.

The package documentation below follows the same internal structure: object layer, solver API, solver methods, subsolvers, then numerical behavior.
