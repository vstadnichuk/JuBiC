# HNDP Solver Overview

This note summarizes the HNDP model-generation functions in
[hndp_model_generation.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/hndp_model_generation.jl)
and the JuBiC solver parameter objects in
[solver_parameters.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/src/solvers/solver_parameters.jl).

The structure is:
- common HNDP concepts
- common JuBiC solve parameters
- solver/model-specific generation options

## Common Concepts

The HNDP pipeline separates two layers of configuration:
- model generation:
  which mathematical reformulation is built from an `HNDPwC` instance
- solver execution:
  which JuBiC solver parameter object is used when calling `solve_instance!`

Typical workflow:
1. Build or load an `HNDPwC` instance.
2. Build a JuBiC instance with one of the `build_hndp_*_instance(...)` functions.
3. Solve it with the matching JuBiC parameter object such as `BLCparam`, `GBCparam`, `MIPparam`, or `MibSparam`.

## Shared HNDP Option Sets

### Big-M Modes

These constants control how user-specific big-M values are derived:

- `HNDP_BIGM_FIXED_NETWORK_PATH`
  Uses the follower cost of the shortest or constrained-shortest path in the fixed network.
  This is usually the tighter and more problem-aware option.

- `HNDP_BIGM_N_MINUS_ONE`
  Uses the `n - 1` most expensive arcs heuristic.
  This is the fallback-style option and does not require a feasible fixed-network path.

These modes are used by:
- `build_hndp_blc_instance`
- `build_hndp_sd_instance`
- `build_hndp_sd_auto_instance`
- `build_hndp_hybrid_instance` for fallback users
- `build_hndp_hybrid_blc_instance` for fallback users
- `build_hndp_gbc_instance` only when `subproblem_method == HNDP_SUBPROBLEM_BLC_JUMP`

### Subproblem Methods

These constants control how follower problems are solved inside decomposition-type models:

- `HNDP_SUBPROBLEM_MIP`
  Classical JuMP-based MIP follower model.

- `HNDP_SUBPROBLEM_ASTAR`
  AStar / labeling based shortest-path or constrained-shortest-path follower solver.

- `HNDP_SUBPROBLEM_BLC_JUMP`
  JuMP-based bilevel-aware follower subsolver that keeps persistent BlC-style cuts.
  This is only supported for GBC generation.

### Hybrid Fallback Modes

Currently implemented:

- `HNDP_HYBRID_FALLBACK_SD`
  Uses strong duality as fallback for users whose path enumeration times out.

This is currently the only supported fallback mode in the generic hybrid path builder.

## Common JuBiC Solve Parameters

The JuBiC parameter structs are independent of HNDP and live in
[src/solvers/solver_parameters.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/src/solvers/solver_parameters.jl).

### Parameters Shared By Most JuBiC Solvers

These fields appear in most JuBiC solver parameter structs:

- `solver`
  The `SolverWrapper`, typically `GurobiSolver()`.

- `debbug_out`
  If `true`, extra logs and model files are written.

- `output_folder_path`
  Folder for logs, exported models, and solution files.

- `file_format_output`
  Usually `"lp"` or `"mps"` for JuMP model exports.

- `stats`
  A `RunStats` object where run information is collected.

- `runtime`
  Requested wall-clock limit in seconds.

- `seed`
  Random seed forwarded to the underlying Gurobi-based MIP solves.
  The convenience constructors currently default this to `42`.

### Thread Settings

For GBC, BlC, and BlCLag, JuBiC distinguishes:

- `threads_master`
  Threads for the master MIP

- `threads_sub_con`
  Threads for JuMP-based follower subproblems and connector LP routines

For `MIPparam`, only `threads_master` is relevant because the model is solved as a single MIP.

### Practical Reproducibility Note

If reproducibility matters, use:
- a fixed `seed`
- fixed thread counts
- ideally small thread counts, often `1`, if you want the most stable run-to-run behavior

The seed is now propagated to:
- master Gurobi MIPs
- JuMP-based follower MIPs
- JuMP-based bilevel follower oracles

## Solver Parameter Objects

### `GBCparam`

Use with instances returned by `build_hndp_gbc_instance(...)`.

Important fields:
- `runtime`
- `seed`
- `threads_master`
- `threads_sub_con`
- `pareto`
  One of:
  - `PARETO_NONE`
  - `PARETO_OPTIMALITY_ONLY`
  - `PARETO_OPTIMALITY_AND_FEASIBILITY`
- `warmstart`
  If `false`, connector LP objects are reset more aggressively.
- `bigMwithLC`
  Enables stronger big-M handling through additional connector-LP work.
- `trim_coeff`
  Allows coefficient trimming in generated cuts.
- `infinity_num`
  Numeric surrogate for infinity in some bounded reformulations.
- `g_round_digit`
  Numeric rounding parameter used in connector-related logic.

Typical convenience usage:

```julia
param = GBCparam(solver, false, outdir, "lp", PARETO_OPTIMALITY_ONLY, 3600)
```

This uses:
- `seed = 42`
- `threads_master = 8`
- `threads_sub_con = 8`

### `BLCparam`

Use with instances returned by `build_hndp_blc_instance(...)` or
`build_hndp_hybrid_blc_instance(...)`.

Important fields:
- `runtime`
- `seed`
- `threads_master`
- `threads_sub_con`

Typical convenience usage:

```julia
param = BLCparam(solver, false, outdir, "lp", 3600)
```

### `BlCLagparam`

Use with BlCLag-compatible bilevel instances. This is not currently the main
HNDP route, but it remains part of JuBiC core.

Important fields:
- `runtime`
- `seed`
- `threads_master`
- `threads_sub_con`
- `pareto`
- `warmstart`
- `infinity_num`

### `MIPparam`

Use with compact reformulations such as:
- strong duality
- path model
- hybrid path model

Important fields:
- `runtime`
- `seed`
- `threads_master`

Typical convenience usage:

```julia
param = MIPparam(solver, false, outdir, "lp", 3600)
```

### `MibSparam`

Use with instances returned by `build_hndp_mibs_instance(...)`.

Important fields:
- `debbug_out`
- `output_folder_path`
- `stats`
- `runtime`

Current behavior:
- JuBiC uses its own MiBS execution pipeline, not the old plain BilevelJuMP call.
- JuBiC writes a MiBS parameter file and forwards `runtime` via `Alps_timeLimit`.
- Before exporting, JuBiC validates that the transformed bilevel model is MIP-MIP and raises a more informative error if non-integer variables are found.

Typical convenience usage:

```julia
param = MibSparam(false, outdir, 3600.0)
```

## HNDP Model Builders

## `build_hndp_blc_instance`

Signature:

```julia
build_hndp_blc_instance(
    hndp,
    solver;
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    subproblem_method=HNDP_SUBPROBLEM_MIP,
)
```

Purpose:
- Builds the classical Benders-like decomposition reformulation.

Supported options:
- `big_m_mode`
  - `HNDP_BIGM_FIXED_NETWORK_PATH`
  - `HNDP_BIGM_N_MINUS_ONE`
- `subproblem_method`
  - `HNDP_SUBPROBLEM_MIP`
  - `HNDP_SUBPROBLEM_ASTAR`

Use with:
- `BLCparam`

## `build_hndp_gbc_instance`

Signature:

```julia
build_hndp_gbc_instance(
    hndp,
    solver;
    partial_decomposition=true,
    include_objL2=false,
    subproblem_method=HNDP_SUBPROBLEM_MIP,
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    enforce_integer_construction_cost=false,
)
```

Purpose:
- Builds the GBC reformulation.

Shared options:
- `partial_decomposition`
  If `true`, follower flow constraints are injected into the master through the master callback setup.
- `include_objL2`
  Exposes second-level objective variables in the master.
  Only valid when `partial_decomposition=true`.

Follower options:
- `subproblem_method`
  - `HNDP_SUBPROBLEM_MIP`
  - `HNDP_SUBPROBLEM_ASTAR`
  - `HNDP_SUBPROBLEM_BLC_JUMP`

Big-M relevance:
- `big_m_mode` matters only for `HNDP_SUBPROBLEM_BLC_JUMP`

Special MiBS-related option:
- `enforce_integer_construction_cost`
  Intended for the MiBS path.
  Forces `construction_cost_var` to integer so the transformed bilevel model remains MIP-MIP.

Use with:
- `GBCparam`

## `build_hndp_mibs_instance`

Signature:

```julia
build_hndp_mibs_instance(hndp, solver; partial_decomposition=false)
```

Purpose:
- Builds an HNDP instance for MiBS by:
  - first constructing a GBC-style JuMP decomposition model
  - then merging followers if needed
  - then transforming the result into a single-follower MiBS bilevel model

Important behavior:
- `partial_decomposition` is forwarded to the underlying GBC construction
- integer construction cost is enforced automatically for MiBS compatibility
- multi-follower instances are merged with `merge_subproblems(...)`

Use with:
- `MibSparam`

## `build_hndp_sd_instance`

Signature:

```julia
build_hndp_sd_instance(
    hndp,
    solver;
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    indicator_constraints=false,
    bound_duals=true,
)
```

Purpose:
- Builds the explicit strong-duality reformulation.

Restrictions:
- only supported for instances without weight bounds

Options:
- `big_m_mode`
  - `HNDP_BIGM_FIXED_NETWORK_PATH`
  - `HNDP_BIGM_N_MINUS_ONE`
- `indicator_constraints`
  If `true`, uses indicator constraints for some dual feasibility logic.
- `bound_duals`
  If `true`, bounds the dual potentials using the derived big-M values.

Use with:
- `MIPparam`

## `build_hndp_sd_auto_instance`

Signature:

```julia
build_hndp_sd_auto_instance(
    hndp,
    solver;
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
)
```

Purpose:
- Experimental automated strong-duality route based on generic JuBiC wrappers.

Restrictions:
- only supported for instances without weight bounds
- currently experimental rather than part of the stable HNDP pipeline

Use with:
- `MIPparam`

## `build_hndp_path_instance`

Signature:

```julia
build_hndp_path_instance(
    hndp,
    solver;
    enumeration_time_limit,
    parallelize=false,
)
```

Purpose:
- Builds the path-based compact reformulation.

Options:
- `enumeration_time_limit`
  Global time budget for path-bound computation and path enumeration.
- `parallelize`
  If `true`, uses the parallel precomputation variant.

Important behavior:
- if no feasible fixed-network path exists, JuBiC falls back to the `n - 1` bound and logs a warning
- the parallel variant interprets `enumeration_time_limit` as a wall-clock budget for precomputation only

Use with:
- `MIPparam`

## `build_hndp_hybrid_instance`

Signature:

```julia
build_hndp_hybrid_instance(
    hndp,
    solver;
    enumeration_time_limit,
    fallback_mode=HNDP_HYBRID_FALLBACK_SD,
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    indicator_constraints=false,
    bound_duals=true,
)
```

Purpose:
- Builds the hybrid path model:
  users solved by path enumeration stay in the compact path formulation;
  users that time out fall back to another formulation.

Options:
- `enumeration_time_limit`
  Global wall-clock deadline for parallel precomputation.
- `fallback_mode`
  Currently only `HNDP_HYBRID_FALLBACK_SD`.
- `big_m_mode`
  Used for fallback users.
- `indicator_constraints`
  Passed to the SD fallback model.
- `bound_duals`
  Passed to the SD fallback model.

Use with:
- `MIPparam`

## `build_hndp_hybrid_blc_instance`

Signature:

```julia
build_hndp_hybrid_blc_instance(
    hndp,
    solver;
    enumeration_time_limit,
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    subproblem_method=HNDP_SUBPROBLEM_MIP,
)
```

Purpose:
- Builds a mixed formulation:
  users with fully enumerated path sets are represented compactly in the HPR,
  and only timed-out users stay as true BlC subproblems.

Options:
- `enumeration_time_limit`
  Global wall-clock deadline for path precomputation.
- `big_m_mode`
  Used only for fallback users that remain subproblems.
- `subproblem_method`
  - `HNDP_SUBPROBLEM_MIP`
  - `HNDP_SUBPROBLEM_ASTAR`

Use with:
- `BLCparam`

## Recommended Pairings

Common model / parameter combinations:

- BlC model:
  `build_hndp_blc_instance(...)` + `BLCparam(...)`

- GBC model:
  `build_hndp_gbc_instance(...)` + `GBCparam(...)`

- MiBS model:
  `build_hndp_mibs_instance(...)` + `MibSparam(...)`

- Strong duality:
  `build_hndp_sd_instance(...)` + `MIPparam(...)`

- Path model:
  `build_hndp_path_instance(...)` + `MIPparam(...)`

- Hybrid SD fallback:
  `build_hndp_hybrid_instance(...)` + `MIPparam(...)`

- Hybrid BlC fallback:
  `build_hndp_hybrid_blc_instance(...)` + `BLCparam(...)`

## Notes And Caveats

- MiBS currently requires a MIP-MIP bilevel model.
  JuBiC now checks this explicitly before export and gives a variable-level error message when possible.

- For MiBS generation in HNDP, JuBiC currently assumes the construction-cost expression is integral and therefore forces the helper variable `construction_cost_var` to be integer.

- AStar-based subproblem methods are available only where explicitly listed above.

- Strong-duality-based formulations are only valid for users without weight bounds unless explicitly modeled through another route.
