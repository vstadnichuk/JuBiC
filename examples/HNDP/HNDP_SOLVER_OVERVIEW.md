# HNDP Solver Overview

This note summarizes the HNDP model-generation functions in
[hndp_model_generation.jl](hndp_model_generation.jl)
and the JuBiC solver parameter objects in
[solver_parameters.jl](../../src/solvers/solver_parameters.jl).

The structure is:
- common HNDP concepts
- common JuBiC solve parameters
- solver/model-specific generation options

## Formulations at a glance

Run examples from the repository root with `julia --project=.`. The HNDP
builders are example code loaded with `include`, rather than exported JuBiC APIs.
They accept an `HNDPwC` instance and a solver wrapper. The runner selects them
using the `model_type` field below.

| Formulation | `model_type` | Builder | Solve parameters | Main restriction or behavior |
| --- | --- | --- | --- | --- |
| Benders-like cuts (BlC) | `blc` | `build_hndp_blc_instance` | `BLCparam` | MIP or A* follower oracle |
| Lagrangian Benders-like cuts (BlCLag) | `blclag` | `build_hndp_blclag_instance` | `BlCLagparam` | Persistent BlC-JuMP or MiBS follower oracle |
| Generalized Benders cuts (GBC) | `gbc` | `build_hndp_gbc_instance` | `GBCparam` | MIP, A*, or BlC-JuMP oracle; optional partial decomposition |
| MiBS | `mibs` | `build_hndp_mibs_instance` | `MibSparam` | Requires an installed MiBS executable and a MIP-MIP model |
| Explicit strong duality | `sd` | `build_hndp_sd_instance` | `MIPparam` | No follower weight bounds |
| Enumerated paths | `path` | `build_hndp_path_instance` | `MIPparam` | Enumerates admissible paths; generation can be expensive |
| Paths with strong-duality fallback | `hybrid_sd` | `build_hndp_hybrid_instance` | `MIPparam` | Non-enumerated users need an unweighted shortest-path fallback |
| Paths with BlC fallback | `hybrid_blc` | `build_hndp_hybrid_blc_instance` | `BLCparam` | Keeps fallback users as BlC subproblems |
| Automated strong duality (experimental) | Direct Julia API only | `build_hndp_sd_auto_instance` | `MIPparam` | No weight bounds; experimental generic reformulation |

The path builder returns `(instance, enumeration_seconds, path_counts)`.
Both hybrid builders return five values (also fallback user lists). The other builders return the JuBiC instance directly. The experiment
runner handles these differences and deducts path-generation time from the
remaining solve budget. Full path enumeration is not a substitute for a
completed generation phase; use a hybrid formulation when fallback is needed.

See [NETWORK_GENERATION.md](NETWORK_GENERATION.md) for instance families,
topologies, and construction examples. A weight constraint is present when a
user's `weighlimit` is non-`nothing`; `length_constrained=false` omits it.

## Small runnable solve example

Start Julia from the repository root with `julia --project=. --threads=2`.
This constructs a small unconstrained competition instance and solves the
explicit strong-duality model. A working Gurobi installation/license is required.

```julia
using JuBiC
include("examples/HNDP/hndp_model_generation.jl")

config = Dict{String,Any}("parameter_seeds" => [17], "instances" => [
    Dict{String,Any}("name" => "small_competition", "instance_type" => "competition",
        "topologies" => ["layered_sioux_falls"], "nusers" => [2],
        "length_constrained" => [false])
])
generated = only(generate_hndp_networks(config))
solver = GurobiSolver()
instance = build_hndp_sd_instance(generated.instance, solver)
mkpath("tmp/hndp_sd_example")
params = MIPparam(solver, false, "tmp/hndp_sd_example", "lp", 30.0)
stats = solve_instance!(instance, params)
```

For BlC on the same instance, use
`build_hndp_blc_instance(generated.instance, solver)` and
`BLCparam(solver, false, "tmp/hndp_blc_example", "lp", 30.0)`.
For GBC, use `build_hndp_gbc_instance` and `GBCparam` with the same arguments.
Create the chosen output directory with `mkpath` before direct solver calls.
The remaining builder signatures and options are documented below.

To run the supplied bike and transit-expansion examples through the generic
runner (or append `--dry-run` to construct them without optimization):

```sh
julia --project=. --threads=2 examples/HNDP/run_hndp_manifest.jl examples/HNDP/multimodal_example.json tmp_compare/hndp_multimodal_example
```

The JSON contains `instances`, `models`, and `params`. Replace the model entry,
for example, with `{"name":"gbc_astar","model_type":"gbc","subproblem_method":"astar"}`
to select GBC with A*. Solver-specific options belong in the `params` entry.

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

- `HNDP_BIGM_FIXED_NETWORK_PATH_CURRENT_COST`
  Uses a fixed-network bound with the current follower cost and a numerical
  margin during cut generation. Select `"fixed_network_path_current_cost"`
  in JSON. This is a dynamic cut callback mode, not a new solver.

These modes are used by:
- `build_hndp_blc_instance`
- `build_hndp_sd_instance`
- `build_hndp_sd_auto_instance`
- `build_hndp_hybrid_instance` for fallback users
- `build_hndp_hybrid_blc_instance` for fallback users
- `build_hndp_gbc_instance` only when `subproblem_method == HNDP_SUBPROBLEM_BLC_JUMP`

For layered competition instances, `HNDP_BIGM_FIXED_NETWORK_PATH` now uses the
two-layer competitor graph. For single-layer `k`-arc instances it uses the
fixed-arc subgraph and falls back to the `n - 1` heuristic when no feasible
fixed-network path exists.

### Subproblem Methods

These constants control how follower problems are solved inside decomposition-type models:

- `HNDP_SUBPROBLEM_MIP`
  Classical JuMP-based MIP follower model.

- `HNDP_SUBPROBLEM_ASTAR`
  AStar / labeling based shortest-path or constrained-shortest-path follower solver.

- `HNDP_SUBPROBLEM_BLC_JUMP`
  JuMP-based bilevel-aware follower subsolver that keeps persistent BlC-style cuts.
  Supported by GBC and BlCLag generation.

- `HNDP_SUBPROBLEM_MIBS`
  MiBS follower oracle for BlCLag; select `"mibs"` as its `subproblem_method`.

### Hybrid Fallback Modes

Currently implemented:

- `HNDP_HYBRID_FALLBACK_SD`
  Uses strong duality as fallback for users whose path enumeration times out.

This is currently the only supported fallback mode in the generic hybrid path builder.

The hybrid builders also support a second fallback trigger:

- users may be pushed to fallback even after successful enumeration when the
  number of enumerated paths exceeds the number of arc-flow variables in the
  fallback formulation

These path-count-triggered fallbacks are tracked separately from timeout
fallbacks in the benchmark metadata.

## Common JuBiC Solve Parameters

The JuBiC parameter structs are independent of HNDP and live in
[src/solvers/solver_parameters.jl](../../src/solvers/solver_parameters.jl).

### Parameters Shared By Most JuBiC Solvers

These fields appear in most JuBiC solver parameter structs:

- `solver`
  The `SolverWrapper`, typically `GurobiSolver()`.

- `debbug_out`
  If `true`, extra logs and model files are written.

- `output_folder_path`
  Folder for logs, exported models, and solution files.

- `enable_output_logs`
  Internal flag propagated through the HNDP runner. When disabled, streamed CSV
  output still remains active, but solver-side per-instance debug logs and
  model exports are skipped.

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
  Gurobi threads assigned to each follower subproblem. For GBC, ConnectorLP
  models always use one Gurobi thread; this setting controls the follower
  subsolver. BlC and BlCLag use this setting for their follower subsolver.

- `parallel_subsolvers`
  Number of follower separation workers allowed to run concurrently. Each
  worker owns one follower subsolver and its ConnectorLP. The worker uses the
  configured `threads_sub_con` threads for the follower subsolver. ConnectorLP
  uses one thread in parallel GBC separation and shares the worker's thread
  budget because it and the follower subsolver are not solved simultaneously.

- `parallel_separation`
  If `true`, follower workers are scheduled concurrently. If `false`, one
  follower is separated at a time. The configured follower thread count is
  still used; GBC ConnectorLP remains fixed at one thread.

The requested budget for parallel separation is calculated as
`threads_master + parallel_subsolvers * max(threads_sub_con, connector_threads)`.
If this exceeds the Julia thread count, JuBiC emits a warning with the
calculation but preserves the requested solver and worker settings.

For `MIPparam`, only `threads_master` is relevant because the model is solved as a single MIP.

### Branching, bounds, and cut settings

For GBC, BlC, and BlCLag, select `branching_rule` as `"default"`,
`"linking_first"`, or `"fixed_linking_order"` in the solver parameter JSON.
These control the master branching priorities; they do not change the instance
or select a different formulation. `parallel_subsolvers` controls the number of
follower workers separately from MIP thread counts.

GBC additionally accepts `connector_s_bound` to set the ConnectorLP `s` bound.
`integer_obj` controls stabilization of integral objective-derived cut values.
The HNDP runner defaults this to `true` for GBC and BlC, with an explicit JSON
override available. Big-M selection belongs to the model configuration's
`big_m_mode`, using `"fixed_network_path"`, `"n_minus_one_most_expensive"`, or
`"fixed_network_path_current_cost"`.

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
- `integer_obj`: if enabled, objective-derived values that are within `10^-4`
  of an integer are rounded to the nearest integer. Nonnegative cut
  coefficients use the conservative coefficient-rounding rule; values that
  are not close to integers are retained.
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

Use with `build_hndp_blclag_instance(...)`. The runner selects this route with
`"model_type": "blclag"`.

Important fields:
- `runtime`
- `seed`
- `threads_master`
- `threads_sub_con`
- `integer_obj`: if enabled, objective-derived values that are within `10^-4`
  of an integer are rounded to the nearest integer. Nonnegative cut
  coefficients use the conservative coefficient-rounding rule; values that
  are not close to integers are retained.
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
  - `HNDP_BIGM_FIXED_NETWORK_PATH_CURRENT_COST`
- `subproblem_method`
  - `HNDP_SUBPROBLEM_MIP`
  - `HNDP_SUBPROBLEM_ASTAR`

Use with:
- `BLCparam`

## `build_hndp_blclag_instance`

```julia
build_hndp_blclag_instance(hndp, solver;
    big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    subproblem_method=HNDP_SUBPROBLEM_BLC_JUMP)
```

Builds a Lagrangian Benders-like master. Use `BlCLagparam` to solve it.
The persistent BlC-JuMP oracle supports the big-M options above. The alternative
`HNDP_SUBPROBLEM_MIBS` oracle uses MiBS; it does not use the big-M callback.
MIP and A* are not accepted as direct BlCLag oracle selections.

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
  - `HNDP_BIGM_FIXED_NETWORK_PATH_CURRENT_COST`
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
    use_decision_arc_dominance=true,
    fixed_path_bound_method=:labeling,
)
```

Purpose:
- Builds the path-based compact reformulation.

Options:
- `enumeration_time_limit`
  Global time budget for path-bound computation and path enumeration.
- `parallelize`
  If `true`, uses the parallel precomputation variant.
- `fixed_path_bound_method`
  Selects `:labeling` (default) or `:mip` for the fixed-network bound in
  sequential or parallel precomputation.
- `use_decision_arc_dominance`
  Enables the decision-arc dominance pruning rule during path enumeration.

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
    fallback_if_paths_exceed_flow_vars=false,
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
- `fallback_if_paths_exceed_flow_vars`
  If `true`, users can also fall back after successful enumeration when the
  path set is larger than the fallback arc-flow variable count.

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
    fallback_if_paths_exceed_flow_vars=false,
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
- `fallback_if_paths_exceed_flow_vars`
  If `true`, users can also remain as BlC subproblems after successful
  enumeration when the path set is larger than the fallback arc-flow variable
  count.

Use with:
- `BLCparam`

## Recommended Pairings

Common model / parameter combinations:

- BlCLag model:
  `build_hndp_blclag_instance(...)` + `BlCLagparam(...)`

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

- For single-layer competition generation, the current arc-mode split is:
  - `single_layer_k_competition`
    - `beta` scaling on sampled decision arcs
    - only sampled decision arcs contribute operator-side risk
  - `single_layer_k_decision_only`
    - no `beta` scaling
    - sampled arcs are only the controllable arcs
    - risk remains on all arcs and is positive in this mode

## Multimodal example

The strong-duality example is configured in [multimodal_example.json](multimodal_example.json)
alongside the HNDP implementation. From the repository root, run:

```sh
julia --project=. examples/HNDP/run_hndp_manifest.jl examples/HNDP/multimodal_example.json tmp_compare/hndp_multimodal_example
```

Add `--dry-run` to construct and validate the instances without solving.
The solver reformulation remains in `hndp_model_generation.jl`; the JSON file
only selects its options. See [NETWORK_GENERATION.md](NETWORK_GENERATION.md)
for generation parameters.
