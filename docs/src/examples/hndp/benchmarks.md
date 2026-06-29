# HNDP Benchmark Pipeline

The HNDP benchmark pipeline is designed for long experiment sweeps where partial
results must be preserved. It separates instance generation, model generation,
and solver parameters into three JSON files and streams one generated instance
at a time.

The central runner is `examples/HNDP/hndp_experiment_runner.jl`.

## Core Workflow

The function

```julia
run_hndp_experiments!(
    instance_config_path,
    model_config_path,
    param_config_path;
    output_root,
    resume = true,
)
```

performs the following steps:

1. load the network-generation config,
2. load the model config,
3. load the solver-parameter config,
4. generate one `HNDPGeneratedNetwork`,
5. build one JuBiC `Instance`,
6. solve it through `solve_instance!`,
7. append one row to `results/batch_summary.csv`,
8. continue with the next configuration.

Rows are appended after each instance solve. If a solve throws an exception, the
runner writes an error row instead of losing the whole batch.

## Output Folder

Each run writes an output root with the following structure:

```text
output_root/
  manifests/
    instances.json
    models.json
    params.json
    run_manifest.json
  results/
    batch_summary.csv
  runs/
    ...
```

The copied manifests are the main reproducibility record. They store the exact
network, model, and parameter configurations used for the sweep.

When `write_run_logs = false`, per-instance solver debug logs and model exports
are disabled. The CSV and launcher stdout/stderr are still preserved. Transient
run folders are created under the explicit `output_root`, not in an operating
system temp directory.

## Resume Behavior

With `resume = true`, the runner reads completed experiment ids from
`results/batch_summary.csv` and skips those configurations. This allows a sweep
to continue after interruption without restarting completed instances.

Experiment ids are generated from:

- generated instance name,
- model name,
- parameter name,
- a short hash of the full combination.

## Instance Configuration

The instance config is passed to the network generator. A typical layered
competition sweep looks like:

```json
{
  "parameter_seeds": [42, 43, 44],
  "instances": [
    {
      "name": "layered_path_screen",
      "instance_type": "competition",
      "topologies": ["layered_sioux_falls", "layered_ema"],
      "nusers": [100, 200, 500, 1000],
      "length_constrained": [false],
      "alpha": [1.0],
      "beta": [1.0, 0.8, 0.6],
      "construction_cost_min": 0,
      "construction_cost_max": 0,
      "availability_budget_fraction": [0.05, 0.1, 0.25, 0.5]
    }
  ]
}
```

The generator expands list-valued fields into a Cartesian product.

## Model Configuration

The model config selects HNDP builders. Example:

```json
{
  "models": [
    {
      "name": "path_dom",
      "model_type": "path",
      "enumeration_time_limit": 300.0,
      "parallelize": true,
      "use_decision_arc_dominance": true
    },
    {
      "name": "sd",
      "model_type": "sd",
      "big_m_mode": "fixed_network_path",
      "indicator_constraints": false,
      "bound_duals": true
    },
    {
      "name": "hybrid_sd_dom",
      "model_type": "hybrid_sd",
      "enumeration_time_limit": 300.0,
      "use_decision_arc_dominance": true,
      "fallback_if_paths_exceed_flow_vars": true
    }
  ]
}
```

Supported `model_type` values include:

- `blc`
- `blclag`
- `gbc`
- `mibs`
- `sd`
- `path`
- `hybrid_sd`
- `hybrid_blc`

## Parameter Configuration

The parameter config controls the solver runtime and logging behavior:

```json
{
  "write_run_logs": false,
  "params": [
    {
      "name": "gurobi_10min",
      "mip_solver": "Gurobi",
      "runtime": 600.0,
      "threads": 8
    }
  ]
}
```

The runner forwards the resolved solver name, output folder, logging flag, and
runtime to JuBiC's parameter constructors.

For path and hybrid models, enumeration time is counted against the total
runtime. The model builder records `enum_runtime`, and the remaining MIP runtime
is reduced accordingly.

## Reproducing a Sweep

Most benchmark scripts in `examples/HNDP/` are thin wrappers around
`run_hndp_experiments!`. They create temporary JSON config files, write them
into the benchmark output root, and call the runner.

Examples include:

- `hndp_path_sd_blc_layered_benchmark.jl`
- `hndp_path_sd_blc_layered_benchmark_u100_1000_with_dom.jl`
- `hndp_path_sd_budgeted_layered_ema_sioux_alpha1_benchmark.jl`
- `hndp_gbc_single_layer_k_sweep.jl`
- `hndp_gbc_single_layer_k_decision_only_sweep.jl`

A typical script call has the form:

```powershell
julia --project=. examples\HNDP\hndp_path_sd_blc_layered_benchmark_u100_1000_with_dom.jl `
    tmp_compare\runs\my_hndp_run `
    100,200,500,1000 `
    1.0,0.8,0.6 `
    600 `
    42
```

The exact positional arguments depend on the benchmark script. The reproducible
record is the output folder: keep `manifests/` and `results/batch_summary.csv`.

## Result Columns

`batch_summary.csv` combines:

- instance metadata with prefix `instance_`,
- model metadata with prefix `model_`,
- parameter metadata with prefix `param_`,
- JuBiC `RunStats` fields.

For path and hybrid models, additional useful columns include:

- `model_enum_runtime`,
- `model_path_counts`,
- `model_fallback_users`,
- `model_path_count_fallback_users`,
- `model_path_count_fallback_user_count`.

These columns make it possible to compare pure path, strong-duality, BlC, and
hybrid formulations from the same CSV file.
