# HNDP Benchmark Pipeline

The HNDP benchmark pipeline is designed for long experiment sweeps where partial
results must be preserved. It separates instance generation, model generation,
and solver parameters into three JSON files (or one combined file) and streams
one generated instance at a time.

The central runner is `examples/HNDP/hndp_osm_experiment_runner.jl`.

## Run the supplied OSM example

From the repository root, with a working Gurobi installation/license:

```sh
julia --project=. --threads=2 examples/HNDP/run_osm_hndp_manifest.jl examples/HNDP/osm_example.json tmp_compare/osm_example
```

Append `--dry-run` to construct the configured OSM instances and count jobs without solving;
append `--no-resume` to rerun completed jobs. The combined JSON contains
`instances`, `models`, and `params`. The launcher also accepts a directory
containing `instances.json`, `models.json`, and `params.json` for custom runs.
Paper-specific manifests and seed-hunt/reproduction batches are no longer part
of the generic example distribution.

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

The instance config is passed to the OSM generator. See the complete
`examples/HNDP/osm_network_generation_example.json` example. Each entry names
`osm_path`, `boundary_path`, and `gtfs_dir`, with `instance_type: "osm_gtfs"`.
User counts and parameter seeds can be expanded into a grid. For other cities,
supply their paths and geographic centre. `saved_instances` can instead name
previously generated JSON networks.

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
      "threads_master": 2,
      "threads_sub_con": 1,
      "parallel_subsolvers": 1,
      "parallel_separation": false
    }
  ]
}
```

The runner forwards the resolved solver name, output folder, logging flag, and
runtime to JuBiC's parameter constructors.

For path and hybrid models, enumeration time is counted against the total
runtime. The model builder records `enum_runtime`, and the remaining MIP runtime
is reduced accordingly.

## Result Columns

`batch_summary.csv` combines information from all three stages of the pipeline:
instance-generation parameters, model-generation metrics, and solver run
metrics. The concrete columns depend on the selected instance generator, model
builder, and solver configuration.

## Branching and reproducibility

For GBC, BlC, and BlCLag, `branching_rule` selects `"default"`, `"linking_first"`,
or `"fixed_linking_order"`. `parallel_subsolvers` controls follower worker count
separately from the per-model thread limits. Use fixed seeds and fixed thread
counts when comparing formulations. The [parameter reference](../../generated/hndp/formulations.md)
describes big-M callbacks, integer-cut stabilization, and `connector_s_bound`.
