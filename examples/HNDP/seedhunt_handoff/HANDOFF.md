# HNDP Seed-Hunt Handoff

This folder is a compact handoff bundle for the HNDP parameter-seed runtime
tests that were started earlier and should now be continued on another PC.

## Read This First

If you are an AI model continuing the work, start with exactly these files:

1. [run_seedhunt.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/seedhunt_handoff/run_seedhunt.jl)
2. [hndp_experiment_runner.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/hndp_experiment_runner.jl)
3. [hndp_model_generation.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/hndp_model_generation.jl)

In most cases, `run_seedhunt.jl` plus `hndp_experiment_runner.jl` are enough to
understand what batch is being run and how to resume it.

## What This Batch Does

The seed-hunt batch benchmarks a single HNDP GBC configuration over multiple
instance-generation seeds and three user counts.

Current setup:

- Topology: `sioux_falls`
- Instance type: `constrained_shortest_path`
- User counts: `1200`, `1600`, `2400`
- Parameter seeds: `1:8`
- `length_constrained = true`
- `length_slack = 0.0`
- `two_stage = false`
- `user_parameter_mode = "per_user"`
- Construction / max arc parameters all set to `100`

Model setup:

- `model_type = "gbc"`
- `name = "gbc_mip_no_partial_ws"`
- `partial_decomposition = false`
- `subproblem_method = "mip"`
- `big_m_mode = "n_minus_one_most_expensive"`

Solver parameter setup:

- `mip_solver = "Gurobi"`
- `runtime = 60.0`
- `seed = 42`
- `threads_master = 2`
- `threads_sub_con = 2`
- `pareto = "OPT"`
- `warmstart = true`
- `write_run_logs = false`

## How To Run

From the repository root:

```bash
julia --project=. examples/HNDP/seedhunt_handoff/run_seedhunt.jl
```

This writes results into:

```text
examples/HNDP/seedhunt_handoff/seedhunt_runs
```

To continue on another machine and keep results somewhere else:

```bash
julia --project=. examples/HNDP/seedhunt_handoff/run_seedhunt.jl /path/to/output true
```

Arguments:

- first argument: `output_root`
- second argument: `resume` as `true` or `false`

Examples:

```bash
julia --project=. examples/HNDP/seedhunt_handoff/run_seedhunt.jl /tmp/hndp_seedhunt true
julia --project=. examples/HNDP/seedhunt_handoff/run_seedhunt.jl /tmp/hndp_seedhunt false
```

## Resume Behavior

The script regenerates the JSON config files under:

```text
<output_root>/configs
```

and then calls the streamed HNDP experiment runner with `resume=true` by
default. The runner resumes from:

```text
<output_root>/results/batch_summary.csv
```

So to continue on another PC, you mainly need:

- this folder
- the repo codebase
- the previous `output_root` folder if you want to resume from existing results

## What To Commit

To keep the git footprint small, commit only this handoff folder if the goal is
just to preserve the test recipe:

- [run_seedhunt.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/seedhunt_handoff/run_seedhunt.jl)
- [HANDOFF.md](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/seedhunt_handoff/HANDOFF.md)

Do **not** commit generated `output_root` data unless you explicitly want the
partial results in git as well.

## Notes On The Original Local State

There was an earlier scratch version under `tmp/hndp_seedhunt.jl` and
`tmp/hndp_seedhunt/`, but this handoff folder supersedes it. The scratch files
do not contain additional logic beyond what is captured in this bundle.
