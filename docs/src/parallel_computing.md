# Parallel Computing

JuBiC supports parallel work in two different places:

- inside one solve, where decomposition methods may solve follower-side tasks in parallel
- across benchmark jobs, where the benchmark runner may launch separate Julia processes

These two layers solve different problems and should be configured differently.

## In-Solver Parallelism

For `GBC`, `BlC`, and `BlCLag`, JuBiC can evaluate follower-side work in parallel through the parameter

```julia
parallel_separation = true
```

This means that different subproblems or connector LP tasks are processed concurrently. The goal is to reduce wall-clock time for multi-follower instances.

### Important Gurobi Restriction

When JuBiC uses Gurobi-backed subsolvers in parallel, each worker-side model must run single-threaded. Accordingly:

- `parallel_separation = true` requires `threads_sub_con = 1`
- JuBiC throws an exception if parallel separation is enabled while `threads_sub_con > 1`

The reason is that JuBiC creates separate Gurobi environments for concurrently solved worker models. This is safer than sharing one environment across active worker solves, but it also means that each worker model should use one Gurobi thread only.

The master solve can still use multiple threads through `threads_master`.

## Benchmark Parallelism

Long benchmark sweeps create many short-lived JuMP and Gurobi-backed models. In this setting, the main stability risk is not the solve itself, but repeated native-object cleanup inside one long Julia process.

To reduce this risk, the HNDP benchmark runner now supports

```julia
execution_mode = "subprocess_per_experiment"
```

and this is the default benchmark mode.

In this mode:

- the parent Julia process enumerates the pending benchmark jobs
- one child Julia process is launched for one `(instance, model, parameter)` experiment
- the child solves exactly that one job and appends one row to the summary CSV
- after that, the child process exits and native solver cleanup happens at process exit

This design is more robust for long Gurobi-heavy benchmark sweeps than solving many jobs inside one persistent Julia process.

## Resume Behavior

The subprocess benchmark mode still uses the same `experiment_id`-based resume logic as the in-process runner.

This means:

- completed rows already present in `results/batch_summary.csv` are skipped
- if one child process crashes, only the current experiment is lost
- restarting the same benchmark launcher continues from the existing CSV

## When To Use Which Mode

Use `execution_mode = "subprocess_per_experiment"` when:

- you run large benchmark sweeps
- you use many Gurobi-backed subsolves
- stability matters more than process-start overhead

Use `execution_mode = "in_process"` when:

- you want the older behavior for debugging
- you run only a small number of experiments
- process isolation is not needed

For typical benchmark runs with per-instance time limits in minutes, the extra process startup cost is usually negligible compared with the stability gain.
