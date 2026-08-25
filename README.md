# JuBiC — paper reproduction branch

This branch is associated with the paper “Automated Benders-like Cut Generation and its Application to the Bilevel Network Design Problem”.

## EMA benchmark reproduction

The published benchmark is the EMA topology experiment. Reproduce the complete BlC/BlCLag run with:

```powershell
julia --project=. benchmark_run/fixed_order_big_runs/reproduce_ema_big_runs.jl
```

The script runs EMA negative-profit k-decision instances with 43 decision arcs, 50–200 users in steps of 25, alpha values 0.01/0.05/0.1, seeds 2–4, both length settings, fixed linking-order branching, 8 Gurobi threads, sequential separation, and a 10-minute limit. It runs all three BlC Big-M modes and BlCLag with warm start enabled.

The script prints a Big-M summary after each solver/instance combination and writes the full console stream to `console.log`. Raw output is written under `benchmark_run/runs/`; set `JUBIC_REPRO_OUTPUT` to choose another raw-output location.

## Published results

The reviewer-facing results are in [`benchmark_run/fixed_order_big_runs`](./benchmark_run/fixed_order_big_runs):

- [`all_results.csv`](./benchmark_run/fixed_order_big_runs/all_results.csv) contains 336 results: 252 BlC and 84 BlCLag warm-start runs.
- [`big_m_coefficient_summary.txt`](./benchmark_run/fixed_order_big_runs/big_m_coefficient_summary.txt) contains Big-M coefficient statistics aggregated over all EMA runs.
- [`reproduce_ema_big_runs.jl`](./benchmark_run/fixed_order_big_runs/reproduce_ema_big_runs.jl) reproduces both solver suites.

## Setup and tests

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. test/runtests.jl
```

The supplied configuration uses Gurobi 13.0.1 when installed at `C:\gurobi1301\win64`. Runtime depends on the CPU, Gurobi installation, license, and solver settings.

## Citation

For the general JuBiC decomposition framework, see:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online.

JuBiC is released under the license in [LICENSE](./LICENSE).
