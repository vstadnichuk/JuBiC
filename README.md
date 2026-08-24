# JuBiC — paper reproduction branch

This branch is associated with the paper "Automated Benders-like Cut Generation and its Application to the Bilevel Network Design Problem".

## Reproduction entry point

The complete fixed-order benchmark is reproduced by one script:

```powershell
julia --project=. examples/HNDP/reproduce_fixed_order_big_runs.jl
```

The script regenerates the Sioux Falls layered instances and runs:

- BlC with n−1, fixed-path, and fixed-path/current-cost Big-M modes;
- BlCLag with n−1 Big-M, with warm start and cold start;
- fixed linking-first branching order;
- sequential separation, 8 Gurobi threads for master and subproblems, and a 10-minute limit.

After each instance/solver combination, the script prints a compact Big-M coefficient summary. It writes a fresh result tree under `benchmark_run/runs/`. To choose another output location, set `JUBIC_REPRO_OUTPUT` before launching the script:

```powershell
$env:JUBIC_REPRO_OUTPUT = "tmp_compare/runs/fixed_order_reproduction"
julia --project=. examples/HNDP/reproduce_fixed_order_big_runs.jl
```

## Published benchmark results

The directly inspectable results are in [`benchmark_run/fixed_order_big_runs`](./benchmark_run/fixed_order_big_runs):

- [`all_results.csv`](./benchmark_run/fixed_order_big_runs/all_results.csv) contains all 300 solver results;
- [`big_m_coefficient_summary.txt`](./benchmark_run/fixed_order_big_runs/big_m_coefficient_summary.txt) reports Big-M coefficient ranges aggregated over all instances and for a common instance solved optimally by all five configurations.

## Setup and tests

Instantiate the Julia environment and run the tests with:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
julia --project=. test/runtests.jl
```

The HNDP documentation is available at:

- [HNDP motivation](./docs/src/examples/hndp/motivation.md)
- [HNDP instances](./docs/src/examples/hndp/instances.md)
- [HNDP solver models](./docs/src/examples/hndp/solvers.md)
- [HNDP benchmark pipeline](./docs/src/examples/hndp/benchmarks.md)

The reproduction sweep is computationally intensive. Runtime depends on the CPU, Gurobi version, license, and solver settings. The supplied configuration uses Gurobi 13.0.1 when installed at `C:\gurobi1301\win64`.

## Citation

For the general JuBiC decomposition framework, see:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel Problems with Binary Linking Variables*, Optimization Online.

JuBiC is released under the license in [LICENSE](./LICENSE).
