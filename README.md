# JuBiC: Paper Reproduction Branch

This branch is dedicated to the paper
[*A Catalog of Formulations for the Multi-Follower Discrete Bilevel Network Design Problem*](https://optimization-online.org/?p=35437).
Here, we provide the code for the formulations discussed in the paper. You can
also reproduce the computations reported in the paper.

If you want to work with JuBiC beyond reproducing this paper, we recommend
checking out the latest version of the `main` branch. The broader JuBiC
documentation starts at [`docs/src/index.md`](./docs/src/index.md), and the
HNDP example, i.e., the description of the existing algorithms including the
instance generation used in the paper, is documented under
[`docs/src/examples/hndp/motivation.md`](./docs/src/examples/hndp/motivation.md).

To build the rendered HTML documentation locally, run:

```powershell
julia --project=. docs/make.jl
```

---

## Getting Started

1. Install Julia and instantiate the project environment:

```powershell
julia --project=. -e "using Pkg; Pkg.instantiate()"
```

2. Run the unit tests:

```powershell
julia --project=. test/runtests.jl
```

3. Read the HNDP documentation:

- [HNDP Motivation](./docs/src/examples/hndp/motivation.md)
- [HNDP Instances](./docs/src/examples/hndp/instances.md)
- [HNDP Solver Models](./docs/src/examples/hndp/solvers.md)
- [HNDP Benchmark Pipeline](./docs/src/examples/hndp/benchmarks.md)

---

## Reproducing Results

Paper reproduction manifests are stored in
[`examples/HNDP/paper_manifests`](./examples/HNDP/paper_manifests). Each
experiment is defined by three JSON files:

- `instances.json`: topology and instance-generation grid,
- `models.json`: formulation or hybrid model choices,
- `params.json`: runtime, threads, solver, and logging settings.

The helper script
[`examples/HNDP/run_hndp_manifest.jl`](./examples/HNDP/run_hndp_manifest.jl)
loads such a manifest directory and calls the HNDP benchmark pipeline.

Before launching a long run, validate that the manifest loads:

```powershell
julia --project=. examples/HNDP/run_hndp_manifest.jl `
  examples/HNDP/paper_manifests/solvercomparison_layered_sioux_10min `
  tmp_compare/runs/paper_solvercomparison_layered_sioux_10min `
  --dry-run
```

To run the 10-minute Sioux Falls layered solver comparison:

```powershell
julia --project=. examples/HNDP/run_hndp_manifest.jl `
  examples/HNDP/paper_manifests/solvercomparison_layered_sioux_10min `
  tmp_compare/runs/paper_solvercomparison_layered_sioux_10min
```

This sweep compares path enumeration, strong-duality formulations, and
Benders-like cut formulations on layered Sioux Falls instances with a
10-minute time limit per solve.

To validate the 30-minute path and hybrid-SD manifest:

```powershell
julia --project=. examples/HNDP/run_hndp_manifest.jl `
  examples/HNDP/paper_manifests/path_hybrid_sd_layered_sioux_ema_30min `
  tmp_compare/runs/paper_path_hybrid_sd_layered_sioux_ema_30min `
  --dry-run
```

To run it:

```powershell
julia --project=. examples/HNDP/run_hndp_manifest.jl `
  examples/HNDP/paper_manifests/path_hybrid_sd_layered_sioux_ema_30min `
  tmp_compare/runs/paper_path_hybrid_sd_layered_sioux_ema_30min
```

This sweep compares complete path enumeration with the hybrid strong-duality
model. The hybrid model enumerates paths first and falls back to the
strong-duality formulation when enumeration fails or when the enumerated path
model would exceed the size of the arc-flow formulation.

The benchmark pipeline appends partial results to
`results/batch_summary.csv` after each completed instance and supports
restarting with the same output folder. Per-instance debug logs are disabled in
the supplied manifests, but the CSV and copied manifests are still written.

Hint: running all paper reproduction experiments can take roughly one week
on a typical workstation. The exact runtime depends strongly on CPU count,
Julia thread configuration, Gurobi version, and license settings.

The HNDP instances are randomly generated from the manifest seeds. Reruns should
therefore show trends comparable to the paper, but concrete runtimes can vary
with the generated instance structure.

---

## Disclaimer

JuBiC is still under active development. Results should be interpreted with
care, and we recommend independently verifying computational results whenever
possible.

If you encounter difficulties or have questions, contact:
Vladimir Stadnichuk, `vladimir.stadnichuk@uni-kassel.de`.

---

## Citation

If you use this branch for the HNDP formulation catalog, please cite:

> Stadnichuk and Koster, *A Catalog of Formulations for the Multi-Follower
> Discrete Bilevel Network Design Problem*, Optimization Online.
> [Link](https://optimization-online.org/?p=35437)

For the general JuBiC decomposition framework, see:

> Stadnichuk and Koster (2024), *Solving Multi-Follower Mixed-Integer Bilevel
> Problems with Binary Linking Variables*, Optimization Online.
> [Link](https://optimization-online.org/?p=28877)

---

## License

JuBiC is released as an open-source project. See the [LICENSE](./LICENSE) file
for details.
