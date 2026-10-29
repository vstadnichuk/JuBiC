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

## OSM network experiments

This branch evaluates the path-based models on realistic OSM/GTFS networks.
Start with the [source-data pipeline](examples/HNDP/OSM_GTFS_PIPELINE.md) and
the [OSM generation guide](examples/HNDP/NETWORK_GENERATION.md). The entry point
is `examples/HNDP/hndp_osm_network_generation.jl`.

With the source data cached under `data/osm/`, validate the supplied experiment:

```powershell
julia --project=. --threads=2 examples/HNDP/run_osm_hndp_manifest.jl `
  examples/HNDP/osm_example.json tmp_compare/osm_example --dry-run
```

Remove `--dry-run` to compare path enumeration and strong duality. The runner
supports saved JSON instances, resume, model configurations, and solver limits.
The [archived city comparisons](examples/HNDP/aggregated_network_model_results/README.md)
retain the realistic-network results and reproduction scripts.

For a small offline generation and round-trip test:

```powershell
julia --project=. examples/HNDP/osm_smoke.jl
```

The source-format-specific experiments removed from this branch remain on
`main`. Solver implementations and their small in-memory test fixtures remain
available here.

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
