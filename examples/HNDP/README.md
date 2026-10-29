# HNDP examples

Start Julia from the repository root with `julia --project=. --threads=2`.
The HNDP example APIs are loaded with `include`; they are not package exports.

- [Instance generators and construction examples](NETWORK_GENERATION.md):
  constrained shortest paths, layered/single-layer competition, multimodal bike
  stations, transit expansion, and saved networks.
- [Solver comparison and solve examples](HNDP_SOLVER_OVERVIEW.md): BlC, BlCLag,
  GBC, MiBS, strong duality, path enumeration, and hybrid formulations.
- [General generation configuration](network_generation_example.json).
- [Small multimodal strong-duality configuration](multimodal_example.json).
- [OSM/GTFS city-instance pipeline](AACHEN_OSM_GTFS_PIPELINE.md), retained on
  this branch for realistic road and transit networks.
- [Archived OSM/GTFS model comparisons](aggregated_network_model_results/README.md).

Run the two multimodal examples with a working Gurobi license:

```sh
julia --project=. --threads=2 examples/HNDP/run_hndp_manifest.jl examples/HNDP/multimodal_example.json tmp_compare/hndp_multimodal_example
```

Append `--dry-run` to generate and count jobs without solving. Custom runs can
also supply a directory containing `instances.json`, `models.json`, and
`params.json`. The solver implementations are in `hndp_model_generation.jl`;
the generic batch runner is `hndp_experiment_runner.jl`.

Generated output belongs under `tmp/` or `tmp_compare/`, or another directory
you choose.
