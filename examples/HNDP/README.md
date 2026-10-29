# OSM-based HNDP examples

Run Julia from the repository root with `julia --project=. --threads=2`.

- [OSM/GTFS source-data pipeline](OSM_GTFS_PIPELINE.md).
- [OSM generation API and saved networks](NETWORK_GENERATION.md).
- [Solver/model parameters](HNDP_SOLVER_OVERVIEW.md).
- [OSM generation configuration](osm_network_generation_example.json).
- [OSM path/strong-duality experiment](osm_example.json).
- [Archived city model comparisons](aggregated_network_model_results/README.md).

With cached OSM, boundary, and GTFS inputs:

```sh
julia --project=. --threads=2 examples/HNDP/run_osm_hndp_manifest.jl examples/HNDP/osm_example.json tmp_compare/osm_example --dry-run
```

Remove `--dry-run` to optimize. `--no-resume` reruns completed jobs. Custom
experiments may use a directory containing `instances.json`, `models.json`, and
`params.json`; `saved_instances` reuses generated JSON networks.

`hndp_osm_network_generation.jl` constructs the application networks,
`hndp_model_generation.jl` builds formulations, and
`hndp_osm_experiment_runner.jl` runs model comparisons. Outputs belong under
`tmp/` or `tmp_compare/`; raw data belong under ignored `data/osm/`.

For an offline fixture test, run `julia --project=. examples/HNDP/osm_smoke.jl`.
