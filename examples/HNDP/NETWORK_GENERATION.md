# OSM network generation

This branch constructs HNDP application instances from cached OpenStreetMap
roads, municipal boundaries, and GTFS transit schedules. The entry point is
[hndp_osm_network_generation.jl](hndp_osm_network_generation.jl); the road and
transit importer is [hndp_osm_gtfs_generation.jl](hndp_osm_gtfs_generation.jl).
Shared data containers are in `hndp_instances.jl`, and JSON/GEXF interchange is
in `hndp_osm_network_io.jl`. Solver/model builders remain independent of source
file formats.

## Prepare the source data

Follow [OSM_GTFS_PIPELINE.md](OSM_GTFS_PIPELINE.md) for buffered extracts,
boundary compression, GTFS conversion, and cost assumptions. Raw downloads
belong under ignored `data/osm/`. Supply a complete road extract covering the
municipal study area and its surrounding buffer for comparable city networks.

## Generate a city directly

Run from the repository root:

```julia
include("examples/HNDP/hndp_osm_network_generation.jl")
network = generate_osm_gtfs(
    osm_path="data/osm/aachen_buffer_highways.json",
    boundary_path="data/osm/aachen_boundary.json",
    gtfs_dir="data/osm/avv_gtfs",
    city_name="aachen", city_center=(50.775, 6.083),
    use_boundary_compression=true,
    nusers=100, seed=17, bike_station_count=20, bike_station_budget=10,
)
```

The returned `HNDPGeneratedNetwork` contains `.instance` (an `HNDPwC`) and
`.metadata`. Its four layers are walking, car, transit, and bicycle. Station
transfer arcs retain their decision groups and station budget.

For other cities, replace the three input paths and provide `city_name` and
`city_center`. `source_bbox` can record the source extract bounds. The Aachen
convenience functions remain available for existing OSM benchmark scripts.

## Batch configurations

The complete generation configuration is
[osm_network_generation_example.json](osm_network_generation_example.json).
Every source entry uses `instance_type: "osm_gtfs"` and the three required path
fields. Input paths are relative to the working directory. `nusers` and
`parameter_seeds` accept scalars or arrays of integers; all other options pass
through as importer keyword arguments. Boundary compression defaults to true
in the batch API. Other source-instance types are rejected explicitly.

```julia
networks = generate_hndp_osm_networks(
    "examples/HNDP/osm_network_generation_example.json")
visit_hndp_osm_networks(
    "examples/HNDP/osm_network_generation_example.json",
    network -> println(network.name),
)
```

The visitor holds only one generated network at a time. `instance_output`
optionally writes each network before visiting it. Set `folder` and choose
`graph_format` as `json`, `gexf`, or `both`.

## Save and reuse networks

```julia
write_hndp_generated_networks([network], Dict{String,Any}(
    "folder" => "tmp_compare/osm_instances", "graph_format" => "json"))
restored = read_hndp_instance(
    "tmp_compare/osm_instances/" * network.name * "/instance.json")
```

`instance.json` preserves graph arcs, user costs and objectives, decision
groups, and the decision budget. The reader reconstructs sparse user matrices
so city-sized networks do not require dense quadratic allocations. The
shortest-weight matrix is a derived cache and is rebuilt only when needed.

A configuration with `saved_instances` can supply one JSON path or an array of
paths instead of `instances`. It works in both the batch API and experiment
runner; saved paths are deduplicated.

## Evaluate the models

```sh
julia --project=. --threads=2 examples/HNDP/run_osm_hndp_manifest.jl examples/HNDP/osm_example.json tmp_compare/osm_example --dry-run
```

Remove `--dry-run` to compare path enumeration and strong duality. The example
requires cached input data and a working Gurobi license for optimization. A
directory of `instances.json`, `models.json`, and `params.json` is also accepted.
See [HNDP_SOLVER_OVERVIEW.md](HNDP_SOLVER_OVERVIEW.md) for model parameters and
[aggregated_network_model_results/README.md](aggregated_network_model_results/README.md)
for the retained city benchmark results and reproduction scripts.

## Offline smoke test

```sh
julia --project=. examples/HNDP/osm_smoke.jl
```

This builds a tiny OSM/GTFS fixture locally and checks the renamed generation
API, station groups, source filtering, and saved-instance round trips without
network downloads or an optimization license. Small in-memory toy HNDP graphs
remain solely as solver regression fixtures.
