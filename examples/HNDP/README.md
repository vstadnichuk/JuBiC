# Station location on multimodal OSM/GTFS networks

This folder holds the application of the paper "Bilevel Feasible Paths for the
Station Location Problem": the generator of the multimodal city networks, the demand
and aggregation procedure, and the station-based path model. The scripts to
reproduce the paper are in [station_location_paper/scripts/](station_location_paper/scripts/).

## Guides

| Guide | Content |
|---|---|
| [INSTANCE_GENERATOR.pdf](INSTANCE_GENERATOR.pdf) | Description of the instance generator as in the paper: arc costs, demand, maps of the six cities, comparison with observed modal splits |
| [OSM_NETWORKS.md](OSM_NETWORKS.md) | Sources, network construction, arc costs and physical arc lengths of the six city networks |
| [DEMAND_AND_AGGREGATION.md](DEMAND_AND_AGGREGATION.md) | Full-network demand, aggregation to 50–400 locations, potential stations and budgets, benchmark instances |
| [STATION_PATH_MODEL.md](STATION_PATH_MODEL.md) | Formulation, enumeration of bilevel feasible paths, experiment cases and the dominance audit |

## Pipeline

| Step | Script in `station_location_paper/scripts/` | Library code | Guide |
|---|---|---|---|
| Download and prepare sources | `download_sources.py`, `prepare_sources.py` | [sources/osm_city_sources.json](sources/osm_city_sources.json) | [Networks](OSM_NETWORKS.md#sources) |
| Build city networks | `generate_networks.jl` | [network_generation/](network_generation/) | [Networks](OSM_NETWORKS.md#network-construction) |
| Recover transit lengths | — | [network_generation/osm_transit_arc_lengths.py](network_generation/osm_transit_arc_lengths.py) | [Networks](OSM_NETWORKS.md#physical-arc-lengths) |
| Full-network demand and modal split | `full_network_city_profiles.jl` | [demand/](demand/) | [Demand](DEMAND_AND_AGGREGATION.md#full-network-demand) |
| Aggregate instances | `prepare_instances.jl` | [demand/](demand/) | [Demand](DEMAND_AND_AGGREGATION.md#aggregation) |
| Run experiments | `run_experiments.py`, `worker.jl`, `benchmark_factory.jl`, `audit_station.jl` | [models/hndp_station_path_model.jl](models/hndp_station_path_model.jl) | [Model](STATION_PATH_MODEL.md) |
| Replay designs on the full networks | `replay_designs_full_network.jl`, `resolve_full_network_ties.py` | [demand/](demand/) | — |
| Tables and figures | `analyze_*.py`, `figures_*.py`, `network_topology.py`, `figure_network_maps.py` | — | — |

## Folders

| Folder | Role in the paper |
|---|---|
| [network_generation/](network_generation/) | OSM/GTFS importer (`generate_osm_gtfs`), walking cleanup, road and transit length recovery |
| [demand/](demand/) | All-pairs reference paths, distance-band demand, aggregation and modal-split evaluation |
| [models/](models/) | Station-based path model (`hndp_station_path_model.jl`) and the HNDP instance types and path helpers it builds on |
| [experiments/](experiments/) | Benchmark loader and runner utilities that the experiment worker includes |
| [sources/](sources/) | Source registries: downloads and checksums, distance bands, car availability |
| [tests/](tests/) | Tests; `osm_station_path_tests.jl` covers the paper's model |
| [station_location_paper/](station_location_paper/) | Settings, requirements and reproduction scripts; results go to the git-ignored `study/` folder |

