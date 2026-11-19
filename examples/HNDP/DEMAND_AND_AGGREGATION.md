# Demand, aggregation and benchmark instances

This guide describes how the full city networks of [OSM_NETWORKS.md](OSM_NETWORKS.md)
become the aggregated instances of the experiments (Section 5 of the paper). 

## Two user classes

Demand is modeled for two classes: users with a car and users without one. The car
availability per city is listed in
[osm_car_ownership_sources.json](sources/osm_car_ownership_sources.json). Each class
has its own reference network on the full graph:

- the class without a car also loses every arc incident to the car layer.

## Full-network demand

1. For every ordered pair of distinct physical locations and each class, compute a
   minimum-cost path in the class's reference network. Record its cost and its
   physical length; on cost ties, the physically shorter path is taken.
2. Assign each pair to a distance band by that physical length. The band limits and
   masses per city are in
   [osm_distance_demand_sources.json](sources/osm_distance_demand_sources.json).
3. For class `c` with share `p[c]`, let `N[b,c]` be the number of pairs in band `b`.
   Each pair in band `b` gets demand `p[c] * q[b] / N[b,c]`, where `q[b]` is the
   band's mass. Total demand is one; a band with mass but no pair is an error.

`full_network_city_profiles.jl` builds the all-pairs cost and length caches
(`study/work/fixed_reference_cache/`, about 30 GB) and reports the users and the
modal split with all candidate stations open (Tables 2 and 6). The implementation is
[osm_fixed_network_demand.jl](demand/osm_fixed_network_demand.jl).

## Aggregation

`prepare_instances.jl` aggregates each city to `k` = 50, 75, 100, 200, 300 and 400
retained locations ([osm_aggregation_helpers.jl](demand/osm_aggregation_helpers.jl)).

1. **Metric.** Minimum journey cost between candidate stations on the full network,
   with the car and all candidate stations available; the cheaper direction makes it
   symmetric.
2. **Retained locations.** Choose `k` seeds among the candidate stations by
   farthest-point selection, fill clusters of equal size (±1) around them, and keep
   each cluster's medoid. Only existing candidate stations are retained.
3. **Demand mapping.** For each class, every original origin is mapped to the
   retained location it reaches most cheaply, and every original destination to the
   retained location that reaches it most cheaply. The demand of each original pair
   is added to its mapped pair.
4. **Same-location demand.** Pairs mapped to the same retained location stay on the
   diagonal. They get no follower, earn no sharing revenue and are not redistributed.
   Their mass is recorded in the metadata.
5. **Aggregated graph.** Each retained location has four layer nodes. Movement arcs
   between retained locations are minimum-cost paths within one mode on the full
   network (pure-mode closures), so a station decision cannot be hidden inside a
   mixed-mode arc ([osm_aggregation_mode_share_analysis.jl](demand/osm_aggregation_mode_share_analysis.jl)).
   Access and egress legs to the retained locations are not charged.

## Candidate stations and budgets

Potential stations of an instance are chosen among its retained locations as a
nested sequence: the network medoid first, then repeatedly the location farthest
from those already chosen. Every smaller pool is therefore a subset of every larger
one. The budget allows at most a given number of these stations to open.

| Experiment | Locations | Potential stations | Budget |
|---|---|---|---|
| Dominance audit | 50 | 20 | 10 |
| Aggregation | 50, 75, 100, 200, 300, 400 | 20 | 10 |
| Budget grid | 300 | 20, 80, 140 | 25%, 50%, 75% |
| Station scaling | 300 | 20, 40, 80, 100, …, 300 | 25% |

These values are set in [settings.json](station_location_paper/settings.json).

## Benchmark instances

`build_osm_hour_full_od_benchmark` (`station_location_paper/scripts/benchmark_factory.jl`)
turns an aggregated city into an `HNDPwC` instance when a case starts:

- two followers per ordered pair of distinct retained locations, one per class, each
  weighted by the aggregated demand of its class;
- the leader earns the 61-cent bike entry fee per unit of demand that enters the
  bike layer, and station construction costs zero;
- entry and return arcs of a station open together, within the budget.

## Outputs

- `study/inputs/aggregation/<city>/`: retained locations (`selected_nodes.csv`),
  selection metrics, aggregated graphs (`graph_<k>_edges.csv`) and potential-station
  sequences (`candidate_sets.csv`);
- `study/inputs/demand/<city>/`: class mappings (`mapping_<class>_<k>.csv`) and
  aggregated demand (`demand_<k>.csv`, exact and as floats, diagonal included);
- `study/data/<city>_aggregation_modal.csv`: modal split of each aggregated network
  with all potential stations open.
