# OSM-based HNDP instances

Application instances on this branch are built from cached OSM road geometry,
municipal boundaries, and GTFS schedules. The four layers represent walking,
car travel, transit, and bicycles. Bike-station transfer arcs are grouped and
subject to a station budget.

The `HNDPwC` container retains its graph, users, controllable arcs, construction
costs, optional resource bounds, decision groups, and decision budget. Generic
model builders continue to accept this container independently of input files.

Load `examples/HNDP/network_generation/hndp_osm_network_generation.jl`. Use `generate_osm_gtfs`
for one city, `generate_hndp_osm_networks` for a configured batch, or
`visit_hndp_osm_networks` to stream generation. `read_hndp_instance` restores a
saved JSON network with sparse user matrices.

[City Networks from OSM and GTFS](../../generated/hndp/osm_networks.md) lists the
six-city sources and describes the network construction, the arc costs and the
physical arc lengths. [Demand and Aggregation](../../generated/hndp/demand_and_aggregation.md)
describes the demand model, the aggregation and the selection of potential
stations. The small toy graphs are retained as regression fixtures only.


## Station decisions, construction cost, and where coupling is added

The OSM importer produces a generic `HNDPGeneratedNetwork` containing an
`HNDPwC`; it does not build a JuMP model or add solver constraints. For station
`s`, `decision_groups` stores `[entry_s, return_s]`, and `edgeA` contains both
arcs. `bike_station_construction_cost=K` assigns K to entry and zero to return.
A vector supplies one K per station, in station/group order. The common helper
`hndp_sharing_station_decisions` implements this convention. K defaults to 0 for the current OSM application.

The data flow is:

1. `visit_hndp_osm_networks` / `generate_hndp_osm_networks` calls
   `generate_osm_gtfs` and returns the grouped generic HNDP instance.
2. Optional JSON interchange (`_write_hndp_instance_json`,
   `read_hndp_instance`, format version 2) preserves both `decision_arc` flags,
   the entry/return construction coefficients, `decision_groups`, and budget.
   It preserves data; it does not add coupling constraints or silently repair
   old instances.
3. `run_hndp_experiments!` passes `generated_network.instance` through
   `_build_hndp_model_from_spec` to the chosen SD, BlC, or path model builder.
4. The model builder creates arc availability variables, fixes nondecision
   arcs, calls `_apply_hndp_decision_groups!` to add
   `x[entry_s] == x[return_s]`, and calls `_apply_hndp_decision_budget!` to add
   `sum_s x[entry_s] <= budget`. Construction cost is
   `sum_(a in edgeA) edge_price[a] * x[a]`, hence `sum_s K_s*x[entry_s]`.

There is no intervening multimodal conversion that infers these equalities
from layer indices. The shared model builders consume the explicit groups.
Their current implementation still creates variables for fixed arcs and adds
equalities fixing them to one; eliminating those variables is a separate change.

The same convention is used by `build_osm_full_od_benchmark` in
`osm_full_od_mip_benchmark.jl`, with keyword `station_construction_cost`.
The completed historical full-OD benchmark charged both directions and is kept
unchanged as a record of what was run. Future runs should use the corrected
builder and a fresh output directory.

Older generator instances that list only entry in `edgeA` but include return
in its group must be regenerated or repaired explicitly. The model builder now
rejects a grouped arc missing from `edgeA`, rather than fixing its return open
and consequently forcing the entire station open. No historical network
archives are rewritten automatically.

The full-OD benchmark now takes demand generated on the full fixed network
before graph aggregation. `prepare_osm_fixed_demand` computes class-specific
paths, kilometre bands and directed nearest-site maps; the benchmark builder
requires its `demand_root` and retains local/diagonal mass separately. See
[Demand and Aggregation](../../generated/hndp/demand_and_aggregation.md) for the exact definitions.
