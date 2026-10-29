# OSM-based HNDP instances

Application instances on this branch are built from cached OSM road geometry,
municipal boundaries, and GTFS schedules. The four layers represent walking,
car travel, transit, and bicycles. Bike-station transfer arcs are grouped and
subject to a station budget.

The `HNDPwC` container retains its graph, users, controllable arcs, construction
costs, optional resource bounds, decision groups, and decision budget. Generic
model builders continue to accept this container independently of input files.

Load `examples/HNDP/hndp_osm_network_generation.jl`. Use `generate_osm_gtfs`
for one city, `generate_hndp_osm_networks` for a configured batch, or
`visit_hndp_osm_networks` to stream generation. `read_hndp_instance` restores a
saved JSON network with sparse user matrices.

See [the generation guide](../../generated/hndp/generators.md) for configuration
and [the OSM/GTFS pipeline](../../generated/hndp/osm_pipeline.md) for buffered
extracts, boundary compression, transit conversion, and cost assumptions.
The small toy graphs are retained as regression fixtures only.
