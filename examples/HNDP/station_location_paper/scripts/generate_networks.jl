# Generate the full multimodal city networks of the paper from the prepared OSM/GTFS sources.
#
#   julia --project=../../.. --threads=8 scripts/generate_networks.jl [city ...]
#
# Input  (from prepare_sources.py): <study>/source_data/sources/<city>/{roads_buffered.json,boundary.json,gtfs/}
# Output: <study>/inputs/city_networks_walk_connected/graphs/<city>/{edges.csv,metadata.json,physical_node_mapping.csv}
#         <study>/inputs/city_networks_walk_connected/speed_overrides/<city>.csv  (car-layer arc costs)
#         <study>/inputs/geometry/road_lengths/<city>.csv                        (physical OSM arc lengths)
#         <study>/source_data/graphs/<city>/metadata.json   (node coordinates before the walking cleanup,
#                                                           used by osm_transit_arc_lengths.py)
# All cost parameters are the defaults of generate_osm_gtfs (see examples/HNDP/INSTANCE_GENERATOR.pdf).
using CSV, JSON, Graphs, SparseArrays, SHA, Dates
include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(SRC, "network_generation", "hndp_osm_network_generation.jl"))
include(joinpath(SRC, "network_generation", "osm_road_arc_lengths.jl"))

const CONFIG = JSON.parsefile(SETTINGS_FILE)
const REGISTRY = Dict(c["city"] => c for c in JSON.parsefile(joinpath(SRC, "sources", "osm_city_sources.json"))["cities"])
const SOURCES = joinpath(B, "source_data")
const NETWORKS = joinpath(B, "inputs", "city_networks_walk_connected")

"""Node coordinates after road compression and car-component restriction, before the walking cleanup."""
function pre_cleanup_coords(input)
    coords, car, bike, walk = _read_osm_ways(joinpath(input, "roads_buffered.json"))
    coords, car, bike, walk = _compress_external_car_edges(coords, car, bike, walk, joinpath(input, "boundary.json"))
    coords, car, bike, walk = _restrict_to_car_component(coords, car, bike, walk)
    return coords
end

function generate_city(city)
    folder = joinpath(NETWORKS, "graphs", city)
    isfile(joinpath(folder, "COMPLETE.json")) && (println("CACHED ", city); return)
    mkpath(folder); mkpath(joinpath(NETWORKS, "speed_overrides"))
    input = joinpath(SOURCES, "sources", city)
    spec = REGISTRY[city]
    println("GENERATING ", city); flush(stdout); timer = time()
    generated = generate_osm_gtfs(osm_path=joinpath(input, "roads_buffered.json"),
        boundary_path=joinpath(input, "boundary.json"), gtfs_dir=joinpath(input, "gtfs"),
        use_boundary_compression=true, city_name=city, city_center=Tuple(Float64.(spec["centre"])),
        nusers=1, seed=20261028, bike_station_count=20, bike_station_budget=10)
    h = generated.instance; m = generated.metadata; n = nv(h.mygraph) ÷ 4
    m["all_candidate_stations"] = sort!(Int.(m["all_candidate_stations"]))
    m["buffer_km"] = spec["buffer_km"]; m["boundary_relation_id"] = spec["relation_id"]
    m["source_bbox"] = spec["bbox"]; m["generation_seconds"] = time() - timer
    m["generator_sha256"] = bytes2hex(sha256(read(joinpath(SRC, "network_generation", "hndp_osm_gtfs_generation.jl"))))
    m["source_commit"] = "JuBiC generate_osm_gtfs, generator sha256 " * m["generator_sha256"][1:12]
    m["generated_at"] = string(now())
    groups = Dict(arc => i for (i, arcs) in enumerate(h.decision_groups) for arc in arcs)
    costs = h.users[1].mcost; risks = h.users[1].mrisk
    CSV.write(joinpath(folder, "edges.csv"), ((source=src(e), target=dst(e), cost_cents=costs[src(e), dst(e)],
        leader_objective_coefficient=risks[src(e), dst(e)], decision_group=get(groups, (src(e), dst(e)), 0))
        for e in edges(h.mygraph)))
    # Car-layer costs as an explicit overlay; downstream readers apply it to the car layer.
    CSV.write(joinpath(NETWORKS, "speed_overrides", city * ".csv"), ((source=src(e), target=dst(e),
        cost_cents=costs[src(e), dst(e)]) for e in edges(h.mygraph) if n < src(e) <= 2n && n < dst(e) <= 2n))
    retained = Int.(m["walking_retained_input_ids"])
    CSV.write(joinpath(folder, "physical_node_mapping.csv"),
        ((original_node=old, new_node=new) for (new, old) in enumerate(retained)))
    open(joinpath(folder, "metadata.json"), "w") do io; JSON.print(io, m, 2); end
    # Coordinates before the walking cleanup (GTFS stops are snapped to these nodes).
    pre = joinpath(SOURCES, "graphs", city); mkpath(pre)
    open(joinpath(pre, "metadata.json"), "w") do io
        JSON.print(io, Dict("city" => city, "coords" => pre_cleanup_coords(input)), 2)
    end
    recover_osm_road_lengths(city; raw_root=SOURCES, network_root=NETWORKS,
        output_root=joinpath(B, "inputs", "geometry", "road_lengths"))
    summary = Dict("city" => city, "physical_nodes" => n, "eligible_stations" => length(m["all_candidate_stations"]),
        "arcs" => ne(h.mygraph), "seconds" => time() - timer)
    open(joinpath(folder, "COMPLETE.json"), "w") do io; JSON.print(io, summary, 2); end
    println(JSON.json(summary)); flush(stdout)
end

for city in (isempty(ARGS) ? String.(CONFIG["cities"]) : ARGS)
    generate_city(city)
    GC.gc()
end
