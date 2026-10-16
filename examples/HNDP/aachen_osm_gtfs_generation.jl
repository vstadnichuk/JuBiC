"""Build a JuBiC multimodal instance from cached Aachen OSM and AVV GTFS data.

This is intentionally a small, reproducible importer rather than a full traffic
assignment pipeline. OSM ways are converted to directed arcs, bicycle-only OSM
ways are retained only in the bike layer, and GTFS stop sequences become fixed
transit services. Raw downloads belong in the ignored `data/osm/` directory.
"""

using CSV
using Graphs
using JSON
using Random
using SparseArrays

include("hndp_network_generation.jl")

const AACHEN_BBOX = (50.75, 6.05, 50.82, 6.15)
const AACHEN_BUFFER_BBOX = (50.68, 5.90, 50.88, 6.30)
const AACHEN_CENTER = (50.775, 6.083)

_s(x) = x === missing ? "" : string(x)
_f(x) = x === missing || x === nothing || isempty(_s(x)) ? nothing :
    x isa Number ? Float64(x) : try parse(Float64, strip(_s(x))) catch; nothing end

_time_cost_cents(seconds, value_of_time_eur_per_hour) =
    max(1, round(Int, 100 * seconds * value_of_time_eur_per_hour / 3600))
_mode_cost(seconds, value_of_time_eur_per_hour, generalized_cost) = generalized_cost ?
    _time_cost_cents(seconds, value_of_time_eur_per_hour) : max(1, round(Int, seconds))

function _haversine_km(lat1, lon1, lat2, lon2)
    r = 6371.0088
    p = π / 180
    a = sin((lat2-lat1)*p/2)^2 + cos(lat1*p)*cos(lat2*p)*sin((lon2-lon1)*p/2)^2
    return 2r * asin(min(1.0, sqrt(a)))
end

function _osm_way_kind(tags)
    highway = get(tags, "highway", "")
    bike_only = highway in ("cycleway", "path", "track", "footway", "pedestrian") &&
        get(tags, "bicycle", "yes") != "no"
    explicit_bike = get(tags, "bicycle", "") in ("yes", "designated", "permissive")
    car_usable = highway in ("motorway", "trunk", "primary", "secondary", "tertiary",
        "unclassified", "residential", "living_street", "service", "motorway_link",
        "trunk_link", "primary_link", "secondary_link", "tertiary_link") &&
        get(tags, "motor_vehicle", get(tags, "access", "yes")) ∉ ("no", "private")
    high_speed = highway in ("motorway", "motorway_link", "trunk", "trunk_link")
    bike_usable = highway != "steps" && (bike_only || explicit_bike || (car_usable && !high_speed))
    walk_usable = highway != "steps" && !high_speed &&
        get(tags, "access", "yes") ∉ ("no", "private")
    return car_usable, bike_usable, bike_only, walk_usable
end

const CAR_CLASS_SPEED_KMH = Dict(
    "motorway" => 90, "motorway_link" => 70, "trunk" => 65, "trunk_link" => 55,
    "primary" => 45, "primary_link" => 35, "secondary" => 40, "secondary_link" => 32,
    "tertiary" => 35, "tertiary_link" => 30, "unclassified" => 30,
    "residential" => 25, "living_street" => 15, "service" => 20,
)

const BIKE_CLASS_SPEED_KMH = Dict(
    "cycleway" => 20, "track" => 18, "path" => 15, "footway" => 12,
    "pedestrian" => 12, "living_street" => 18, "residential" => 17,
    "service" => 15, "primary" => 17, "secondary" => 17, "tertiary" => 18,
    "unclassified" => 16,
)

function _tagged_speed_kmh(tags)
    value = get(tags, "maxspeed", "")
    m = match(r"(\\d+(?:\\.\\d+)?)", value)
    m === nothing && return nothing
    speed = parse(Float64, m.captures[1])
    occursin(r"mph", lowercase(value)) && (speed *= 1.609344)
    return speed
end

function _car_operating_speed(tags, kind, fallback)
    legal = _tagged_speed_kmh(tags)
    base = get(CAR_CLASS_SPEED_KMH, kind, fallback)
    # Legal speed limits are reduced to approximate operating/free-flow speed.
    return round(Int, clamp(something(legal, Float64(base)) * (legal === nothing ? 1.0 : 0.60), 10.0, 110.0))
end

function _bike_operating_speed(tags, kind, fallback)
    base = get(BIKE_CLASS_SPEED_KMH, kind, fallback)
    dedicated = kind in ("cycleway", "track") ||
        any(startswith(key, "cycleway") for key in keys(tags))
    return dedicated ? max(base, 19) : base
end

function _read_osm_ways(path::String; car_default_speed_kmh=35, bike_default_speed_kmh=16)
    raw = JSON.parsefile(path)
    points = Dict{Tuple{Float64,Float64},Int}()
    coords = Tuple{Float64,Float64}[]
    signal_points = Set{Tuple{Float64,Float64}}()
    for element in raw["elements"]
        get(element, "type", "") == "node" || continue
        tags = Dict{String,String}(String(k) => _s(v) for (k,v) in get(element, "tags", Dict()))
        get(tags, "highway", "") == "traffic_signals" || continue
        push!(signal_points, (round(Float64(element["lat"]), digits=7),
            round(Float64(element["lon"]), digits=7)))
    end
    car_edges = Dict{Tuple{Int,Int},NamedTuple}()
    bike_edges = Dict{Tuple{Int,Int},NamedTuple}()
    walk_edges = Dict{Tuple{Int,Int},NamedTuple}()
    function point_id(lat, lon)
        key = (round(lat, digits=7), round(lon, digits=7))
        get!(points, key) do
            push!(coords, key)
            length(coords)
        end
    end
    for element in raw["elements"]
        get(element, "type", "") == "way" || continue
        haskey(element, "geometry") || continue
        tags = Dict{String,String}(String(k) => _s(v) for (k,v) in get(element, "tags", Dict()))
        car_ok, bike_ok, bike_only, walk_ok = _osm_way_kind(tags)
        (car_ok || bike_ok || walk_ok) || continue
        geom = element["geometry"]
        length(geom) >= 2 || continue
        # Intermediate OSM geometry points describe shape, not intersections.
        # Contract them into the way length; retaining only way endpoints is
        # essential for keeping city extracts usable by JuBiC's layered model.
        ids = [point_id(Float64(geom[1]["lat"]), Float64(geom[1]["lon"])),
            point_id(Float64(geom[end]["lat"]), Float64(geom[end]["lon"]))]
        oneway = get(tags, "oneway", "no") in ("yes", "1", "true")
        cap = something(_f(get(tags, "lanes", missing)), 1.0)
        kind = get(tags, "highway", "unclassified")
        speed = _car_operating_speed(tags, kind, car_default_speed_kmh)
        for k in 1:length(ids)-1
            i, j = ids[k], ids[k+1]
            d = 0.0
            for q in 1:length(geom)-1
                d += _haversine_km(Float64(geom[q]["lat"]), Float64(geom[q]["lon"]),
                    Float64(geom[q+1]["lat"]), Float64(geom[q+1]["lon"]))
            end
            d > 0 || continue
            signal = coords[i] in signal_points || coords[j] in signal_points
            attr = (distance=d, capacity=cap, speed=speed, kind=kind, bike_only=bike_only, signal=signal)
            if car_ok
                car_edges[(i,j)] = attr
                oneway || (car_edges[(j,i)] = attr)
            end
            if bike_ok
                bike_attr = (distance=d, capacity=cap, speed=_bike_operating_speed(tags, kind, bike_default_speed_kmh), kind=kind, bike_only=bike_only)
                bike_edges[(i,j)] = bike_attr
                # Unless bicycle-specific tagging says otherwise, allow bikes to
                # use both directions even when motor traffic is one-way.  This
                # is the intended abstraction for the current bike layer: it
                # represents bicycle access to a road, not strict motor-vehicle
                # lane directionality.  Explicit oneway:bicycle=yes remains
                # respected.
                bicycle_oneway = get(tags, "oneway:bicycle", "") in ("yes", "1", "true")
                !bicycle_oneway && (bike_edges[(j,i)] = bike_attr)
            end
            if walk_ok
                walk_attr = (distance=d, capacity=cap, speed=5, kind=kind, bike_only=false)
                walk_edges[(i,j)] = walk_attr
                walk_edges[(j,i)] = walk_attr
            end
        end
    end
    return coords, car_edges, bike_edges, walk_edges
end

function _nearest_node(coords, lat, lon)
    best, best_d = 1, Inf
    for i in eachindex(coords)
        d = _haversine_km(lat, lon, coords[i][1], coords[i][2])
        d < best_d && ((best, best_d) = (i, d))
    end
    return best, best_d
end

function _gtfs_time_seconds(value)
    text = strip(_s(value)); isempty(text) && return -1
    parts = split(text, ':'); length(parts) == 3 || return -1
    try
        return 3600 * parse(Int, parts[1]) + 60 * parse(Int, parts[2]) + parse(Int, parts[3])
    catch
        return -1
    end
end

function _integer_median(values::Vector{Int})
    isempty(values) && return -1
    sorted = sort(values)
    return sorted[cld(length(sorted), 2)]
end

function _restrict_to_car_component(coords, car_edges, bike_edges, walk_edges)
    adjacency = [Int[] for _ in coords]
    for (i,j) in keys(car_edges)
        push!(adjacency[i], j); push!(adjacency[j], i)
    end
    seen = falses(length(coords)); largest = Int[]
    for start in eachindex(coords)
        seen[start] && continue
        queue = [start]; seen[start] = true; component = Int[]
        while !isempty(queue)
            v = popfirst!(queue); push!(component, v)
            for w in adjacency[v]
                if !seen[w]; seen[w] = true; push!(queue, w); end
            end
        end
        length(component) > length(largest) && (largest = component)
    end
    remap = Dict(old => new for (new, old) in enumerate(largest))
    newcoords = coords[largest]
    new_car = Dict{Tuple{Int,Int},Any}()
    new_bike = Dict{Tuple{Int,Int},Any}()
    new_walk = Dict{Tuple{Int,Int},Any}()
    for (e,a) in car_edges
        haskey(remap,e[1]) && haskey(remap,e[2]) && (new_car[(remap[e[1]],remap[e[2]])] = a)
    end
    for (e,a) in bike_edges
        haskey(remap,e[1]) && haskey(remap,e[2]) && (new_bike[(remap[e[1]],remap[e[2]])] = a)
    end
    for (e,a) in walk_edges
        haskey(remap,e[1]) && haskey(remap,e[2]) && (new_walk[(remap[e[1]],remap[e[2]])] = a)
    end
    return newcoords, new_car, new_bike, new_walk
end

function _read_boundary_rings(path::String)
    relation = first(JSON.parsefile(path)["elements"])
    segments = Vector{Vector{Tuple{Float64,Float64}}}()
    for member in get(relation, "members", Any[])
        get(member, "role", "") == "outer" || continue
        haskey(member, "geometry") || continue
        push!(segments, [(Float64(p["lat"]), Float64(p["lon"])) for p in member["geometry"]])
    end
    rings = Vector{Vector{Tuple{Float64,Float64}}}()
    while !isempty(segments)
        ring = popfirst!(segments)
        while !isempty(segments)
            last_point = last(ring); match = findfirst(s -> s[1] == last_point || s[end] == last_point, segments)
            match === nothing && break
            segment = splice!(segments, match)
            segment[1] == last_point || reverse!(segment)
            append!(ring, segment[2:end])
            ring[end] == ring[1] && break
        end
        push!(rings, ring)
    end
    return rings
end

function _point_in_ring(lat, lon, ring)
    inside = false
    j = length(ring)
    for i in eachindex(ring)
        lati, loni = ring[i]; latj, lonj = ring[j]
        crosses = ((loni > lon) != (lonj > lon)) &&
            (lat < (latj-lati) * (lon-loni) / (lonj-loni + eps()) + lati)
        crosses && (inside = !inside)
        j = i
    end
    return inside
end

function _compress_external_car_edges(coords, car_edges, bike_edges, walk_edges, boundary_path)
    rings = _read_boundary_rings(boundary_path)
    inside = [any(_point_in_ring(coords[v][1], coords[v][2], ring) for ring in rings) for v in eachindex(coords)]
    boundary_nodes = [v for v in eachindex(coords) if inside[v] &&
        any((i == v && !inside[j]) || (j == v && !inside[i]) for (i,j) in keys(car_edges))]
    result = Dict{Tuple{Int,Int},Any}()
    for (e,a) in car_edges
        inside[e[1]] && inside[e[2]] && (result[e] = a)
    end
    outgoing = [Tuple{Tuple{Int,Int},Any}[] for _ in eachindex(coords)]
    for (edge, attr) in car_edges
        push!(outgoing[edge[1]], (edge, attr))
    end
    for source in boundary_nodes
        distances = fill(Inf, length(coords)); distances[source] = 0.0
        queue = [(0.0, source)]
        while !isempty(queue)
            sort!(queue, by=first); current_distance, v = popfirst!(queue)
            current_distance > distances[v] && continue
            for (edge, attr) in outgoing[v]
                w = edge[2]
                w == source && continue
                if inside[w]
                    current_distance > 0 && !haskey(result, (source,w)) &&
                        (result[(source,w)] = (distance=current_distance + attr.distance, capacity=attr.capacity, speed=attr.speed, kind="external_shortcut", bike_only=false))
                else
                    candidate = current_distance + attr.distance
                    candidate < distances[w] && (distances[w] = candidate; push!(queue, (candidate,w)))
                end
            end
        end
    end
    filtered_bike = Dict(e=>a for (e,a) in bike_edges if inside[e[1]] && inside[e[2]])
    filtered_walk = Dict(e=>a for (e,a) in walk_edges if inside[e[1]] && inside[e[2]])
    # Keep the original indexing until `_restrict_to_car_component` performs
    # the final compact remapping.
    return coords, result, filtered_bike, filtered_walk
end

function _read_gtfs_services(zip_dir::String, coords, max_stop_distance_km)
    stops = Dict{String,Tuple{Int,Float64} }()
    for row in CSV.File(joinpath(zip_dir, "stops.txt"), silencewarnings=true)
        lat, lon = _f(getproperty(row, :stop_lat)), _f(getproperty(row, :stop_lon))
        if lat === nothing || lon === nothing
            continue
        end
        node, d = _nearest_node(coords, lat, lon)
        if d > max_stop_distance_km
            continue
        end
        stops[_s(getproperty(row, :stop_id))] = (node, d)
    end
    trips = Dict{String,Vector{NTuple{4,Int}}}()
    departures = Dict{Int,Vector{Int}}()
    for row in CSV.File(joinpath(zip_dir, "stop_times.txt"), silencewarnings=true)
        sid = _s(getproperty(row, :stop_id))
        tid = _s(getproperty(row, :trip_id))
        seq = something(_f(getproperty(row, :stop_sequence)), 0.0)
        arrival = _gtfs_time_seconds(getproperty(row, :arrival_time))
        departure = _gtfs_time_seconds(getproperty(row, :departure_time))
        # Keep unmatched stops as barriers in the original sequence.  We do
        # not want to connect two in-area stops after silently skipping an
        # out-of-area stop between them.
            node = haskey(stops, sid) ? stops[sid][1] : 0
            push!(get!(trips, tid, NTuple{4,Int}[]),
                (Int(round(seq)), node, arrival, departure))
            node > 0 && departure >= 0 && push!(get!(departures, node, Int[]), departure)
    end
    services = Set{Tuple{Int,Int}}()
    observed_times = Dict{Tuple{Int,Int},Vector{Int}}()
    for sequence in values(trips)
        sort!(sequence, by=first)
        last_matched = nothing
        skipped = 0
        for record in sequence
            if record[2] > 0
                if last_matched !== nothing && last_matched[2] != record[2]
                    pair = (last_matched[2], record[2])
                    push!(services, pair)
                    start_time = last_matched[4] >= 0 ? last_matched[4] : last_matched[3]
                    end_time = record[3] >= 0 ? record[3] : record[4]
                    duration = end_time - start_time
                    duration > 0 && push!(get!(observed_times, pair, Int[]), duration)
                end
                last_matched = record
                skipped = 0
            elseif last_matched !== nothing
                skipped += 1
            end
        end
    end
    matched_nodes = Int[]
    for value in values(stops)
        push!(matched_nodes, value[1])
    end
    service_times = Dict(pair => _integer_median(values) for (pair, values) in observed_times)
    return unique(matched_nodes), services, service_times, observed_times, departures
end

function _transit_wait_seconds(departures, node, service_window_hours, wait_cap_seconds)
    times = get(departures, node, Int[])
    isempty(times) && return 0
    frequency_per_hour = length(times) / max(1.0, service_window_hours)
    return min(wait_cap_seconds, max(1, round(Int, 1800 / frequency_per_hour)))
end

function _reachable_pairs(graph, count, rng)
    candidates = collect(vertices(graph)); pairs = Tuple{Int,Int}[]
    for _ in 1:(20count)
        length(pairs) >= count && break
        o, d = rand(rng, candidates), rand(rng, candidates)
        o != d && has_path(graph, o, d) && push!(pairs, (o,d))
    end
    isempty(pairs) && error("No reachable OSM OD pairs found.")
    return [pairs[mod1(i, length(pairs))] for i in 1:count]
end

"""Generate an Aachen bike/car/transit HNDP instance from cached files."""
function generate_aachen_osm_gtfs(; osm_path="data/osm/aachen_highways.json",
    gtfs_dir="data/osm/avv_gtfs", nusers=100, seed=20260826,
    car_speed_kmh=35, bike_speed_kmh=16, transit_speed_kmh=25,
    bike_station_count=20, bike_station_budget=10,
    all_nodes_bike_stations=false, generalized_cost=true,
    use_boundary_compression=false,
    bike_entry_fee_cents=50, bike_travel_fee_cents_per_minute=5,
    bike_entry_cost_cents=61, bike_cost_per_km_cents=0,
    bike_daily_subscription_cost_cents=nothing,
    car_entry_time_seconds=30, car_parking_search_time_seconds=90,
    car_signal_delay_seconds=15,
    car_value_of_time_eur_per_hour=15, transit_value_of_time_eur_per_hour=10,
    bike_value_of_time_eur_per_hour=12, walk_value_of_time_eur_per_hour=9,
    transit_monthly_cost_cents=6300, transit_rides_per_month=60,
    parking_inner_cost_cents=300, parking_outer_cost_cents=100,
    parking_inner_radius_km=2.0,
    car_local_parking_cost_cents_per_day=71,
    bike_access_time_seconds=60, bike_return_time_seconds=30,
    transit_service_window_hours=16, transit_wait_cap_seconds=600,
    boundary_path="data/osm/aachen_boundary.json")
    transit_rides_per_month > 0 || error("transit_rides_per_month must be positive")
    transit_trip_cost_cents = round(Int, transit_monthly_cost_cents / transit_rides_per_month)
    coords, car_edges, bike_edges, walk_edges = _read_osm_ways(
        osm_path; car_default_speed_kmh=car_speed_kmh, bike_default_speed_kmh=bike_speed_kmh)
    external_shortcut_count = 0
    boundary_compression = use_boundary_compression ||
        basename(osm_path) in ("aachen_city_highways.json", "aachen_buffer_highways.json")
    if boundary_compression
        coords, car_edges, bike_edges, walk_edges = _compress_external_car_edges(coords, car_edges, bike_edges, walk_edges, boundary_path)
        external_shortcut_count = count(a -> a.kind == "external_shortcut", values(car_edges))
    end
    coords, car_edges, bike_edges, walk_edges = _restrict_to_car_component(coords, car_edges, bike_edges, walk_edges)
    bike_time = (e, a) -> max(1, round(Int, 3600*a.distance/a.speed))
    n = length(coords)
    car = DiGraph(n); bike = DiGraph(n)
    for (e, _) in car_edges; add_edge!(car, e...); end
    for (e, _) in bike_edges; add_edge!(bike, e...); end
    # Ways are contracted to arcs, but OSM signal nodes at way endpoints are
    # retained in the arc attributes. Apply delay only to those tagged arcs.
    car_time = (e, a) -> max(1, round(Int, 3600*a.distance/a.speed)) +
        (a.kind == "external_shortcut" || !a.signal ? 0 : car_signal_delay_seconds)
    stop_nodes, transit_services, transit_times, _, stop_departures = _read_gtfs_services(gtfs_dir, coords, 0.25)
    rng = MersenneTwister(seed)
    degrees = [indegree(car,v)+outdegree(car,v) for v in vertices(car)]
    station_candidates = sort(stop_nodes, by=v -> degrees[v], rev=true)
    stations = all_nodes_bike_stations ? collect(1:n) :
        station_candidates[1:min(bike_station_count, length(station_candidates))]
    # GTFS services are independent transit-layer links; they need not be
    # direct road arcs (the whole point is that transit follows its own route).
    transit_services = [(i,j) for (i,j) in transit_services if i != j && haskey(transit_times, (i,j))]
    car_layer = n
    transit_layer = 2n
    bike_layer = 3n
    graph = DiGraph(4n)
    for (e,a) in car_edges; add_edge!(graph,e[1]+car_layer,e[2]+car_layer); end
    for (e,a) in bike_edges; add_edge!(graph,e[1]+bike_layer,e[2]+bike_layer); end
    # Walking is represented directly on the physical base layer.
    for (e,a) in walk_edges; add_edge!(graph,e...); end
    for (i,j) in transit_services; add_edge!(graph,i+transit_layer,j+transit_layer); end
    decision_arcs = Tuple{Int,Int}[]; groups = Vector{Vector{Tuple{Int,Int}}}()
    for station in 1:n
        add_edge!(graph,station,station+car_layer); add_edge!(graph,station+car_layer,station)
    end
    for station in stations
        up=(station,station+bike_layer); down=(station+bike_layer,station)
        add_edge!(graph,up...); add_edge!(graph,down...)
        push!(decision_arcs,up); push!(groups,[up,down]); add_edge!(graph,station,station+transit_layer); add_edge!(graph,station+transit_layer,station)
    end
    for station in unique(vcat(stations, first.(transit_services), last.(transit_services)))
        add_edge!(graph,station,station+transit_layer); add_edge!(graph,station+transit_layer,station)
    end
    m=4n; cost=spzeros(Int,m,m); risk=spzeros(Int,m,m)
    for (e,a) in car_edges
        cost[e[1]+car_layer,e[2]+car_layer] = _mode_cost(
            car_time(e, a), car_value_of_time_eur_per_hour, generalized_cost)
    end
    for (e,a) in bike_edges
        cost[e[1]+bike_layer,e[2]+bike_layer] = _mode_cost(
            bike_time(e, a), bike_value_of_time_eur_per_hour, generalized_cost) +
            (generalized_cost ? (bike_daily_subscription_cost_cents === nothing ?
                round(Int, bike_cost_per_km_cents * a.distance) :
                bike_travel_fee_cents_per_minute * max(1, cld(bike_time(e, a), 60))) : 0)
    end
    for (e,a) in walk_edges
        cost[e...] = _mode_cost(3600*a.distance/5, walk_value_of_time_eur_per_hour, generalized_cost)
    end
    for (i,j) in transit_services
        observed = transit_times[(i,j)]
        cost[i+transit_layer,j+transit_layer] = _mode_cost(
            observed, transit_value_of_time_eur_per_hour, generalized_cost)
    end
    for station in unique(vcat(stations, first.(transit_services), last.(transit_services)))
        wait = _transit_wait_seconds(stop_departures, station, transit_service_window_hours, transit_wait_cap_seconds)
        cost[station,station+transit_layer]=_mode_cost(
            wait, transit_value_of_time_eur_per_hour, generalized_cost) +
            (generalized_cost ? transit_trip_cost_cents : 0)
        cost[station+transit_layer,station]=0
    end
    for station in 1:n
        cost[station,station+car_layer]=_mode_cost(
            car_entry_time_seconds, car_value_of_time_eur_per_hour, generalized_cost)
        parking = car_local_parking_cost_cents_per_day === nothing ?
            (_haversine_km(AACHEN_CENTER..., coords[station]...) <= parking_inner_radius_km ?
                parking_inner_cost_cents : parking_outer_cost_cents) :
            car_local_parking_cost_cents_per_day
        cost[station+car_layer,station]=_mode_cost(
            car_parking_search_time_seconds, car_value_of_time_eur_per_hour, generalized_cost) +
            (generalized_cost ? parking : 0)
    end
    for group in groups
        risk[group[1]...] = -(bike_daily_subscription_cost_cents === nothing ?
            bike_entry_cost_cents : bike_daily_subscription_cost_cents)
        cost[group[1]...] = _mode_cost(
            bike_access_time_seconds, bike_value_of_time_eur_per_hour, generalized_cost) +
            (generalized_cost ? (bike_daily_subscription_cost_cents === nothing ?
                bike_entry_cost_cents : bike_daily_subscription_cost_cents) : 0)
        cost[group[2]...] = _mode_cost(
            bike_return_time_seconds, bike_value_of_time_eur_per_hour, generalized_cost)
    end
    for (e, a) in bike_edges
        bike_arc = (e[1] + bike_layer, e[2] + bike_layer)
        minutes = max(1, cld(bike_time(e, a), 60))
        risk[bike_arc...] = bike_daily_subscription_cost_cents === nothing ?
            -round(Int, bike_cost_per_km_cents * a.distance) : 0
    end
    transit_nodes = unique(vcat(first.(transit_services), last.(transit_services)))
    transit_waits = [_transit_wait_seconds(stop_departures, node, transit_service_window_hours, transit_wait_cap_seconds) for node in transit_nodes]
    users=User[]
    for (k,(o,d)) in enumerate(_reachable_pairs(car,nusers,rng)); push!(users,User("AachenU$(k)",o,d,risk,cost,nothing,nothing)); end
    prices=Dict(a=>1 for a in decision_arcs)
    observed_service_count = length(transit_services)
    city_extract = basename(osm_path) == "aachen_city_highways.json"
    buffered_extract = basename(osm_path) == "aachen_buffer_highways.json"
    instance_name = city_extract ? "aachen_city_boundary_osm_gtfs_bike" : buffered_extract ? "aachen_buffered_boundary_osm_gtfs_bike" : "aachen_osm_gtfs_bike"
    extract_description = city_extract ? "OSM Aachen administrative boundary area" : buffered_extract ? "buffered OSM extract with Aachen-boundary external-road compression" : "rectangular OSM bounding box"
    source_bbox = buffered_extract ? AACHEN_BUFFER_BBOX : AACHEN_BBOX
    meta=Dict{String,Any}("name"=>instance_name, "osm_path"=>osm_path, "gtfs_dir"=>gtfs_dir,
        "osm_extract"=>extract_description,
        "layers"=>["physical/walk", "car", "transit", "bike"], "source_bbox"=>collect(source_bbox), "nnodes"=>4n, "coords"=>coords, "car_arcs"=>length(car_edges), "bike_arcs"=>length(bike_edges), "walk_arcs"=>length(walk_edges), "external_car_shortcuts"=>external_shortcut_count, "transit_services"=>length(transit_services), "transit_services_with_observed_times"=>observed_service_count, "transit_time_source"=>"median GTFS stop-time duration; boundary-spanning services included; untimed services omitted", "time_unit"=>"integer euro-cents generalized cost", "walk_speed_kmh"=>5, "car_entry_time_seconds"=>car_entry_time_seconds, "car_parking_search_time_seconds"=>car_parking_search_time_seconds, "car_signal_delay_seconds"=>car_signal_delay_seconds, "car_value_of_time_eur_per_hour"=>car_value_of_time_eur_per_hour, "transit_value_of_time_eur_per_hour"=>transit_value_of_time_eur_per_hour, "bike_value_of_time_eur_per_hour"=>bike_value_of_time_eur_per_hour, "walk_value_of_time_eur_per_hour"=>walk_value_of_time_eur_per_hour, "parking_inner_cost_cents"=>parking_inner_cost_cents, "parking_outer_cost_cents"=>parking_outer_cost_cents, "parking_inner_radius_km"=>parking_inner_radius_km, "bike_access_time_seconds"=>bike_access_time_seconds, "bike_return_time_seconds"=>bike_return_time_seconds, "transit_service_window_hours"=>transit_service_window_hours, "transit_wait_cap_seconds"=>transit_wait_cap_seconds, "transit_wait_seconds_min"=>(isempty(transit_waits) ? 0 : minimum(transit_waits)), "transit_wait_seconds_median"=>(isempty(transit_waits) ? 0 : _integer_median(transit_waits)), "transit_wait_seconds_mean"=>(isempty(transit_waits) ? 0.0 : sum(transit_waits)/length(transit_waits)), "transit_wait_seconds_max"=>(isempty(transit_waits) ? 0 : maximum(transit_waits)), "bike_entry_fee_cents"=>bike_entry_fee_cents, "bike_travel_fee_cents_per_minute"=>bike_travel_fee_cents_per_minute, "stations"=>stations, "nusers"=>nusers, "seed"=>seed)
    meta["car_local_parking_cost_cents_per_day"] = car_local_parking_cost_cents_per_day
    meta["bike_daily_subscription_cost_cents"] = bike_daily_subscription_cost_cents
    meta["bike_entry_cost_cents"] = bike_entry_cost_cents
    meta["bike_cost_per_km_cents"] = bike_cost_per_km_cents
    meta["transit_monthly_cost_cents"] = transit_monthly_cost_cents
    meta["transit_rides_per_month"] = transit_rides_per_month
    meta["transit_trip_cost_cents"] = transit_trip_cost_cents
    return HNDPGeneratedNetwork(instance_name,HNDPwC(graph,users,decision_arcs,prices,nothing,groups,bike_station_budget),meta)
end

"""Generate the same instance using OSM's Aachen municipal boundary extract."""
function generate_aachen_city_osm_gtfs(; kwargs...)
    return generate_aachen_osm_gtfs(; osm_path="data/osm/aachen_city_highways.json", kwargs...)
end

"""Generate an Aachen instance from a buffered OSM extract with boundary compression."""
function generate_aachen_buffered_osm_gtfs(; kwargs...)
    return generate_aachen_osm_gtfs(; osm_path="data/osm/aachen_buffer_highways.json", kwargs...)
end
