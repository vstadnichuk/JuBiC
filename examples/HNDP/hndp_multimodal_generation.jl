"""Realistic, parameterized multimodal HNDP instance generation.

The generator keeps the source network topology and derives car, bicycle, and
transit times from the distance/capacity fields in the TNTP/SF files.  The
defaults are deliberately exposed as keyword/configuration parameters: they
are useful city-scale starting points, but should be calibrated when a
particular study has better local data.
"""

struct HNDPArcAttributes
    graph::DiGraph
    distance::Dict{Tuple{Int,Int},Float64}
    capacity::Dict{Tuple{Int,Int},Float64}
end

function _load_hndp_arc_attributes(topology_id::String)
    if topology_id == "sioux_falls"
        path = "examples/data/SF_DNDP_10_base.txt"
        rows = _parse_hndp_network_rows(path, true)
        n = 24
    else
        path = topology_id == "anaheim" ? "examples/data/Anaheim_net.tntp" :
            topology_id == "friedrichshain_center" ? "examples/data/friedrichshain-center_net.tntp" :
            topology_id == "berlin_mitte_center" ? "examples/data/berlin-mitte-center_net.tntp" :
            topology_id == "ema" ? "examples/data/EMA_net.tntp" : nothing
        path === nothing && throw(ArgumentError("Unsupported HNDP topology '$topology_id'."))
        rows = _parse_hndp_network_rows(path, false)
        n = parse_tntp_edges(path)[2]
    end
    graph = DiGraph(n)
    distance = Dict{Tuple{Int,Int},Float64}()
    capacity = Dict{Tuple{Int,Int},Float64}()
    for row in rows
        i, j, cap, dist, fft = row
        add_edge!(graph, i, j)
        # Some urban TNTP files contain zero-length connector arcs.  Preserve
        # them as small positive links so shortest-path costs remain useful.
        effective_dist = dist > 0 ? dist : (fft > 0 ? fft : 1.0)
        distance[(i, j)] = effective_dist
        capacity[(i, j)] = max(cap, 0.0)
    end
    return HNDPArcAttributes(graph, distance, capacity)
end

function _parse_hndp_network_rows(path::String, sioux::Bool)
    lines = readlines(path)
    marker = sioux ? findfirst(x -> occursin("Init node", x), lines) :
        findfirst(x -> occursin("<END OF METADATA>", x), lines)
    marker === nothing && throw(ArgumentError("Could not locate network data in $path."))
    rows = NTuple{5,Float64}[]
    for line in lines[(marker + 1):end]
        tokens = split(strip(replace(line, ";" => "")), r"\s+")
        length(tokens) < 5 && continue
        try
            i, j = parse(Int, tokens[1]), parse(Int, tokens[2])
            cap, dist, fft = parse.(Float64, tokens[3:5])
            push!(rows, (Float64(i), Float64(j), cap, dist, fft))
        catch
            continue
        end
    end
    return [(Int(r[1]), Int(r[2]), r[3], r[4], r[5]) for r in rows]
end

function _multimodal_config_value(spec, key, default)
    return haskey(spec, key) ? spec[key] : default
end

function _integer_config_value(spec, key::String, default::Int)
    value = _multimodal_config_value(spec, key, default)
    if value isa Integer
        return Int(value)
    elseif value isa AbstractFloat && isfinite(value) && isinteger(value)
        return Int(round(value))
    end
    throw(ArgumentError("Multimodal parameter '$key' must be integer-valued."))
end

function _basis_point_config(spec, bp_key::String, legacy_key::String, default::Int)
    if haskey(spec, bp_key)
        value = _integer_config_value(spec, bp_key, default)
        0 <= value <= 10_000 || throw(ArgumentError("$bp_key must lie in [0, 10000]."))
        return value
    elseif haskey(spec, legacy_key)
        value = Float64(spec[legacy_key])
        0.0 <= value <= 1.0 || throw(ArgumentError("$legacy_key must lie in [0, 1]."))
        return round(Int, 10_000 * value)
    end
    return default
end

function _build_multimodal_generated_network(
    base_name::String,
    topology_id::String,
    instance_type::String,
    spec::Dict{String,Any},
)
    instance_type in ("multimodal_bike", "multimodal_expansion") ||
        throw(ArgumentError("Unsupported multimodal instance type '$instance_type'."))
    attrs = _load_hndp_arc_attributes(topology_id)
    n = nv(attrs.graph)
    seed = Int(spec["parameter_seed"])
    nusers = Int(spec["nusers"])
    station_rng = MersenneTwister(seed + 1)
    transit_rng = MersenneTwister(seed + 2)
    user_rng = MersenneTwister(seed + 3)

    car_speed = _integer_config_value(spec, "car_speed_kmh", 40)
    bike_speed = _integer_config_value(spec, "bike_speed_kmh", 15)
    degree_penalty_bp = _basis_point_config(spec, "car_degree_penalty_bp", "car_degree_penalty", 1_500)
    capacity_boost_bp = _basis_point_config(spec, "car_capacity_boost_bp", "car_capacity_boost", 1_000)
    bike_capacity_boost_bp = _basis_point_config(spec, "bike_low_capacity_boost_bp", "bike_low_capacity_boost", 1_000)
    delta_min_bp = _basis_point_config(spec, "transit_delta_min_bp", "transit_delta_min", 2_000)
    delta_max_bp = _basis_point_config(spec, "transit_delta_max_bp", "transit_delta_max", 5_000)
    time_scale = _integer_config_value(spec, "time_scale", 1_000)
    station_limit = if haskey(spec, "station_fraction_percent")
        fraction = _integer_config_value(spec, "station_fraction_percent", 10)
        0 <= fraction <= 100 || throw(ArgumentError("station_fraction_percent must lie in [0, 100]."))
        floor(Int, fraction * n / 100)
    else
        _integer_config_value(spec, "station_limit", 8)
    end
    station_limit = min(max(station_limit, 0), n)
    station_limit >= 0 || throw(ArgumentError("station_limit must be nonnegative."))
    delta_min_bp <= delta_max_bp || throw(ArgumentError("transit delta minimum must not exceed maximum."))
    car_speed > 0 && bike_speed > 0 && time_scale > 0 || throw(ArgumentError("Speeds and time_scale must be positive."))

    degrees = [indegree(attrs.graph, v) + outdegree(attrs.graph, v) for v in 1:n]
    stations = _select_spread_stations(attrs.graph, degrees, station_limit, station_rng)
    additional_station_count = if instance_type == "multimodal_expansion"
        if haskey(spec, "expansion_station_fraction_percent")
            fraction = _integer_config_value(spec, "expansion_station_fraction_percent", 10)
            0 <= fraction <= 100 || throw(ArgumentError("expansion_station_fraction_percent must lie in [0, 100]."))
            floor(Int, fraction * (n - length(stations)) / 100)
        else
            _integer_config_value(spec, "expansion_station_count", 4)
        end
    else
        0
    end
    additional_station_count = min(additional_station_count, n - length(stations))
    additional_station_count >= 0 || throw(ArgumentError("expansion_station_count must be nonnegative."))
    additional_stations = instance_type == "multimodal_expansion" ?
        _select_spread_stations(attrs.graph, degrees, additional_station_count, station_rng; initial=stations) : Int[]
    transit_stations = vcat(stations, additional_stations)
    cap_values = collect(values(attrs.capacity))
    cap_min, cap_max = isempty(cap_values) ? (0.0, 1.0) : (minimum(cap_values), maximum(cap_values))
    max_degree = max(maximum(degrees; init=0), 1)

    car_time = Dict{Tuple{Int,Int},Int}()
    bike_time = Dict{Tuple{Int,Int},Int}()
    for edge in keys(attrs.distance)
        d = attrs.distance[edge]
        cap_norm = _normalise_capacity(attrs.capacity[edge], cap_min, cap_max)
        degree_norm = (degrees[edge[1]] + degrees[edge[2]]) / (2max_degree)
        car_value = time_scale * d / car_speed *
            (1 + (degree_penalty_bp / 10_000) * degree_norm) *
            (1 - (capacity_boost_bp / 10_000) * cap_norm)
        bike_value = time_scale * d / bike_speed *
            (1 - (bike_capacity_boost_bp / 10_000) * (1 - cap_norm))
        car_time[edge] = max(1, round(Int, car_value))
        bike_time[edge] = max(1, round(Int, bike_value))
    end

    # The fixed public-transport network is a complete directed subgraph on
    # the preselected base stations.  Expansion adds every directed service
    # involving at least one additional station.
    fixed_transit_pairs = Tuple{Int,Int}[]
    for i in stations, j in stations
        i == j && continue
        push!(fixed_transit_pairs, (i, j))
    end
    expansion_pairs = Tuple{Int,Int}[]
    for i in transit_stations, j in transit_stations
        i == j && continue
        (i in additional_stations || j in additional_stations) || continue
        push!(expansion_pairs, (i, j))
    end
    transit_pairs = vcat(fixed_transit_pairs, expansion_pairs)
    transit_time = Dict{Tuple{Int,Int},Int}()
    for pair in transit_pairs
        # Public transport follows a high-quality version of the car-mode
        # connection between the two stations, not merely a direct source arc.
        baseline = _shortest_car_time(attrs.graph, car_time, pair[1], pair[2])
        delta_bp = rand(transit_rng, delta_min_bp:delta_max_bp)
        transit_time[pair] = max(1, round(Int, baseline * (10_000 - delta_bp) / 10_000))
    end

    layer_count = instance_type == "multimodal_bike" ? 3 : 2
    graph = DiGraph(layer_count * n)
    car_layer = 0
    transit_layer = n
    bike_layer = 2n
    for edge in keys(attrs.distance)
        add_edge!(graph, edge...)
        if instance_type == "multimodal_bike"
            add_edge!(graph, edge[1] + bike_layer, edge[2] + bike_layer)
        end
    end
    for pair in fixed_transit_pairs
        add_edge!(graph, pair[1] + transit_layer, pair[2] + transit_layer)
    end

    decision_arcs = Tuple{Int,Int}[]
    decision_groups = Vector{Vector{Tuple{Int,Int}}}()
    if instance_type == "multimodal_bike"
        for station in stations
            up = (station, station + bike_layer)
            down = (station + bike_layer, station)
            add_edge!(graph, up...); add_edge!(graph, down...)
            append!(decision_arcs, (up, down))
            push!(decision_groups, [up, down])
        end
        # Car-to-transit transfers are fixed and available only at stations.
        for station in stations
            add_edge!(graph, station, station + transit_layer)
            add_edge!(graph, station + transit_layer, station)
        end
    else
        # Public transport is part of the fixed base. Expansion candidates are
        # additional public-transport-layer links, not car-layer links.
        for station in transit_stations
            add_edge!(graph, station, station + transit_layer)
            add_edge!(graph, station + transit_layer, station)
        end
        for pair in expansion_pairs
            transit_arc = (pair[1] + transit_layer, pair[2] + transit_layer)
            add_edge!(graph, transit_arc...)
            push!(decision_arcs, transit_arc)
            push!(decision_groups, [transit_arc])
        end
    end

    m = layer_count * n
    rcost = zeros(Int, m, m)
    rrisk = zeros(Int, m, m)
    for edge in keys(attrs.distance)
        rcost[edge...] = car_time[edge]
        rcost[edge[1] + transit_layer, edge[2] + transit_layer] =
            car_time[edge]
        if instance_type == "multimodal_bike"
            rcost[edge[1] + bike_layer, edge[2] + bike_layer] = bike_time[edge]
        end
    end
    for pair in transit_pairs
        transit_arc = (pair[1] + transit_layer, pair[2] + transit_layer)
        rcost[transit_arc...] = transit_time[pair]
        if pair in expansion_pairs
            rrisk[transit_arc...] = -_integer_config_value(spec, "expansion_arc_profit", 1)
        end
    end
    if instance_type == "multimodal_bike"
        bike_profit = _integer_config_value(spec, "bike_station_profit", 1)
        for group in decision_groups
            rrisk[group[1]...] = -bike_profit
            rcost[group[1]...] = 0
            rcost[group[2]...] = 0
        end
    end
    for station in transit_stations
        for edge in ((station, station + transit_layer), (station + transit_layer, station))
            rcost[edge...] = 0
        end
    end

    feasible_pairs = Tuple{Int,Int}[]
    for origin in 1:n, destination in 1:n
        origin == destination && continue
        # Sample from OD pairs reachable on the fixed car topology.  The path
        # reformulation computes its bound after decision arcs are removed, so
        # reachability through a candidate expansion/bike transfer is not
        # sufficient here.
        _hndp_has_directed_path(attrs.graph, origin, destination) && push!(feasible_pairs, (origin, destination))
    end
    isempty(feasible_pairs) && throw(ArgumentError("Multimodal topology '$topology_id' has no reachable car-layer OD pairs."))

    users = User[]
    for user_index in 1:nusers
        origin, destination = rand(user_rng, feasible_pairs)
        push!(users, User("U$(user_index)", origin, destination, rrisk, rcost, nothing, nothing))
    end
    edge_price = Dict(arc => _integer_config_value(
        spec,
        instance_type == "multimodal_bike" ? "bike_station_cost" : "expansion_arc_cost",
        1,
    ) for arc in decision_arcs)

    budget = if instance_type == "multimodal_bike"
        if haskey(spec, "station_budget_fraction_percent")
            fraction = _integer_config_value(spec, "station_budget_fraction_percent", 100)
            0 <= fraction <= 100 || throw(ArgumentError("station_budget_fraction_percent must lie in [0, 100]."))
            floor(Int, fraction * station_limit / 100)
        else
            _integer_config_value(spec, "station_budget", station_limit)
        end
    else
        if haskey(spec, "arc_budget_fraction_percent")
            fraction = _integer_config_value(spec, "arc_budget_fraction_percent", 100)
            0 <= fraction <= 100 || throw(ArgumentError("arc_budget_fraction_percent must lie in [0, 100]."))
            floor(Int, fraction * length(decision_arcs) / 100)
        else
            _integer_config_value(spec, "arc_budget", length(decision_arcs))
        end
    end
    budget = min(max(budget, 0), instance_type == "multimodal_bike" ? station_limit : length(decision_arcs))
    budget >= 0 || throw(ArgumentError("The decision budget must be nonnegative."))
    metadata = Dict{String,Any}(
        "name" => "$(base_name)_$(topology_id)_$(instance_type)_u$(nusers)_s$(seed)",
        "instance_type" => instance_type,
        "topology_family" => topology_id,
        "scenario" => instance_type == "multimodal_bike" ? "bike" : "expansion",
        "nusers" => nusers,
        "parameter_seed" => seed,
        "stations" => stations,
        "additional_stations" => additional_stations,
        "station_limit" => station_limit,
        "station_fraction_percent" => haskey(spec, "station_fraction_percent") ? _integer_config_value(spec, "station_fraction_percent", 0) : nothing,
        "expansion_station_count" => additional_station_count,
        "expansion_station_fraction_percent" => haskey(spec, "expansion_station_fraction_percent") ? _integer_config_value(spec, "expansion_station_fraction_percent", 0) : nothing,
        "decision_budget" => budget,
        "availability_budget_count" => instance_type == "multimodal_bike" ? 2budget : budget,
        "decision_arcs" => length(decision_arcs),
        "nnodes" => nv(graph),
        "narcs" => ne(graph),
        "car_speed_kmh" => car_speed,
        "bike_speed_kmh" => bike_speed,
        "time_scale" => time_scale,
        "transit_delta_min_bp" => delta_min_bp,
        "transit_delta_max_bp" => delta_max_bp,
        "car_degree_penalty_bp" => degree_penalty_bp,
        "car_capacity_boost_bp" => capacity_boost_bp,
        "bike_low_capacity_boost_bp" => bike_capacity_boost_bp,
        "fixed_transit_pairs" => length(fixed_transit_pairs),
        "expansion_pairs" => length(expansion_pairs),
    )
    return HNDPGeneratedNetwork(metadata["name"], HNDPwC(graph, users, decision_arcs, edge_price, nothing, decision_groups, budget), metadata)
end

function _select_spread_stations(
    graph::DiGraph,
    degrees::Vector{Int},
    count::Int,
    rng::AbstractRNG;
    initial::Vector{Int}=Int[],
)
    count <= 0 && return Int[]
    initial_set = Set(initial)
    selected = copy(initial)
    candidates = [node for node in vertices(graph) if !(node in initial_set)]
    shuffle!(rng, candidates)
    sort!(candidates; by=node -> -degrees[node])

    function is_separated(node)
        neighbours = Set{Int}(vcat(collect(outneighbors(graph, node)), collect(inneighbors(graph, node))))
        return all(!(selected_node in neighbours) for selected_node in selected)
    end

    # First build a high-degree independent-set-like selection.
    for node in candidates
        length(selected) - length(initial) >= count && break
        is_separated(node) || continue
        push!(selected, node)
    end
    # If the topology cannot support the requested separation, fill the
    # remaining slots by degree. This preserves the requested station count.
    for node in candidates
        length(selected) - length(initial) >= count && break
        node in selected && continue
        push!(selected, node)
    end
    return selected[(length(initial) + 1):end]
end

function _shortest_car_time(graph::DiGraph, weights::Dict{Tuple{Int,Int},Int}, source::Int, target::Int)
    n = nv(graph)
    for undirected in (false, true)
        distances = fill(typemax(Int), n)
        used = falses(n)
        distances[source] = 0
        for _ in 1:n
            node = 0
            best = typemax(Int)
            for candidate in 1:n
                !used[candidate] && distances[candidate] < best || continue
                node = candidate
                best = distances[candidate]
            end
            node == 0 && break
            node == target && return distances[node]
            used[node] = true
            neighbours = undirected ?
                union(Set(outneighbors(graph, node)), Set(inneighbors(graph, node))) :
                Set(outneighbors(graph, node))
            for next_node in neighbours
                edge = (node, next_node)
                weight = if haskey(weights, edge)
                    weights[edge]
                elseif undirected && haskey(weights, (next_node, node))
                    weights[(next_node, node)]
                else
                    continue
                end
                candidate_distance = distances[node] + weight
                candidate_distance < distances[next_node] && (distances[next_node] = candidate_distance)
            end
        end
    end
    return 1
end

function _hndp_has_directed_path(graph::DiGraph, source::Int, target::Int)
    source == target && return true
    visited = Set{Int}([source])
    queue = Int[source]
    while !isempty(queue)
        node = popfirst!(queue)
        for next_node in outneighbors(graph, node)
            next_node == target && return true
            if !(next_node in visited)
                push!(visited, next_node)
                push!(queue, next_node)
            end
        end
    end
    return false
end

function _normalise_capacity(value::Float64, lower::Float64, upper::Float64)
    upper <= lower && return 0.5
    return clamp((value - lower) / (upper - lower), 0.0, 1.0)
end

function _pair_distance(attrs::HNDPArcAttributes, pair::Tuple{Int,Int})
    if haskey(attrs.distance, pair)
        return attrs.distance[pair]
    end
    # A direct transit connection is based on the shortest car-network metric;
    # fall back to one kilometre only for disconnected/degenerate inputs.
    weights = zeros(Float64, nv(attrs.graph), nv(attrs.graph))
    for edge in keys(attrs.distance)
        weights[edge...] = attrs.distance[edge]
    end
    paths = floyd_warshall_shortest_paths(attrs.graph, weights)
    value = paths.dists[pair...]
    return isfinite(value) && value > 0 ? value : 1.0
end
