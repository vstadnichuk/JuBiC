#= Demand-weighted trip distances by mode on the full city networks.

Scenario, routing, demand weights and tie rule are those of
full_network_city_profiles.jl and resolve_full_network_ties.py: all existing
sharing candidates open, every ordered nonself OD for both availability
classes, survey mass uniform within distance bands, and ties resolved by
bike priority, then car, then transit. The trip distance of an OD is the
physical length of its class-specific fixed-network reference route, i.e.,
the length that assigns the OD to a distance band. No optimization is solved.

Usage: julia --project=<repo> --threads=8 mode_distance_profiles.jl <output-dir> [car-cache-root] [noncar-cache-root] [city ...]
Both cache roots default to the cache shared with full_network_city_profiles.jl.
Reads the study folder (SLP_STUDY_DIR) and the fixed reference caches; writes only to <output-dir>.
=#
using CSV, JSON, Base.Threads
include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(SRC, "demand", "osm_mode_share_analysis.jl"))
include(joinpath(SRC, "demand", "osm_fixed_network_demand.jl"))
const MS = OSMModeShares
const FD = OSMFixedNetworkDemand
const FAMILIES = ["walking", "car", "transit", "bike_sharing", "multi_powered"]
const STEP_KM = 0.1          # resolution of the distance histogram
const NBINS = 501            # last bin collects all trips longer than 50 km
const REPORT_BINS = [1.0, 2.0, 3.0, 5.0, 10.0]   # inclusive upper bounds; last bin > 10 km

family(mask) = (nw = mask & 14; nw == 0 ? 1 : nw == 2 ? 2 : nw == 4 ? 3 : nw == 8 ? 4 : 5)

# Same rule as resolve_full_network_ties.py: one optimal route, bike priority, then car, then transit.
function select_mask(bits)
    masks = [m for m in 0:15 if bits & (1 << m) != 0]
    isempty(masks) && error("Unreachable OD")
    choice = filter(m -> m & 8 != 0, masks)
    isempty(choice) && (choice = filter(m -> m & 2 != 0, masks))
    isempty(choice) && (choice = filter(m -> m & 4 != 0, masks))
    isempty(choice) && (choice = masks)
    return argmin(m -> (count_ones(m & 14), m & 14, m), choice)
end

function smoke()
    @assert FAMILIES[family(select_mask((1 << 3) | (1 << 9)))] == "bike_sharing"
    @assert FAMILIES[family(select_mask((1 << 3) | (1 << 13)))] == "multi_powered"
    @assert FAMILIES[family(select_mask((1 << 1) | (1 << 3)))] == "car"
    @assert FAMILIES[family(select_mask(1 << 1))] == "walking"
end

distance_bin(km) = clamp(ceil(Int, km / STEP_KM - 1e-9), 1, NBINS)

# hist[family, band, distance bin] = number of ODs
function route_histogram(g, lengths, bounds, lookup)
    n = g.physical_nodes; nb = length(bounds) + 1; workers = min(n, Threads.nthreads())
    hists = [zeros(Int64, 5, nb, NBINS) for _ in 1:workers]
    Threads.@threads :static for worker in 1:workers
        w = MS.Workspace(g); h = hists[worker]
        for origin in worker:workers:n
            MS.shortest_modes!(w, g, origin)
            @inbounds for dest in 1:n
                origin == dest && continue
                km = lengths[dest, origin]
                isfinite(km) || error("Unreachable fixed reference $origin->$dest")
                h[lookup[Int(w.labels[dest]) + 1], searchsortedfirst(bounds, km), distance_bin(km)] += 1
            end
        end
    end
    return reduce(+, hists)
end

function run_city(city, cache_roots)
    root = joinpath(B, "inputs", "city_networks_walk_connected")
    profile = only(filter(x -> x["city"] == city, JSON.parsefile(joinpath(PROVENANCE, "osm_distance_demand_sources.json"))["cities"]))
    owner = only(filter(x -> x["city"] == city, JSON.parsefile(joinpath(PROVENANCE, "osm_car_ownership_sources.json"))["cities"]))
    p = owner["numerator"] / owner["denominator"]
    units = Float64.(profile["mass_units"]); mass = units ./ sum(units)
    bounds = Float64.(profile["upper_bounds_km"])
    g, info = MS.load_city(root, city; station_scope = "gtfs_candidates")
    g = MS.with_car_costs(g, joinpath(root, "speed_overrides", city * ".csv"))
    lookup = [bits == 0 ? 0 : family(select_mask(bits)) for bits in 0:65535]
    weights = Dict{String,Matrix{Float64}}()   # demand mass per family, band-weighted, per distance bin
    for (class, available, weight) in [("car", true, p), ("noncar", false, 1 - p)]
        fixed, _ = FD.load_fixed_graph(joinpath(root, "graphs", city),
            [joinpath(B, "inputs", "geometry", "road_lengths", city * ".csv"), joinpath(B, "inputs", "geometry", "transit_lengths", city * ".csv")];
            car_available = available, car_override = joinpath(root, "speed_overrides", city * ".csv"))
        cache = joinpath(cache_roots[class], city, "paths_" * class)
        marker = JSON.parsefile(joinpath(cache, "complete.json"))
        @assert marker["graph_sha256"] == FD.fingerprint(fixed) && marker["physical_nodes"] == g.physical_nodes
        lengths = FD.matrix_file(joinpath(cache, "lengths_km.f64"), Float64, g.physical_nodes)
        graph = available ? g : MS.without_car(g)
        started = time()
        h = route_histogram(graph, lengths, bounds, lookup)
        println(city, " ", class, " routed in ", round(time() - started, digits = 1), " s"); flush(stdout)
        counts = vec(sum(h, dims = (1, 3)))
        @assert sum(counts) == g.physical_nodes * (g.physical_nodes - 1)
        all(counts[band] > 0 || iszero(mass[band]) for band in eachindex(mass)) || error("Unsupported positive-mass distance band")
        wt = zeros(Float64, 5, NBINS)
        for band in eachindex(mass), f in 1:5, b in 1:NBINS
            h[f, band, b] == 0 && continue
            wt[f, b] += weight * mass[band] * h[f, band, b] / counts[band]
        end
        weights[class] = wt
        lengths = nothing; fixed = nothing; graph = nothing; GC.gc()
    end
    return weights, info
end

centers() = [(b - 0.5) * STEP_KM for b in 1:NBINS]

function quantile_km(wt, q)
    total = sum(wt); acc = 0.0
    for b in eachindex(wt)
        acc += wt[b]
        acc >= q * total - 1e-15 && return b * STEP_KM
    end
    return NBINS * STEP_KM
end

function main()
    smoke()
    output = abspath(ARGS[1]); mkpath(output)
    shared = joinpath(B, "work", "fixed_reference_cache")   # shared with full_network_city_profiles.jl
    cache_roots = Dict("car" => length(ARGS) > 1 ? abspath(ARGS[2]) : shared,
        "noncar" => length(ARGS) > 2 ? abspath(ARGS[3]) : shared)
    cities = length(ARGS) > 3 ? ARGS[4:end] : ["aachen", "cologne", "heidelberg", "karlsruhe", "nuremberg", "dresden"]
    summary = NamedTuple[]; splitrows = NamedTuple[]
    upper = vcat(REPORT_BINS, Inf)
    bin_of(b) = searchsortedfirst(upper, b * STEP_KM - 1e-9)
    for city in cities
        weights, info = run_city(city, cache_roots)
        for class in ["all", "car", "noncar"]
            wt = class == "all" ? weights["car"] .+ weights["noncar"] : weights[class]
            classtotal = sum(wt)
            for f in 1:5
                row = wt[f, :]; total = sum(row)
                total == 0 && continue
                push!(summary, (; city, class, family = FAMILIES[f], demand_share_pct = 100total / classtotal,
                    mean_km = sum(row .* centers()) / total, median_km = quantile_km(row, 0.5),
                    p90_km = quantile_km(row, 0.9),
                    (Symbol("within_$(Int(u))km_pct") => 100sum(row[1:round(Int, u / STEP_KM)]) / total for u in REPORT_BINS)...))
            end
            if class == "all"
                for (k, u) in enumerate(upper)
                    mask = [bin_of(b) == k for b in 1:NBINS]
                    bintotal = sum(wt[:, mask])
                    push!(splitrows, (; city, upper_km = u, demand_share_pct = 100bintotal / classtotal,
                        (Symbol(FAMILIES[f] * "_pct") => 100sum(wt[f, mask]) / bintotal for f in 1:5)...))
                end
                println("CHECK ", city, " shares ", [round(100sum(wt[f, :]) / classtotal, digits = 4) for f in 1:5]); flush(stdout)
            end
        end
        CSV.write(joinpath(output, "mode_distance_summary.csv"), summary)
        CSV.write(joinpath(output, "modal_split_by_distance.csv"), splitrows)
        CSV.write(joinpath(output, "modal_split_by_distance_$(city).csv"), filter(r -> r.city == city, splitrows))
    end
end
main()
