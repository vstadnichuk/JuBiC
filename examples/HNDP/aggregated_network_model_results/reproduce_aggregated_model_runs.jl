using CSV
using DataFrames
using Graphs
using JuBiC

include(joinpath(@__DIR__, "reduced_station_balanced_spatial_sweep.jl"))
include(joinpath(@__DIR__, "run_fixed_node_bike_expansion_benchmark.jl"))

const HARD_CITIES = [
    (name="cologne", osm="data/osm/comparative_cities/cologne_buffered.json", boundary="data/osm/comparative_cities/cologne_boundary.json", gtfs="data/osm/comparative_cities/cologne_gtfs"),
    (name="heidelberg", osm="data/osm/comparative_cities/heidelberg.json", boundary="data/osm/comparative_cities/heidelberg_boundary.json", gtfs="data/osm/comparative_cities/heidelberg_gtfs"),
    (name="karlsruhe", osm="data/osm/comparative_cities/karlsruhe_buffered.json", boundary="data/osm/comparative_cities/karlsruhe_boundary.json", gtfs="data/osm/comparative_cities/karlsruhe_gtfs"),
]
const HARD_TARGETS = [50, 100, 150, 200, 250, 300]
const HARD_USERS = [100, 300, 500, 750, 1000]
const HARD_BUDGET = 0.50
const HARD_LIMIT = 600.0
const HARD_SEED = 20260829
# Reproduce the published/archived grid from scratch. Results are written to
# tmp_compare so reruns remain outside the tracked result package.
const HARD_OLD = String[]
const HARD_OUT = "tmp_compare/network_balanced_hardness_sweep_reproduced.csv"

const HARD_MODELS = [
    Dict{String,Any}("name"=>"sd_fixed_network_path", "model_type"=>"sd", "big_m_mode"=>"fixed_network_path", "indicator_constraints"=>false, "bound_duals"=>true),
    Dict{String,Any}("name"=>"sd_n_minus_1", "model_type"=>"sd", "big_m_mode"=>"n_minus_one_most_expensive", "indicator_constraints"=>false, "bound_duals"=>true),
    Dict{String,Any}("name"=>"blc_astar_fixed_network_path", "model_type"=>"blc", "big_m_mode"=>"fixed_network_path", "subproblem_method"=>"astar"),
    Dict{String,Any}("name"=>"blc_astar_n_minus_1", "model_type"=>"blc", "big_m_mode"=>"n_minus_one_most_expensive", "subproblem_method"=>"astar"),
    Dict{String,Any}("name"=>"path_dominance", "model_type"=>"path", "parallelize"=>false, "use_decision_arc_dominance"=>true, "enumeration_time_limit"=>HARD_LIMIT, "fixed_path_bound_method"=>"labeling"),
    Dict{String,Any}("name"=>"path_no_dominance", "model_type"=>"path", "parallelize"=>false, "use_decision_arc_dominance"=>false, "enumeration_time_limit"=>HARD_LIMIT, "fixed_path_bound_method"=>"labeling"),
]
const HARD_NAMES = Set(String(m["name"]) for m in HARD_MODELS)

function hard_write(rows)
    isempty(rows) && return
    cols = sort!(collect(Set(vcat([collect(keys(r)) for r in rows]...))))
    mkpath(dirname(HARD_OUT))
    CSV.write(HARD_OUT, DataFrame([Symbol(c)=>[get(r,c,missing) for r in rows] for c in cols]))
end

function hard_key(r)
    return (String(r.city), Int(r.retained_stations), Int(r.nusers), round(Int, 100Float64(r.budget_fraction)), String(r.model_name))
end

function reuse_hard_rows(rows, seen)
    for source in HARD_OLD
        isfile(source) || continue
        old = CSV.read(source, DataFrame)
        for r in eachrow(old)
            try
                city = String(r.city); target = Int(r.retained_stations); users = Int(r.nusers)
                budget = round(Int, 100Float64(r.budget_fraction)); model = String(r.model_name)
                target in HARD_TARGETS && users in HARD_USERS && budget == 50 && model in HARD_NAMES || continue
                key = (city, target, users, budget, model)
                key in seen && continue
                old_status = hasproperty(r, :Opt_status) ? String(r.Opt_status) : ""
                old_status in ("Optimal", "Timelimit", "Memory_Limit") || continue
                d = Dict{String,Any}(String(k)=>r[k] for k in names(old))
                old_status == "Timelimit" && (d["run_status"] = "timeout")
                old_status == "Memory_Limit" && (d["run_status"] = "memory_limit")
                d["reused_previous_result"] = true
                d["sweep_name"] = "network_balanced_hardness_sweep"
                push!(rows, d); push!(seen, key)
            catch
                continue
            end
        end
    end
end

function hard_metric(generated)
    h = generated.instance; n = generated.metadata["nnodes"] ÷ 4
    stations = generated.metadata["stations"]; costs = h.users[1].mcost
    car_allowed = (u,v) -> (u <= n && n < v <= 2n) || (n < u <= 2n && v <= n) || (n < u <= 2n && n < v <= 2n)
    car = station_distances_to_targets(h.mygraph, costs, stations, stations, car_allowed)
    metric = copy(car); inf = typemax(Int)
    for i in axes(metric,1), j in axes(metric,2)
        metric[i,j] = metric[i,j] == inf ? metric[j,i] : min(metric[i,j], metric[j,i])
    end
    stations, metric
end

function main()
    rows = Dict{String,Any}[]; seen = Set{Tuple{String,Int,Int,Int,String}}()
    reuse_hard_rows(rows, seen)
    println("Reused ", length(rows), " completed rows")
    hard_write(rows)
    all_threads = Sys.CPU_THREADS
    println("Julia threads available: ", Threads.nthreads(), "; solver threads configured: ", all_threads)
    params = Dict{String,Any}("name"=>"gurobi_10min", "mip_solver"=>"Gurobi", "runtime"=>HARD_LIMIT,
        "threads_master"=>all_threads, "threads_sub_con"=>all_threads,
        "parallel_separation"=>true, "seed"=>HARD_SEED, "write_run_logs"=>false)
    for city in HARD_CITIES
        println("Generating full network ", city.name)
        generated = generate_aachen_osm_gtfs(osm_path=city.osm, boundary_path=city.boundary,
            use_boundary_compression=true, gtfs_dir=city.gtfs, nusers=1, seed=HARD_SEED,
            bike_station_count=10000, all_nodes_bike_stations=false, generalized_cost=true,
            bike_entry_cost_cents=50, bike_cost_per_km_cents=5,
            bike_daily_subscription_cost_cents=nothing)
        stations, metric = hard_metric(generated)
        for target in HARD_TARGETS
            selected = stations[balanced_medoids(metric, target)[1]]
            for nusers in HARD_USERS
                nusers == 100 && target in (50,100,150,200) && continue
                # These combinations already caused memory exhaustion at full
                # solver parallelism. Exclude them so memory failures do not
                # dominate the sweep; retain the adjacent hard cases.
                target >= 200 && nusers >= 750 && continue
                agg = build_aggregated(generated, target, length(selected), HARD_BUDGET,
                    nusers, HARD_SEED + target + nusers; selected_override=selected)
                for model in HARD_MODELS
                    key = (city.name, target, nusers, 50, String(model["name"]))
                    key in seen && continue
                    case_id = "$(city.name)_network_balanced_n$(target)_b50_u$(nusers)_$(model["name"])"
                    println(case_id)
                    row = Dict{String,Any}(); started = time()
                    try
                        row = run_case(agg, model, params, case_id)
                    catch err
                        row["Error"] = sprint(showerror, err); row["run_status"] = "error"
                    end
                    row["city"] = city.name; row["retained_stations"] = target
                    row["station_selection"] = "network_balanced"; row["budget_fraction"] = HARD_BUDGET
                    row["nusers"] = nusers; row["reused_previous_result"] = false
                    row["sweep_name"] = "network_balanced_hardness_sweep"
                    row["wall_time_seconds"] = get(row, "wall_time_seconds", time()-started)
                    push!(rows, row); push!(seen, key); hard_write(rows)
                end
            end
        end
    end
    println("Wrote ", length(rows), " rows to ", HARD_OUT)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
