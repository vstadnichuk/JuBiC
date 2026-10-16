using CSV
using DataFrames
using Graphs
using Random
using JuBiC

include(joinpath(@__DIR__, "hndp_experiment_runner.jl"))
include(joinpath(@__DIR__, "aachen_osm_gtfs_generation.jl"))

const CITIES = [
    (name="aachen", osm="data/osm/aachen_buffer_highways.json", boundary="data/osm/aachen_boundary.json", gtfs="data/osm/avv_gtfs"),
    (name="cologne", osm="data/osm/comparative_cities/cologne_buffered.json", boundary="data/osm/comparative_cities/cologne_boundary.json", gtfs="data/osm/comparative_cities/cologne_gtfs"),
    (name="heidelberg", osm="data/osm/comparative_cities/heidelberg.json", boundary="data/osm/comparative_cities/heidelberg_boundary.json", gtfs="data/osm/comparative_cities/heidelberg_gtfs"),
]
const NODE_TARGETS = [100, 200, 300, 400, 500]
const STATION_COUNTS = [10, 20, 30, 40, 50, 75, 100]
const USER_COUNTS = [100, 250, 500, 1000]
const BUDGET_FRACTIONS = [0.10, 0.25, 0.50]
const SEED = 20260827
const LIMIT = 600.0
const BENCHMARK_KIND = get(ENV, "JUBIC_BENCHMARK_KIND", "sdblc")
const SUBSET_USERS = parse(Int, get(ENV, "JUBIC_SUBSET_USERS", "100"))
const OUT = BENCHMARK_KIND == "path" ?
    "tmp_compare/fixed_node_bike_expansion_path_user_budget_benchmark.csv" :
    BENCHMARK_KIND == "path_subset" ?
    "tmp_compare/fixed_node_bike_expansion_path_cologne_heidelberg_u$(SUBSET_USERS)_subset.csv" :
    "tmp_compare/fixed_node_bike_expansion_user_budget_benchmark.csv"

function heap_push!(h, x)
    push!(h, x); i=length(h)
    while i > 1
        p=fld(i,2); h[p][1] <= h[i][1] && break
        h[p],h[i]=h[i],h[p]; i=p
    end
end
function heap_pop!(h)
    x=h[1]; y=pop!(h); isempty(h) && return x; h[1]=y; i=1
    while true
        l,r,s=2i,2i+1,i
        l<=length(h) && h[l][1]<h[s][1] && (s=l)
        r<=length(h) && h[r][1]<h[s][1] && (s=r)
        s==i && break
        h[i],h[s]=h[s],h[i]; i=s
    end
    x
end

function mode_closure(graph, costs, risks, sources, allowed)
    inf=typemax(Int); k=length(sources); n=nv(graph)
    dmat=fill(inf,k,k); rmat=fill(0,k,k); adjacency=[Tuple{Int,Int,Int}[] for _ in 1:n]
    for u in 1:n, v in outneighbors(graph,u)
        allowed(u,v) || continue
        costs[u,v]==inf || push!(adjacency[u],(v,costs[u,v],risks[u,v]))
    end
    for (si,source) in enumerate(sources)
        dist=fill(inf,n); pred=fill(0,n); dist[source]=0; q=[(0,source)]
        while !isempty(q)
            cur,u=heap_pop!(q); cur==dist[u] || continue
            for (v,w,_) in adjacency[u]
                z=cur+w; z<dist[v] || continue
                dist[v]=z; pred[v]=u; heap_push!(q,(z,v))
            end
        end
        pathrisk=fill(0,n)
        for u in sortperm(dist)
            pred[u]==0 && continue
            pathrisk[u]=pathrisk[pred[u]]+risks[pred[u],u]
        end
        for (di,destination) in enumerate(sources)
            dmat[si,di]=dist[destination]
            dmat[si,di]==inf || (rmat[si,di]=pathrisk[destination])
        end
    end
    dmat,rmat
end

function path_redundancy_ranks(graph, costs, sources, allowed)
    inf=typemax(Int); k=length(sources); n=nv(graph); result=fill(inf,k,k)
    adjacency=[Tuple{Int,Int}[] for _ in 1:n]
    for u in 1:n, v in outneighbors(graph,u)
        allowed(u,v) || continue
        costs[u,v]==inf || push!(adjacency[u],(v,costs[u,v]))
    end
    ranks=Dict(v=>i for (i,v) in enumerate(sources))
    for (si,source) in enumerate(sources)
        dist=fill(inf,n); pred=fill(0,n); best=fill(inf,n); dist[source]=0; q=[(0,source)]
        while !isempty(q)
            cur,u=heap_pop!(q); cur==dist[u] || continue
            for (v,w) in adjacency[u]
                z=cur+w; z<dist[v] || continue
                dist[v]=z; pred[v]=u; heap_push!(q,(z,v))
            end
        end
        for u in sortperm(dist)
            pred[u]==0 && continue
            best[u]=best[pred[u]]
            haskey(ranks,u) && (best[u]=min(best[u],ranks[u]))
        end
        for (di,destination) in enumerate(sources)
            destination==source && continue
            pred[destination]==0 || (result[si,di]=best[pred[destination]])
        end
    end
    result
end

function select_nodes(metadata, n, target)
    k=min(target,n); stations=[s for s in metadata["stations"] if s<=n]
    selected=unique(vcat(stations,collect(1:n)))
    selected[1:k]
end

function build_aggregated(generated, target_nodes, possible_stations, budget_fraction, nusers, seed; selected_override=nothing)
    base=generated.instance; meta=generated.metadata; n=meta["nnodes"] ÷ 4
    selected=isnothing(selected_override) ? select_nodes(meta,n,target_nodes) : collect(selected_override); k=length(selected)
    original=base.mygraph; costs=base.users[1].mcost; risks=base.users[1].mrisk
    car_ok=(u,v)->n<u<=2n && n<v<=2n
    walk_ok=(u,v)->u<=n && v<=n
    transit_ok=(u,v)->2n<u<=3n && 2n<v<=3n
    bike_ok=(u,v)->3n<u<=4n && 3n<v<=4n
    source_layers=[selected,selected .+ n,selected .+ 2n,selected .+ 3n]
    allowed=[walk_ok,car_ok,transit_ok,bike_ok]
    closures=[mode_closure(original,costs,risks,source_layers[i],allowed[i]) for i in 1:4]
    ranks=[path_redundancy_ranks(original,costs,source_layers[i],allowed[i]) for i in 1:4]
    graph=DiGraph(4k); aggcost=spzeros(Int,4k,4k); aggrisk=spzeros(Int,4k,4k); inf=typemax(Int)
    id=(layer,i)->(layer-1)*k+i
    for layer in 1:4, i in 1:k, j in 1:k
        i==j && continue
        ranks[layer][i,j] <= k && continue
        c,r=closures[layer]
        c[i,j]==inf && continue
        add_edge!(graph,id(layer,i),id(layer,j)); aggcost[id(layer,i),id(layer,j)]=c[i,j]; aggrisk[id(layer,i),id(layer,j)]=r[i,j]
    end
    station_set=Set([s for s in meta["stations"] if s in selected])
    possible=collect(Iterators.take(station_set, min(possible_stations,length(station_set))))
    possible_set=Set(possible); decision_arcs=Tuple{Int,Int}[]; groups=Vector{Vector{Tuple{Int,Int}}}()
    for (i,s) in enumerate(selected)
        phys=id(1,i); car=id(2,i); transit=id(3,i); bike=id(4,i)
        if has_edge(original,s,s+n)
            add_edge!(graph,phys,car); aggcost[phys,car]=costs[s,s+n]; aggrisk[phys,car]=risks[s,s+n]
            add_edge!(graph,car,phys); aggcost[car,phys]=costs[s+n,s]; aggrisk[car,phys]=risks[s+n,s]
        end
        if has_edge(original,s,s+2n)
            add_edge!(graph,phys,transit); aggcost[phys,transit]=costs[s,s+2n]; aggrisk[phys,transit]=risks[s,s+2n]
            add_edge!(graph,transit,phys); aggcost[transit,phys]=costs[s+2n,s]; aggrisk[transit,phys]=risks[s+2n,s]
        end
        if s in possible_set && has_edge(original,s,s+3n)
            up=(phys,bike); down=(bike,phys)
            add_edge!(graph,up...); add_edge!(graph,down...)
            aggcost[up...]=costs[s,s+3n]; aggcost[down...]=costs[s+3n,s]
            aggrisk[up...]=risks[s,s+3n]; aggrisk[down...]=risks[s+3n,s]
            # Both directions must be decision arcs; otherwise the return
            # arc is fixed open and the station group cannot be budgeted.
            append!(decision_arcs,(up,down)); push!(groups,[up,down])
        end
    end
    nstations=length(groups)
    budget=ceil(Int,budget_fraction*nstations)
    rng=MersenneTwister(seed); users=User[]
    for u in 1:nusers
        oi,di=rand(rng,1:k),rand(rng,1:k); while oi==di; di=rand(rng,1:k); end
        push!(users,User("U$(u)",id(1,oi),id(1,di),aggrisk,aggcost,nothing,nothing))
    end
    prices=Dict(a=>1 for a in decision_arcs)
    inst=HNDPwC(graph,users,decision_arcs,prices,nothing,groups,budget)
    bike_layer_risks=[aggrisk[u,v] for u in (3k+1):4k for v in (3k+1):4k if has_edge(graph,u,v)]
    validation=(bike_transfer_risk_count=count(a->aggrisk[a...]<0,decision_arcs), bike_layer_risk_count=count(x->x<0,bike_layer_risks))
    metadata2=Dict{String,Any}("name"=>"$(meta["name"])_agg$(k)_p$(length(decision_arcs)÷2)_b$(round(Int,100*budget_fraction))_u$(nusers)", "nnodes"=>4k, "narcs"=>ne(graph), "decision_arcs"=>length(decision_arcs), "decision_budget"=>budget, "bike_possible_stations"=>length(decision_arcs)÷2, "bike_transfer_risk_count"=>validation.bike_transfer_risk_count, "bike_layer_risk_count"=>validation.bike_layer_risk_count, "aggregation_target_nodes"=>target_nodes, "budget_fraction"=>budget_fraction, "nusers"=>nusers)
    HNDPGeneratedNetwork(String(metadata2["name"]),inst,metadata2)
end

function run_case(generated, model_spec, params, case_id)
    t=time(); row=Dict{String,Any}(); status="ok"
    try
        run_root=BENCHMARK_KIND in ("path", "path_subset") ? "tmp_compare/aggregated_bike_expansion_path_runs" : "tmp_compare/aggregated_bike_expansion_runs"
        _,row=_run_hndp_experiment(generated,model_spec,params,run_root,false)
    catch err
        status="error"; row["Error"]=sprint(showerror,err)
    end
    # Preserve solver termination information in the streamed benchmark status.
    # In particular, a memory-limited solve is a valid completed experiment
    # outcome, not an execution error.
    opt_status = get(row, "Opt_status", nothing)
    opt_status == "Memory_Limit" && (status = "memory_limit")
    opt_status == "Timelimit" && (status = "timeout")
    row["case_id"]=case_id; row["wall_time_seconds"]=time()-t; row["run_status"]=status
    row
end

function write_rows(rows)
    isempty(rows) && return
    columns=sort!(collect(Set(vcat([collect(keys(r)) for r in rows]...))))
    table=DataFrame([Symbol(c) => [get(r,c,missing) for r in rows] for c in columns])
    CSV.write(OUT,table)
end

function main()
    rows=Dict{String,Any}[]
    models=if BENCHMARK_KIND in ("path", "path_subset")
        [Dict{String,Any}("name"=>"path_all_accelerations", "model_type"=>"path", "parallelize"=>false, "use_decision_arc_dominance"=>true, "enumeration_time_limit"=>LIMIT)]
    elseif BENCHMARK_KIND == "sdblc"
        [Dict{String,Any}("name"=>"strong_duality_fixed_path", "model_type"=>"sd", "big_m_mode"=>"fixed_network_path", "indicator_constraints"=>false, "bound_duals"=>true), Dict{String,Any}("name"=>"blc_fixed_path", "model_type"=>"blc", "big_m_mode"=>"fixed_network_path", "subproblem_method"=>"mip")]
    else
        throw(ArgumentError("JUBIC_BENCHMARK_KIND must be 'sdblc', 'path', or 'path_subset', got $(BENCHMARK_KIND)."))
    end
    params=Dict{String,Any}("name"=>"gurobi_10min", "mip_solver"=>"Gurobi", "runtime"=>LIMIT, "threads_master"=>8, "threads_sub_con"=>8, "parallel_separation"=>true, "seed"=>SEED)
    cities = BENCHMARK_KIND == "path_subset" ? [city for city in CITIES if city.name in ("cologne", "heidelberg")] : CITIES
    node_targets = BENCHMARK_KIND == "path_subset" ? [100, 200] : NODE_TARGETS
    user_counts = BENCHMARK_KIND == "path_subset" ? [SUBSET_USERS] : USER_COUNTS
    for city in cities
        println("Generating ",city.name)
        g=generate_aachen_osm_gtfs(osm_path=city.osm,boundary_path=city.boundary,use_boundary_compression=true,gtfs_dir=city.gtfs,nusers=1,seed=SEED,bike_station_count=10000,all_nodes_bike_stations=false,generalized_cost=true,bike_entry_cost_cents=50,bike_cost_per_km_cents=5,bike_daily_subscription_cost_cents=nothing)
        n=g.metadata["nnodes"] ÷ 4
        for target in node_targets
            target>n && continue
            for p in STATION_COUNTS
                p>target && continue
                for budget_fraction in BUDGET_FRACTIONS
                    for nusers in user_counts
                        println(city.name," nodes=",target," possible_stations=",p," budget=",budget_fraction," users=",nusers)
                        agg=build_aggregated(g,target,p,budget_fraction,nusers,SEED+target+p+round(Int,100*budget_fraction)+nusers)
                        for model in models
                            case_id="$(city.name)_n$(target)_p$(p)_b$(round(Int,100*budget_fraction))_u$(nusers)_$(model["name"] )"
                            row=run_case(agg,model,params,case_id)
                            row["city"]=city.name; row["target_nodes"]=target; row["possible_stations_requested"]=p; row["budget_fraction"]=budget_fraction; row["nusers"]=nusers
                            append!(rows,[row]); write_rows(rows)
                        end
                    end
                end
            end
        end
    end
    println("Wrote ",length(rows)," rows to ",OUT)
end
if abspath(PROGRAM_FILE)==abspath(@__FILE__)
    main()
end
