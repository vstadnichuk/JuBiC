using CSV,JSON,Graphs,SparseArrays,SHA,Mmap,Dates,Base.Threads
driver=joinpath(@__DIR__,"full_network_city_profiles.jl")
include_string(Main,replace(read(driver,String),r"\nmain\(\)\s*$"=>"\n"),driver)
include(joinpath(SRC, "demand", "osm_aggregation_mode_share_analysis.jl"))
const CONFIG=JSON.parsefile(SETTINGS_FILE)
const TARGETS=Int.(CONFIG["targets"])
const CACHE=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(B,"work","fixed_reference_cache")

function preferred_family(bits)
    masks=[m for m in 0:15 if bits&(1<<m)!=0]
    isempty(masks) && error("Unreachable route")
    chosen=argmin(m->(m&8!=0 ? 0 : m&2!=0 ? 1 : m&4!=0 ? 2 : 3,count_ones(m&14),m&14,m),masks)
    p=chosen&14
    p==0 ? 1 : p==2 ? 2 : p==4 ? 3 : p==8 ? 4 : 5
end
function sharing_subset(g,allowed)
    n=g.physical_nodes
    arcs=Tuple{Int,Int,Int64}[]
    for u in 1:4n,e in g.offsets[u]:g.offsets[u+1]-1
        v=g.targets[e]
        bike_transfer=(u<=n && v==u+3n)||(v<=n && u==v+3n)
        bike_transfer && !(min(u,v) in allowed) && continue
        push!(arcs,(u,v,g.costs[e]))
    end
    MS.mode_graph(n,arcs)
end
function evaluate_aggregate(g,conditional_matrices,p)
    shares=zeros(Rational{BigInt},5);diagonal=zero(p)
    for (cl,available,cm) in [("car",true,p),("noncar",false,1-p)]
        graph=available ? g : MS.without_car(g);demand=conditional_matrices[cl]
        diagonal+=cm*sum(demand[i,i] for i in axes(demand,1))
        partial=[zeros(Rational{BigInt},5) for _ in 1:Threads.nthreads()]
        Threads.@threads :static for worker in 1:Threads.nthreads()
            w=MS.Workspace(graph);s=partial[worker]
            for o in worker:Threads.nthreads():g.physical_nodes
                MS.shortest_modes!(w,graph,o)
                for d in 1:g.physical_nodes
                    o==d && continue
                    f=preferred_family(Int(w.labels[d]));s[f]+=cm*demand[o,d]
                end
            end
        end
        shares .+= reduce(+,partial)
    end
    @assert sum(shares)+diagonal==1
    shares,diagonal
end
function aggregate_city(city)
    started=time();root=joinpath(B,"inputs/city_networks_walk_connected")
    folder=joinpath(B,"inputs/aggregation",city);mkpath(folder)
    demandfolder=joinpath(B,"inputs/demand",city);mkpath(demandfolder)
    meta=JSON.parsefile(joinpath(root,"graphs",city,"metadata.json"))
    g,_=MS.load_city(root,city);g=MS.with_car_costs(g,joinpath(root,"speed_overrides",city*".csv"))
    n=g.physical_nodes;stations=MS.city_candidate_stations(root,city,meta)
    fg=DiGraph(4n);fc=spzeros(Int,4n,4n)
    for u in 1:4n,e in g.offsets[u]:g.offsets[u+1]-1
        v=g.targets[e];add_edge!(fg,u,v);fc[u,v]=g.costs[e]
    end
    println(city," FULL MULTIMODAL SELECTION ",length(stations)," eligible sites");flush(stdout)
    metric=ExistingAggregation.multimodal_site_metric(fg,fc,stations)
    selection=NamedTuple[];ids=Dict{Int,Vector{Int}}();submetrics=Dict()
    for k in TARGETS
        @assert k<=length(stations)
        medoids,sizes=ExistingAggregation.balanced_medoids(metric,k);ids[k]=stations[medoids]
        submetrics[k]=metric[medoids,medoids]
        ExistingAggregation.save_selection_metric(folder,k,submetrics[k])
        for (i,s) in enumerate(ids[k])
            push!(selection,(;city,target=k,aggregate_node=i,original_node=s,latitude=meta["coords"][s][1],longitude=meta["coords"][s][2],cluster_size=sizes[i]))
        end
    end
    CSV.write(joinpath(folder,"selected_nodes.csv"),selection)
    profile=only(filter(x->x["city"]==city,JSON.parsefile(joinpath(PROVENANCE,"osm_distance_demand_sources.json"))["cities"]))
    p=frac_local(JSON.parsefile(joinpath(B,"raw/full_network_city_profiles",city,"result.json"))["car_probability_exact"])
    matrices=Dict(k=>Dict{String,Matrix{Rational{BigInt}}}() for k in TARGETS)
    demandmeta=Dict()
    for cl in ["car","noncar"]
        costs,lengths=FD.load_cache(cl=="noncar" ? joinpath(NONCAR_CACHE,city,"paths_noncar") : joinpath(CACHE,city,"paths_car"))
        for k in TARGETS
            println(city," DEMAND ",cl," ",k);flush(stdout)
            r=FD.aggregate_demand(costs,lengths,ids[k],Float64.(profile["upper_bounds_km"]),Int.(profile["mass_units"]);car_probability=1//1)
            matrices[k][cl]=r.total
            demandmeta["$(cl)_$(k)"]=Dict("band_counts"=>r.band_counts,"conditional_diagonal_exact"=>string(r.diagonal_mass))
            CSV.write(joinpath(demandfolder,"mapping_$(cl)_$(k).csv"),((physical_node=i,access_site=r.access[i],egress_site=r.egress[i]) for i in 1:n))
        end
        costs=nothing;lengths=nothing;GC.gc()
    end
    graphs,costs=split_layers(g);rows=NamedTuple[];candidate_rows=NamedTuple[]
    for k in TARGETS
        println(city," COMPACT ",k);flush(stdout)
        ag=aggregate_existing(g,ids[k],graphs,costs)
        CSV.write(joinpath(folder,"graph_$(k)_edges.csv"),((source=u,target=ag.targets[e],cost_cents=ag.costs[e]) for u in 1:4k for e in ag.offsets[u]:ag.offsets[u+1]-1))
        owner=p.*matrices[k]["car"];nonowner=(1-p).*matrices[k]["noncar"]
        CSV.write(joinpath(demandfolder,"demand_$(k).csv"),((origin=i,destination=j,
            car_numerator=string(numerator(owner[i,j])),car_denominator=string(denominator(owner[i,j])),
            noncar_numerator=string(numerator(nonowner[i,j])),noncar_denominator=string(denominator(nonowner[i,j])),
            car_demand=Float64(owner[i,j]),noncar_demand=Float64(nonowner[i,j]),intrazonal=i==j) for j in 1:k for i in 1:k))
        for m in unique(vcat(filter(x->x<=k,Int.(CONFIG["candidate_counts"])),k))
            order=ExistingAggregation.potential_station_sequence(submetrics[k],m)
            for (rank,i) in enumerate(order)
                push!(candidate_rows,(;city,locations=k,candidates=m,rank,aggregate_node=i,original_node=ids[k][i]))
            end
            shares,diagonal=evaluate_aggregate(sharing_subset(ag,Set(order)),matrices[k],p)
            push!(rows,(;city,locations=k,candidates_open=m,scenario=m==k ? "all_retained_eligible_open" : "potential_pool_all_open",
                diagonal_pct=100Float64(diagonal),interlocation_mass_exact=string(1-diagonal),
                (Symbol(NAMES[f]*"_original_pct")=>100Float64(shares[f]) for f in 1:5)...,
                (Symbol(NAMES[f]*"_conditional_pct")=>100Float64(shares[f]/(1-diagonal)) for f in 1:5)...,
                (Symbol(NAMES[f]*"_exact")=>string(shares[f]) for f in 1:5)...))
            println(city," MODAL ",k," candidates=",m," ",[round(100Float64(x/(1-diagonal)),digits=3) for x in shares]);flush(stdout)
        end
        CSV.write(joinpath(B,"data",city*"_aggregation_modal.csv"),rows)
        CSV.write(joinpath(folder,"candidate_sets.csv"),candidate_rows)
    end
    ref=Dict("policy"=>CONFIG["station_selection_reference"],"candidate_order"=>CONFIG["candidate_order"],
        "network_sha256"=>bytes2hex(sha256(read(joinpath(root,"graphs",city,"edges.csv")))),"eligible_sites"=>stations,
        "targets"=>TARGETS,"seconds"=>time()-started,"demand_details"=>demandmeta,
        "intrazonal_policy"=>"retain diagonal; do not classify as walking or redistribute; conditional shares exclude it",
        "same_existing_candidate_universe"=>true,"no_MIP_solved"=>true)
    open(joinpath(folder,"selection_reference.json"),"w") do io;JSON.print(io,ref,2);end
    println("FINISHED AGGREGATES ",city," ",time()-started,"s");flush(stdout)
end
frac_local(x)=parse(BigInt,split(x,"//")[1])//parse(BigInt,split(x,"//")[2])

const NONCAR_CACHE=get(ENV,"JUBIC_NONCAR_CACHE",CACHE)
# Arguments: [cache_folder [city ...]]
for city in (length(ARGS) > 1 ? ARGS[2:end] : CONFIG["cities"])
    marker=joinpath(B,"inputs/aggregation",city,"COMPLETE.json")
    isfile(marker) && (println("CACHED PREPARATION ",city);continue)
    root=joinpath(B,"inputs/city_networks_walk_connected")
    for (cl,available) in [("car",true),("noncar",false)]
        fixed,_=FD.load_fixed_graph(joinpath(root,"graphs",city),[joinpath(B,"inputs/geometry/road_lengths",city*".csv"),joinpath(B,"inputs/geometry/transit_lengths",city*".csv")];car_available=available,car_override=joinpath(root,"speed_overrides",city*".csv"))
        hash=FD.fingerprint(fixed);cache=available ? joinpath(CACHE,city,"paths_car") : joinpath(NONCAR_CACHE,city,"paths_noncar")
        if !isfile(joinpath(cache,"complete.json"))
            println("BUILDING REFERENCE ",city," ",cl);flush(stdout)
            FD.cache_all_pairs(fixed,cache;progress=(d,n)->begin
                if d%1000==0 || d==n;println(city," fixed ",cl," ",d,"/",n);flush(stdout);end
            end)
        end
        @assert JSON.parsefile(joinpath(cache,"complete.json"))["graph_sha256"]==hash
        fixed=nothing;GC.gc()
    end
    aggregate_city(city)
    open(marker,"w") do io;JSON.print(io,Dict("city"=>city,"complete"=>true,"car_eur_per_km"=>0.22,"targets"=>TARGETS,"selection_reference"=>CONFIG["selection_reference"]),2);end
    GC.gc()
end
