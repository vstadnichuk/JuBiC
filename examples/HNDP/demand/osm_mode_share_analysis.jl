module OSMModeShares
using CSV, JSON
using Base.Threads

struct ModeGraph
    physical_nodes::Int
    offsets::Vector{Int}
    targets::Vector{Int}
    costs::Vector{Int64}
    flags::Vector{UInt8}
end

"""Build compact adjacency arrays. Zero-cost arcs must be transit exits.

Travel arcs mark walking/car/transit/bike as bits 1/2/4/8. Transfers do not
count as travel. This permits reporting walking access separately from trips
combining multiple powered modes. Integer costs preserve exact ties.
"""
function mode_graph(n::Int, arcs)
    n>0 || throw(ArgumentError("physical node count must be positive"))
    ordered=sort!(collect(arcs),by=e->(e[1],e[2]))
    offsets=ones(Int,4n+1); targets=Int[]; costs=Int64[]; flags=UInt8[]
    counts=zeros(Int,4n)
    seen=Set{Tuple{Int,Int}}()
    for (u,v,w) in ordered
        1<=u<=4n && 1<=v<=4n || throw(ArgumentError("invalid vertex"))
        w>=0 && isinteger(w) && w<typemax(Int64)÷(4n) || throw(ArgumentError("costs must be bounded nonnegative integers"))
        (u,v) in seen && throw(ArgumentError("duplicate arc ($u,$v)"))
        push!(seen,(u,v))
        lu,lv=cld(u,n),cld(v,n)
        w==0 && !(lu==3 && lv==1) && throw(ArgumentError("zero-cost arcs must be transit-to-physical exits"))
        push!(targets,v);push!(costs,Int64(w));counts[u]+=1
        push!(flags,lu==lv ? UInt8(1<<(lu-1)) : UInt8(0))
    end
    for u in 1:4n;offsets[u+1]=offsets[u]+counts[u];end
    ModeGraph(n,offsets,targets,costs,flags)
end

# A UInt16 records which of the 16 mode combinations occur on optimal paths.
# Apply one travel mode to every represented combination without path copying.
@inline function add_mode(s::UInt16,f::UInt8)
    f==0x01 && return ((s&0x5555)<<1)|(s&0xaaaa)
    f==0x02 && return ((s&0x3333)<<2)|(s&0xcccc)
    f==0x04 && return ((s&0x0f0f)<<4)|(s&0xf0f0)
    f==0x08 && return ((s&0x00ff)<<8)|(s&0xff00)
    return s
end

mutable struct Workspace
    distance::Vector{Int64}
    labels::Vector{UInt16}
    heap::Vector{Int}
    position::Vector{Int}
end
Workspace(g::ModeGraph)=Workspace(fill(typemax(Int64),4g.physical_nodes),zeros(UInt16,4g.physical_nodes),Int[],zeros(Int,4g.physical_nodes))

# At equal distance, transit vertices precede physical vertices. Together with
# positive travel/entry costs, this is a topological order of shortest-path
# arcs, including the zero-cost transit exits. Every tied predecessor is thus
# processed before its successor; no path-mask update is lost after settlement.
@inline function before(a,b,w::Workspace,n)
    da,db=w.distance[a],w.distance[b]
    da!=db && return da<db
    ta,tb=2n<a<=3n,2n<b<=3n
    ta!=tb && return ta
    return a<b
end
@inline function rise!(w::Workspace,v,n)
    i=w.position[v]
    if i==0
        push!(w.heap,v);i=length(w.heap)
    end
    while i>1
        p=i>>1;parent=w.heap[p]
        before(v,parent,w,n) || break
        w.heap[i]=parent;w.position[parent]=i;i=p
    end
    w.heap[i]=v;w.position[v]=i
end
@inline function pop_vertex!(w::Workspace,n)
    v=w.heap[1];last=pop!(w.heap);w.position[v]=-1
    if !isempty(w.heap)
        i=1;len=length(w.heap)
        while 2i<=len
            child=2i
            child<len && before(w.heap[child+1],w.heap[child],w,n) && (child+=1)
            before(w.heap[child],last,w,n) || break
            w.heap[i]=w.heap[child];w.position[w.heap[i]]=i;i=child
        end
        w.heap[i]=last;w.position[last]=i
    end
    return v
end

function shortest_modes!(w::Workspace,g::ModeGraph,source::Int)
    n=g.physical_nodes
    fill!(w.distance,typemax(Int64));fill!(w.labels,0);fill!(w.position,0);empty!(w.heap)
    w.distance[source]=0;w.labels[source]=1;rise!(w,source,n)
    @inbounds while !isempty(w.heap)
        u=pop_vertex!(w,n);du=w.distance[u];lu=w.labels[u]
        for e in g.offsets[u]:g.offsets[u+1]-1
            v=g.targets[e];candidate=du+g.costs[e]
            if candidate<w.distance[v]
                w.distance[v]=candidate;w.labels[v]=add_mode(lu,g.flags[e]);rise!(w,v,n)
            elseif candidate==w.distance[v]
                w.labels[v] |= add_mode(lu,g.flags[e])
            end
        end
    end
    return w
end

"""Count all ordered physical-node pairs, excluding self-pairs.

One indexed-heap Dijkstra per origin; sources run in independent worker chunks.
Memory is O(workers*(V+2^16)+E), with no dense all-pairs distance matrix.
"""
function count_mode_shares(g::ModeGraph; workers=Threads.nthreads(), progress=true)
    n=g.physical_nodes;workers=min(workers,Threads.nthreads(),n)
    workers>0 || throw(ArgumentError("workers must be positive"))
    counters=[zeros(Int64,65536) for _ in 1:workers]
    completed=Threads.Atomic{Int}(0);started=time()
    Threads.@threads :static for worker in 1:workers
        w=Workspace(g);counts=counters[worker]
        for source in worker:workers:n
            shortest_modes!(w,g,source)
            @inbounds for target in 1:n
                source==target && continue
                counts[Int(w.labels[target])+1]+=1
            end
            done=Threads.atomic_add!(completed,1)+1
            if progress && (done%1000==0 || done==n)
                println("  origins ",done,"/",n,"; ",round(time()-started,digits=1)," seconds");flush(stdout)
            end
        end
    end
    return reduce(+,counters)
end

haversine(lat1,lon1,lat2,lon2)=2*6371.0088*asin(min(1.0,sqrt(sind((lat2-lat1)/2)^2+cosd(lat1)*cosd(lat2)*sind((lon2-lon1)/2)^2)))

"""Same nearest-road-node / 250m rule as the importer, using spatial buckets."""
function candidate_stations(coords,stops_path; radius_km=0.25)
    width=0.005
    bucket(lat,lon)=(floor(Int,lat/width),floor(Int,lon/width))
    grid=Dict{Tuple{Int,Int},Vector{Int}}()
    for (i,p) in enumerate(coords);push!(get!(grid,bucket(p[1],p[2]),Int[]),i);end
    stations=Set{Int}()
    for row in CSV.File(stops_path;silencewarnings=true)
        (ismissing(row.stop_lat)||ismissing(row.stop_lon)) && continue
        lat,lon=Float64(row.stop_lat),Float64(row.stop_lon)
        isfinite(lat)&&isfinite(lon) || continue
        a,b=bucket(lat,lon)
        # Conservative spherical longitude bound; only nodes within radius
        # matter. Ties use the first coordinate, matching the linear importer.
        delta_lat=rad2deg(radius_km/6371.0088)
        lat_span=ceil(Int,delta_lat/width)+1
        lon_span=ceil(Int,delta_lat/(max(0.001,cosd(abs(lat)+delta_lat))*width))+1
        best=0;bestdist=Inf
        for x in a-lat_span:a+lat_span,y in b-lon_span:b+lon_span,i in get(grid,(x,y),Int[])
            p=coords[i];d=haversine(lat,lon,p[1],p[2])
            if d<bestdist || (d==bestdist && i<best)
                best=i;bestdist=d
            end
        end
        bestdist<=radius_km && push!(stations,best)
    end
    return sort!(collect(stations))
end

function city_candidate_stations(bundle,city,meta)
    if haskey(meta,"all_candidate_stations")
        stations=sort!(unique(Int.(meta["all_candidate_stations"])))
        all(s->1<=s<=length(meta["coords"]),stations) || error("invalid retained candidate station")
        return stations
    end
    path=joinpath(bundle,"sources",city,"gtfs","stops.txt")
    isfile(path) || (path=joinpath(meta["gtfs_dir"],"stops.txt"))
    candidate_stations(meta["coords"],path)
end

function load_city(bundle,city; station_scope="gtfs_candidates")
    folder=joinpath(bundle,"graphs",city)
    meta=JSON.parsefile(joinpath(folder,"metadata.json"))
    n=Int(meta["nnodes"])÷4
    coords=meta["coords"]
    length(coords)==n || error("coordinate count mismatch")
    stations=station_scope=="all_nodes" ? collect(1:n) : station_scope=="gtfs_candidates" ?
        city_candidate_stations(bundle,city,meta) :
        throw(ArgumentError("station_scope must be gtfs_candidates or all_nodes"))
    Set(Int.(meta["stations"])) ⊆ Set(stations) || error("candidate set omitted an existing station")
    arcs=Dict{Tuple{Int,Int},Int64}()
    for e in CSV.File(joinpath(folder,"edges.csv"))
        isinteger(e.cost_cents) || error("nonintegral cost")
        arcs[(Int(e.source),Int(e.target))]=Int64(e.cost_cents)
    end
    time_cost(seconds)=max(1,round(Int,100*seconds*meta["bike_value_of_time_eur_per_hour"]/3600))
    fee=meta["bike_daily_subscription_cost_cents"]===nothing ? meta["bike_entry_cost_cents"] : meta["bike_daily_subscription_cost_cents"]
    entry=time_cost(meta["bike_access_time_seconds"])+Int(fee)
    exit=time_cost(meta["bike_return_time_seconds"])
    for s in stations
        arcs[(s,s+3n)]=entry;arcs[(s+3n,s)]=exit
    end
    graph=mode_graph(n,[(u,v,w) for ((u,v),w) in arcs])
    info=Dict{String,Any}("city"=>city,"station_scope"=>station_scope,"available_bike_stations"=>length(stations),
        "previous_bike_stations"=>length(meta["stations"]),"station_budget_enforced"=>false,
        "physical_nodes"=>n,"layered_nodes"=>4n,"edges"=>length(arcs),"bike_entry_cost_cents"=>entry,
        "bike_exit_cost_cents"=>exit,"buffer_km"=>meta["buffer_km"],"source_commit"=>meta["source_commit"],
        "source_metadata_file"=>abspath(joinpath(folder,"metadata.json")))
    return graph,info
end

"""Remove all travel/transfer arcs touching the car layer; keep vertex IDs."""
function without_car(g::ModeGraph)
    n=g.physical_nodes
    arcs=Tuple{Int,Int,Int64}[]
    for u in 1:4n,e in g.offsets[u]:g.offsets[u+1]-1
        v=g.targets[e]
        (n<u<=2n || n<v<=2n) && continue
        push!(arcs,(u,v,g.costs[e]))
    end
    return mode_graph(n,arcs)
end

"""Apply a complete, topology-preserving car-cost overlay from the importer."""
function with_car_costs(g::ModeGraph,path)
    n=g.physical_nodes
    replacements=Dict{Tuple{Int,Int},Int64}()
    for row in CSV.File(path)
        u,v=Int(row.source),Int(row.target)
        n<u<=2n && n<v<=2n || throw(ArgumentError("override outside car travel layer"))
        haskey(replacements,(u,v)) && throw(ArgumentError("duplicate override"))
        replacements[(u,v)]=Int64(row.cost_cents)
    end
    arcs=Tuple{Int,Int,Int64}[];used=0
    for u in 1:4n,e in g.offsets[u]:g.offsets[u+1]-1
        v=g.targets[e];cost=g.costs[e]
        if n<u<=2n && n<v<=2n
            haskey(replacements,(u,v)) || throw(ArgumentError("missing car-cost override"))
            cost=replacements[(u,v)];used+=1
        end
        push!(arcs,(u,v,cost))
    end
    used==length(replacements) || throw(ArgumentError("override changes topology"))
    return mode_graph(n,arcs)
end

"""Exact weights for the two copies of EVERY OD pair, with no random sampling.

An owner copy has weight numerator(p)/denominator(p), and a non-owner copy
has the complementary weight. Integer weight units retain this ratio exactly,
including when the number of OD pairs is not divisible by the denominator.
"""
function ownership_mix(owner_counts,nonowner_counts,p::Rational)
    0<=p<=1 || throw(ArgumentError("ownership share must be between zero and one"))
    length(owner_counts)==length(nonowner_counts)==65536 || throw(ArgumentError("invalid histogram size"))
    sum(owner_counts)==sum(nonowner_counts) || throw(ArgumentError("class OD counts differ"))
    for i in eachindex(nonowner_counts)
        nonowner_counts[i]>0 && ((i-1)&0xcccc)!=0 && throw(ArgumentError("non-owner path uses a car"))
    end
    a,b=numerator(p),denominator(p)
    weighted=a .* owner_counts .+ (b-a) .* nonowner_counts
    total=sum(owner_counts)
    @assert sum(weighted)==b*total
    return (weight_units=weighted,weight_denominator=b,owner_units=a*total,
        nonowner_units=(b-a)*total,owner_fraction=p)
end

const FAMILY_NAMES=["walking","car","transit","bike","multiple_powered_modes","tie_between_families"]
const MODE_NAMES=["walking","car","transit","bike"]
mode_name(mask)=join([MODE_NAMES[k] for k in 1:4 if mask&(1<<(k-1))!=0],"+")
function summarize(counts,info)
    total=sum(counts);unreachable=counts[1];reachable=total-unreachable
    families=zeros(Int64,6);patterns=zeros(Int64,16)
    any_multimodal=0;only_multimodal=0;any_powered_mix=0;only_powered_mix=0
    exact_mode_ties=0
    for bits in 1:65535
        c=counts[bits+1];c==0 && continue
        fam=Set{Int}();multi=Bool[];powered=Bool[]
        for mask in 0:15
            bits&(1<<mask)==0 && continue
            patterns[mask+1]+=c
            p=mask>>1
            family=p==0 ? 1 : p==1 ? 2 : p==2 ? 3 : p==4 ? 4 : 5
            push!(fam,family);push!(multi,count_ones(mask)>=2);push!(powered,count_ones(p)>=2)
        end
        families[length(fam)==1 ? only(fam) : 6]+=c
        any_multimodal+=any(multi)*c;only_multimodal+=all(multi)*c
        any_powered_mix+=any(powered)*c;only_powered_mix+=all(powered)*c
        exact_mode_ties+=(count_ones(bits)>1)*c
    end
    out=copy(info)
    out["ordered_pairs"]=total;out["reachable_pairs"]=reachable;out["unreachable_pairs"]=unreachable
    out["family_counts"]=Dict(FAMILY_NAMES[i]=>families[i] for i in 1:6)
    out["family_shares_all_pairs"]=Dict(FAMILY_NAMES[i]=>families[i]/total for i in 1:6)
    out["optimal_mode_combination_pair_counts"]=Dict(mode_name(i)=>patterns[i+1] for i in 1:15)
    out["any_multimodal_optimum_including_walking"]=any_multimodal
    out["only_multimodal_optima_including_walking"]=only_multimodal
    out["any_multiple_powered_modes_optimum"]=any_powered_mix
    out["only_multiple_powered_modes_optima"]=only_powered_mix
    out["pairs_with_tied_exact_mode_combinations"]=exact_mode_ties
    @assert sum(families)+unreachable==total
    return out
end

summary_row(s) = (;city=s["city"],physical_nodes=s["physical_nodes"],
            bike_stations=s["available_bike_stations"],ordered_pairs=s["ordered_pairs"],unreachable_pairs=s["unreachable_pairs"],
            (Symbol(k*"_share")=>s["family_shares_all_pairs"][k] for k in FAMILY_NAMES)...,
            any_multimodal_share=s["any_multimodal_optimum_including_walking"]/s["ordered_pairs"],
            any_multiple_powered_modes_share=s["any_multiple_powered_modes_optimum"]/s["ordered_pairs"],
            seconds=s["analysis_seconds"])

function run(bundle,output;cities=["aachen","cologne","heidelberg","karlsruhe","nuremberg","dresden"],station_scope="gtfs_candidates")
    mkpath(output);summaries=Dict{String,Any}[]
    for city in cities
        println("Loading ",city," (",station_scope,")");flush(stdout)
        graph,info=load_city(bundle,city;station_scope)
        println("  ",graph.physical_nodes," physical nodes; ",info["available_bike_stations"]," bike stations");flush(stdout)
        started=time();counts=count_mode_shares(graph)
        info["analysis_seconds"]=time()-started;info["threads"]=Threads.nthreads()
        summary=summarize(counts,info);push!(summaries,summary)
        open(joinpath(output,city*".json"),"w") do io;JSON.print(io,summary,2);end
        CSV.write(joinpath(output,city*"_optimal_mode_sets.csv"),
            ((optimal_mode_set_bits=i-1,pairs=counts[i]) for i in eachindex(counts) if counts[i]>0))
        CSV.write(joinpath(output,"summary.csv"),summary_row.(summaries))
        println("Finished ",city," in ",round(info["analysis_seconds"],digits=1)," seconds");flush(stdout)
    end
    return summaries
end

end # module

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)>=2 || error("Usage: julia --project=. --threads=4 osm_mode_share_analysis.jl <city_bundle> <output_dir> [city ...]")
    OSMModeShares.run(ARGS[1],ARGS[2];cities=length(ARGS)>2 ? ARGS[3:end] : ["aachen","cologne","heidelberg","karlsruhe","nuremberg","dresden"])
end
