module OSMFixedNetworkDemand
using Mmap, JSON, SHA, CSV

"""CSR graph with distinct generalized costs and physical path lengths in km."""
struct FixedGraph
    physical_nodes::Int
    vertices::Int
    offsets::Vector{Int}
    targets::Vector{Int}
    costs::Vector{Int64}
    lengths_km::Vector{Float64}
end
function FixedGraph(n::Int, edges; vertices=maximum((max(e[1],e[2]) for e in edges);init=n))
    rows=sort!(collect(edges),by=e->(e[1],e[2]))
    offsets=ones(Int,vertices+1)
    for (u,v,c,l) in rows
        1<=u<=vertices && 1<=v<=vertices || throw(ArgumentError("Invalid arc endpoint."))
        isinteger(c) && c>=0 || throw(ArgumentError("Fixed costs must be nonnegative integer cents."))
        isfinite(l) && l>=0 || throw(ArgumentError("Physical lengths must be finite, nonnegative kilometres."))
        offsets[u+1]+=1
    end
    for i in 2:length(offsets);offsets[i]+=offsets[i-1]-1;end
    FixedGraph(n,vertices,offsets,Int[e[2] for e in rows],Int64[e[3] for e in rows],Float64[e[4] for e in rows])
end
mutable struct Workspace
    cost::Vector{Int64}
    length::Vector{Float64}
    heap::Vector{Tuple{Int64,Float64,Int}}
end
Workspace(g::FixedGraph)=Workspace(fill(typemax(Int64),g.vertices),fill(Inf,g.vertices),Tuple{Int64,Float64,Int}[])
function heap_push!(heap,x)
    push!(heap,x);i=length(heap)
    while i>1
        p=i÷2;isless(heap[i],heap[p]) || break
        heap[i],heap[p]=heap[p],heap[i];i=p
    end
end
function heap_pop!(heap)
    first=heap[1];last=pop!(heap)
    if !isempty(heap)
        heap[1]=last;i=1
        while true
            a=2i;b=a+1;j=i
            a<=length(heap) && isless(heap[a],heap[j]) && (j=a)
            b<=length(heap) && isless(heap[b],heap[j]) && (j=b)
            j==i && break
            heap[i],heap[j]=heap[j],heap[i];i=j
        end
    end
    first
end
"""Minimize cost; among equal-cost paths minimize physical length."""
function shortest_paths!(w::Workspace,g::FixedGraph,source::Int)
    fill!(w.cost,typemax(Int64));fill!(w.length,Inf);empty!(w.heap)
    w.cost[source]=0;w.length[source]=0;heap_push!(w.heap,(Int64(0),0.0,source))
    while !isempty(w.heap)
        cost,length,u=heap_pop!(w.heap)
        cost==w.cost[u] && length==w.length[u] || continue
        for e in g.offsets[u]:g.offsets[u+1]-1
            v=g.targets[e];newcost=Base.Checked.checked_add(cost,g.costs[e]);newlength=length+g.lengths_km[e]
            if newcost<w.cost[v] || (newcost==w.cost[v] && newlength<w.length[v])
                w.cost[v]=newcost;w.length[v]=newlength
                heap_push!(w.heap,(newcost,newlength,v))
            end
        end
    end
    w
end
function matrix_file(path,::Type{T},n;write=false) where T
    io=open(path,write ? "w+" : "r")
    write && truncate(io,n*n*sizeof(T))
    matrix=Mmap.mmap(io,Matrix{T},(n,n));close(io);matrix
end
"""
Persist exact fixed-network shortest cost and selected-path physical length.
Disk matrices are indexed [destination, origin], one contiguous source column.
Only one Dijkstra workspace is held in memory. No geometric detour factor.
"""
function cache_all_pairs(g::FixedGraph,folder;progress=(done,total)->nothing)
    mkpath(folder);n=g.physical_nodes
    marker=joinpath(folder,"complete.json")
    isfile(marker) && rm(marker)
    costs=matrix_file(joinpath(folder,"costs.i64"),Int64,n;write=true)
    lengths=matrix_file(joinpath(folder,"lengths_km.f64"),Float64,n;write=true)
    done=Threads.Atomic{Int}(0);workers=min(Threads.nthreads(),n)
    Threads.@threads :static for worker in 1:workers
        w=Workspace(g)
        for origin in worker:workers:n
            shortest_paths!(w,g,origin)
            for destination in 1:n
                isfinite(w.length[destination]) || throw(ArgumentError("Unreachable fixed-network OD $origin → $destination."))
                costs[destination,origin]=w.cost[destination]
                lengths[destination,origin]=w.length[destination]
            end
            progress(Threads.atomic_add!(done,1)+1,n)
        end
    end
    Mmap.sync!(costs);Mmap.sync!(lengths)
    open(marker,"w") do io
        JSON.print(io,Dict("physical_nodes"=>n,"orientation"=>"destination,origin","cost_type"=>"Int64","length_type"=>"Float64","path_rule"=>"minimum cost, then minimum physical length","graph_sha256"=>fingerprint(g),"schema"=>1),2)
    end
    costs,lengths
end
function fingerprint(g::FixedGraph)
    io=IOBuffer();write(io,g.physical_nodes,g.vertices)
    for column in (g.offsets,g.targets,g.costs,g.lengths_km);write(io,column);end
    bytes2hex(sha256(take!(io)))
end
function load_cache(folder;expected_graph_hash=nothing)
    meta=JSON.parsefile(joinpath(folder,"complete.json"));n=Int(meta["physical_nodes"])
    expected_graph_hash===nothing || get(meta,"graph_sha256",nothing)==expected_graph_hash || throw(ArgumentError("Stale fixed-network path cache."))
    matrix_file(joinpath(folder,"costs.i64"),Int64,n),matrix_file(joinpath(folder,"lengths_km.f64"),Float64,n)
end
"""
Access uses origin→site cost; egress uses site→destination cost. Directed
asymmetry is respected. Equal costs are resolved by physical length, then ID.
"""
function nearest_sites(costs,lengths,sites)
    n=size(costs,1);access=zeros(Int,n);egress=zeros(Int,n)
    isempty(sites) && throw(ArgumentError("At least one retained location is required."))
    all(1<=s<=n for s in sites) && length(unique(sites))==length(sites) || throw(ArgumentError("Invalid retained locations."))
    for node in 1:n
        access[node]=argmin(i->(costs[sites[i],node],lengths[sites[i],node],sites[i]),eachindex(sites))
        egress[node]=argmin(i->(costs[node,sites[i]],lengths[node,sites[i]],sites[i]),eachindex(sites))
    end
    access,egress
end
"""
Construct full-network band counts, then sum original-pair mass into retained
ODs. Diagonal mass is retained separately, without renormalizing it away.
"""
function aggregate_demand(costs,lengths,sites,bounds,mass_units;
    car_probability::Rational=71//100,progress=(done,total)->nothing)
    n=size(costs,1);k=length(sites);nb=length(bounds)+1
    length(mass_units)==nb || throw(ArgumentError("Distance-band masses mismatch."))
    issorted(bounds) || throw(ArgumentError("Distance-band bounds must be sorted."))
    all(x->isfinite(x)&&x>=0,bounds) || throw(ArgumentError("Distance bounds must be finite and nonnegative."))
    all(x->isinteger(x)&&x>=0,mass_units) && sum(mass_units)>0 || throw(ArgumentError("Demand masses must be nonnegative integer units with positive total."))
    0<=car_probability<=1 || throw(ArgumentError("Invalid car probability."))
    access,egress=nearest_sites(costs,lengths,sites)
    counts=zeros(Int64,nb);hist=zeros(Int64,k,k,nb)
    for origin in 1:n
        for destination in 1:n
            origin==destination && continue
            b=searchsortedfirst(bounds,lengths[destination,origin])
            counts[b]+=1;hist[access[origin],egress[destination],b]+=1
        end
        progress(origin,n)
    end
    all(counts[b]>0 || mass_units[b]==0 for b in 1:nb) || throw(ArgumentError("A positive target band has no original OD support; consult its definition."))
    q=BigInt.(mass_units).//BigInt(sum(mass_units));total=zeros(Rational{BigInt},k,k)
    for b in 1:nb,j in 1:k,i in 1:k
        counts[b]==0 && continue
        total[i,j]+=q[b]*(BigInt(hist[i,j,b])//BigInt(counts[b]))
    end
    p=BigInt(numerator(car_probability))//BigInt(denominator(car_probability))
    owner=total.*p;nonowner=total.*(1-p)
    @assert sum(total)==1 && sum(owner)==p && sum(nonowner)==1-p
    diagonal=sum(total[i,i] for i in 1:k)
    (;total,owner,nonowner,band_counts=counts,band_aggregate_pair_counts=hist,
      access,egress,diagonal_mass=diagonal,inter_site_mass=1-diagonal)
end
"""Load the fixed walk/car/transit graph with explicit length provenance.

All bike-layer arcs are omitted: with sharing transfers closed that layer is
unreachable from physical origins. Non-car profiles also omit car-layer arcs.
Transfers between colocated layers have physical length zero. No missing
travel length is silently replaced by geometric endpoint distance here.
"""
function load_fixed_graph(network_folder,length_paths;car_available=true,car_override=nothing)
    meta=JSON.parsefile(joinpath(network_folder,"metadata.json"));n=Int(meta["nnodes"])÷4
    lengths=Dict{Tuple{Int,Int},Float64}();sources=Dict{Tuple{Int,Int},String}()
    for path in length_paths,row in CSV.File(path)
        arc=(Int(row.source),Int(row.target));value=Float64(row.length_km)
        haskey(lengths,arc) && throw(ArgumentError("Duplicate physical length for $arc."))
        lengths[arc]=value;sources[arc]=String(row.length_source)
    end
    overrides=Dict{Tuple{Int,Int},Int64}()
    if car_override!==nothing
        for row in CSV.File(car_override);overrides[(Int(row.source),Int(row.target))]=Int64(row.cost_cents);end
    end
    edges=Tuple{Int,Int,Int64,Float64}[];provenance=Dict{String,Int}()
    for row in CSV.File(joinpath(network_folder,"edges.csv"))
        u=Int(row.source);v=Int(row.target)
        u<=3n && v<=3n || continue
        !car_available && (n<u<=2n || n<v<=2n) && continue
        arc=(u,v);cost=get(overrides,arc,Int64(row.cost_cents))
        transfer=cld(u,n)!=cld(v,n)
        if transfer
            mod1(u,n)==mod1(v,n) || throw(ArgumentError("Noncolocated transfer needs an explicit physical length."))
            length=0.0;source="colocated_transfer"
        else
            haskey(lengths,arc) || throw(ArgumentError("Missing physical arc length for $arc."))
            length=lengths[arc];source=sources[arc]
        end
        push!(edges,(u,v,cost,length));provenance[source]=get(provenance,source,0)+1
    end
    FixedGraph(n,edges;vertices=3n),provenance
end
end
