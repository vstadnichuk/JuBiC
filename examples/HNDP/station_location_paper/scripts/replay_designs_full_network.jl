# Station-based rerun: replays the designs in raw/aggregation_cases on the full networks.
using CSV,JSON,SHA,Mmap,Dates,Base.Threads
include(joinpath(@__DIR__,"paths.jl"))
include(joinpath(SRC, "demand", "osm_mode_share_analysis.jl"));include(joinpath(SRC, "demand", "osm_fixed_network_demand.jl"))
const MS=OSMModeShares;const FD=OSMFixedNetworkDemand
const CAR_CACHE=get(ENV,"JUBIC_CAR_CACHE",joinpath(B,"work/fixed_reference_cache"))
const NONCAR_CACHE=get(ENV,"JUBIC_NONCAR_CACHE",joinpath(B,"work/fixed_reference_cache"))
config=JSON.parsefile(SETTINGS_FILE);cities=String.(config["cities"])
complete(city,k)=let f=joinpath(B,"raw/aggregation_cases","$(city)_$(k)_s20_b10/result.json");isfile(f)&&get(JSON.parsefile(f),"status","")=="OPTIMAL";end
common=filter(k->all(complete(c,k) for c in cities),Int.(config["targets"]))
if isempty(common)
    mkpath(joinpath(B,"raw/station_full_network_design_replay"))
    open(joinpath(B,"raw/station_full_network_design_replay/COMPLETE.json"),"w") do io;JSON.print(io,Dict("complete"=>true,"skipped"=>true,"reason"=>"No common optimal aggregation size"),2);end
    exit(0)
end
large=300 in common ? 300 : maximum(common);sizes=sort!(unique([minimum(common),large]))
function family(bits)
 bits==0 && return 0
 masks=[m for m in 0:15 if bits&(1<<m)!=0];isempty(masks)&&error("Unreachable OD")
 m=argmin(m->(m&8!=0 ? 0 : m&2!=0 ? 1 : m&4!=0 ? 2 : 3,count_ones(m&14),m&14,m),masks)
 p=m&14;p==0 ? 1 : p==2 ? 2 : p==4 ? 3 : p==8 ? 4 : 5
end
const FAMILY_LOOKUP=[family(bits) for bits in 0:65535]
function count_design_modes(graph,lengths,bounds,nb,tag)
 n=graph.physical_nodes;workers=Threads.nthreads()
 histograms=[zeros(Int64,5,nb) for _ in 1:workers]
 completed=Threads.Atomic{Int}(0)
 # Keep task-local scratch and histogram bindings separate from the reduced
 # result. Reassigning a captured binding after @threads can share it across
 # tasks and lose counts even though the histogram arrays were separate.
 Threads.@threads :static for worker in 1:workers
  local workspace=MS.Workspace(graph)
  local histogram=histograms[worker]
  for origin in worker:workers:n
   MS.shortest_modes!(workspace,graph,origin)
   for destination in 1:n
    origin==destination&&continue
    band=searchsortedfirst(bounds,lengths[destination,origin])
    chosen=FAMILY_LOOKUP[Int(workspace.labels[destination])+1]
    chosen>0||error("Unreachable OD")
    histogram[chosen,band]+=1
   end
   progress=Threads.atomic_add!(completed,1)+1
   if progress%2000==0||progress==n;println("FULL REPLAY ",tag," ",progress,"/",n);flush(stdout);end
  end
 end
 totals=reduce(+,histograms)
 @assert sum(totals)==n*(n-1)
 totals
end
function subset(g,allowed)
 n=g.physical_nodes;arcs=Tuple{Int,Int,Int64}[]
 for u in 1:4n,e in g.offsets[u]:g.offsets[u+1]-1
  v=g.targets[e];transfer=(u<=n&&v==u+3n)||(v<=n&&u==v+3n)
  transfer&&!(min(u,v) in allowed)&&continue
  push!(arcs,(u,v,g.costs[e]))
 end
 MS.mode_graph(n,arcs)
end
rows=NamedTuple[]
for city in cities
 root=joinpath(B,"inputs/city_networks_walk_connected");g,_=MS.load_city(root,city);g=MS.with_car_costs(g,joinpath(root,"speed_overrides",city*".csv"));n=g.physical_nodes
 owner=only(filter(x->x["city"]==city,JSON.parsefile(joinpath(PROVENANCE,"osm_car_ownership_sources.json"))["cities"]))
 p=BigInt(owner["numerator"])//BigInt(owner["denominator"])
 profile=only(filter(x->x["city"]==city,JSON.parsefile(joinpath(PROVENANCE,"osm_distance_demand_sources.json"))["cities"]))
 bounds=Float64.(profile["upper_bounds_km"]);units=BigInt.(profile["mass_units"]);q=units.//sum(units)
 for k in sizes
  started=time();directory=joinpath(B,"raw/station_full_network_design_replay","$(city)_$(k)");mkpath(directory)
  if isfile(joinpath(directory,"result.json"))
   old=JSON.parsefile(joinpath(directory,"result.json"))
   push!(rows,(;city,aggregation_locations=k,opened_stations=old["opened_stations"],seconds=old["seconds"],(Symbol(f*"_pct")=>old[f*"_pct"] for f in ["walking","car","transit","bike_sharing","multi_powered"])...))
   CSV.write(joinpath(B,"data/station_full_network_design_replay.csv"),rows)
   println("CACHED FULL REPLAY ",city," ",k);continue
  end
  result=JSON.parsefile(joinpath(B,"raw/aggregation_cases","$(city)_$(k)_s20_b10/result.json"))
  maps=Dict(Int(r.aggregate_node)=>Int(r.original_node) for r in CSV.File(joinpath(B,"inputs/aggregation",city,"selected_nodes.csv")) if r.target==k)
  originals=sort!([maps[Int(i)] for i in result["open_candidate_ids"]]);@assert all(i->i in JSON.parsefile(joinpath(root,"graphs",city,"metadata.json"))["all_candidate_stations"],originals)
  design=subset(g,Set(originals));shares=zeros(Rational{BigInt},5);classrows=NamedTuple[]
  for (cl,car,weight) in [("car",true,p),("noncar",false,1-p)]
   fixed,_=FD.load_fixed_graph(joinpath(root,"graphs",city),[joinpath(B,"inputs/geometry/road_lengths",city*".csv"),joinpath(B,"inputs/geometry/transit_lengths",city*".csv")];car_available=car,car_override=joinpath(root,"speed_overrides",city*".csv"))
   cache=joinpath(car ? CAR_CACHE : NONCAR_CACHE,city,"paths_"*cl)
   if !isfile(joinpath(cache,"complete.json"))
    FD.cache_all_pairs(fixed,cache;progress=(d,n)->begin
     if d%1000==0||d==n;println("REBUILD REFERENCE ",city," ",cl," ",d,"/",n);flush(stdout);end
    end)
   end
   _,lengths=FD.load_cache(cache;expected_graph_hash=FD.fingerprint(fixed))
   graph=car ? design : MS.without_car(design);nb=length(q)
   h=count_design_modes(graph,lengths,bounds,nb,"$city k=$k $cl")
   counts=vec(sum(h,dims=1));@assert sum(counts)==n*(n-1)
   conditional=[sum(q[b]*(BigInt(h[f,b])//BigInt(counts[b])) for b in eachindex(q)) for f in 1:5];@assert sum(conditional)==1
   shares .+=weight.*conditional
   for f in 1:5;push!(classrows,(class=cl,mode=["walking","car","transit","bike_sharing","multi_powered"][f],conditional_pct=100Float64(conditional[f]),conditional_exact=string(conditional[f])));end
   lengths=nothing;fixed=nothing;GC.gc()
  end
  @assert sum(shares)==1
  row=(;city,aggregation_locations=k,opened_stations=length(originals),seconds=time()-started,(Symbol(f*"_pct")=>100Float64(shares[i]) for (i,f) in enumerate(["walking","car","transit","bike_sharing","multi_powered"]))...)
  push!(rows,row);CSV.write(joinpath(directory,"classes.csv"),classrows)
  open(joinpath(directory,"result.json"),"w") do io;JSON.print(io,merge(Dict(string(key)=>value for (key,value) in pairs(row)),Dict("opened_original_ids"=>originals,"exact_shares"=>string.(shares),"same_full_network_demand"=>true,"diagnostic_tie_rule"=>"bike then car; normal mixed-mode assignment","optimization"=>"Replay fixed optimized station choices; no reoptimization")),2);end
  println("FINISHED FULL REPLAY ",city," ",k);flush(stdout)
  CSV.write(joinpath(B,"data/station_full_network_design_replay.csv"),rows)
 end
end
open(joinpath(B,"raw/station_full_network_design_replay/COMPLETE.json"),"w") do io;JSON.print(io,Dict("complete"=>true,"common_sizes"=>sizes,"cities"=>cities),2);end
