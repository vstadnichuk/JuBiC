# Routing histogram helpers retained for regression tests; demand is built by station_location_paper/scripts/prepare_instances.jl.
if !isdefined(@__MODULE__, :OSMModeShares)
    include("osm_mode_share_analysis.jl")
end
if !isdefined(@__MODULE__, :OSMDemand)
    include("../network_generation/hndp_osm_demand.jl")
end
using CSV, JSON, Base.Threads

const DEMAND_METRICS=vcat(OSMModeShares.FAMILY_NAMES,["unreachable","any_multimodal","any_multiple_powered_modes"])

function demand_label_metrics()
    table=zeros(Bool,9,65536);table[7,1]=true
    for bits in 1:65535
        families=Set{Int}()
        for mask in 0:15
            bits&(1<<mask)==0 && continue
            powered=mask>>1
            push!(families,powered==0 ? 1 : powered==1 ? 2 : powered==2 ? 3 : powered==4 ? 4 : 5)
            table[8,bits+1] |= count_ones(mask)>=2
            table[9,bits+1] |= count_ones(powered)>=2
        end
        table[length(families)==1 ? only(families) : 6,bits+1]=true
    end
    table
end

"""Two Dijkstras per origin, reused for every demand and distance scenario.

Only distance-band counters are retained; no OD matrix or pairwise routing.
Both availability copies retain every non-self OD, including unreachable ones.
"""
function count_distance_demand(owner,nonowner,coords,p;factors=[1.0,1.3,1.5])
    OSMDemand.validate(p)
    all(f->isfinite(f)&&f>0,factors) || error("invalid distance factors")
    n=owner.physical_nodes;@assert nonowner.physical_nodes==length(coords)==n
    edges=Float64.(p["upper_bounds_km"]);nb=length(edges)+1
    workers=min(Threads.nthreads(),n);lookup=demand_label_metrics()
    counters=[zeros(Int64,9,nb,length(factors),2) for _ in 1:workers]
    totals=[zeros(Int64,nb,length(factors)) for _ in 1:workers]
    distance_sums=[zeros(Float64,nb,length(factors)) for _ in 1:workers]
    done=Threads.Atomic{Int}(0);started=time()
    Threads.@threads :static for worker in 1:workers
        a=OSMModeShares.Workspace(owner);b=OSMModeShares.Workspace(nonowner)
        h=counters[worker];t=totals[worker];ds=distance_sums[worker]
        for o in worker:workers:n
            OSMModeShares.shortest_modes!(a,owner,o)
            OSMModeShares.shortest_modes!(b,nonowner,o)
            @inbounds for d in 1:n
                o==d && continue
                km=OSMDemand.distance_km(coords[o],coords[d])
                labels=(Int(a.labels[d])+1,Int(b.labels[d])+1)
                for f in eachindex(factors)
                    dist=km*factors[f];band=OSMDemand.band(dist,edges)
                    t[band,f]+=1;ds[band,f]+=dist
                    for c in 1:2,k in 1:9
                        h[k,band,f,c]+=lookup[k,labels[c]]
                    end
                end
            end
            completed=Threads.atomic_add!(done,1)+1
            if completed%2000==0 || completed==n
                println("  ",completed,"/",n," origins; ",round(time()-started,digits=1)," seconds");flush(stdout)
            end
        end
    end
    return reduce(+,counters),reduce(+,totals),reduce(+,distance_sums)
end

