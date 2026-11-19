# Shared closure and histogram functions. Legacy standalone aggregation runs are archived.
# Exhaustive aggregation sensitivity; uses the existing balanced medoid selector
# and existing shortest-path closure / redundancy-pruning implementations.
include("osm_distance_demand_analysis.jl")
using Graphs, SparseArrays
if !isdefined(@__MODULE__, :OSMConnectivity)
    include("../network_generation/hndp_osm_connectivity.jl")
end
if !isdefined(@__MODULE__, :OSMAggregation)
    include("osm_aggregation_helpers.jl")
end
# Preserve the namespace used by archived callers while maintaining helpers here.
const ExistingAggregation = OSMAggregation
const AGG_TARGETS=[50,75,100,200,300,400]
const AGG_COST_CONFIGURATION="car_eur0p22_speed0p50_fallback20_vot13"
const AGG_BIKE_SETTING="sharing"
const AGG_PRIMARY_WEIGHTING="distance_weighted"

function sharing_marker_compatible(marker)
    get(marker,"bike_setting",nothing)==AGG_BIKE_SETTING &&
        get(marker,"primary_weighting",nothing)==AGG_PRIMARY_WEIGHTING &&
        get(marker,"analysis_schema",nothing)==5 &&
        get(marker,"cost_configuration",nothing)==AGG_COST_CONFIGURATION &&
        get(marker,"selection_metric",nothing)=="full_multimodal_all_existing_sharing_open" &&
        get(marker,"targets",nothing)==AGG_TARGETS &&
        get(marker,"walking_connectivity_policy",nothing)=="largest_walking_scc"
end

function sharing_cache_complete(folder)
    path=joinpath(folder,"COMPLETE.json")
    isfile(path) || return false
    sharing_marker_compatible(JSON.parsefile(path)) ||
        throw(ArgumentError("$folder contains a legacy/different analysis; use a fresh output directory for sharing-only distance-weighted results"))
    true
end

function split_layers(g)
    n=g.physical_nodes
    graphs=[DiGraph(n) for _ in 1:4];costs=[spzeros(Int,n,n) for _ in 1:4]
    for u in 1:4n,e in g.offsets[u]:g.offsets[u+1]-1
        v=g.targets[e];l=cld(u,n)
        l==cld(v,n) || continue
        a=mod1(u,n);b=mod1(v,n);add_edge!(graphs[l],a,b);costs[l][a,b]=g.costs[e]
    end
    graphs,costs
end

function aggregate_existing(g,selected,graphs,costs)
    n=g.physical_nodes;k=length(selected);arcs=Tuple{Int,Int,Int64}[]
    for l in 1:4
        closure,_=ExistingAggregation.mode_closure(graphs[l],costs[l],spzeros(Int,n,n),selected,(u,v)->true)
        ranks=ExistingAggregation.path_redundancy_ranks(graphs[l],costs[l],selected,(u,v)->true)
        for i in 1:k,j in 1:k
            i==j && continue
            closure[i,j]==typemax(Int) && continue
            ranks[i,j]<=k && continue
            push!(arcs,((l-1)*k+i,(l-1)*k+j,closure[i,j]))
        end
    end
    index=Dict(s=>i for (i,s) in enumerate(selected))
    for s in selected,l in 1:4
        u=s+(l-1)*n
        for e in g.offsets[u]:g.offsets[u+1]-1
            v=g.targets[e];lv=cld(v,n)
            lv==l && continue
            @assert mod1(v,n)==s
            push!(arcs,((l-1)*k+index[s],(lv-1)*k+index[s],g.costs[e]))
        end
    end
    result=OSMModeShares.mode_graph(k,arcs)
    # Redundancy pruning must preserve every pure-mode shortest distance.
    ag,ac=split_layers(result)
    for l in 1:4
        before=ExistingAggregation.station_distances_to_targets(graphs[l],costs[l],selected,selected,(u,v)->true)
        after=ExistingAggregation.station_distances_to_targets(ag[l],ac[l],collect(1:k),collect(1:k),(u,v)->true)
        @assert before==after
    end
    result
end

# Cache labels/costs for selected endpoints. When exhaustive=true, count all
# physical OD pairs too. Histograms preserve all exact optimum-mode ties.
function evaluate_graph(g,coords,profile,selected;exhaustive=false)
    n=g.physical_nodes;m=length(selected);nb=length(profile["upper_bounds_km"])+1
    workers=min(Threads.nthreads(),n);nonowner=OSMModeShares.without_car(g)
    labels=zeros(UInt16,m,m,2);dist=fill(typemax(Int64),m,m,2)
    index=zeros(Int,n);index[selected]=1:m
    hist=[zeros(Int64,65536,nb,2) for _ in 1:workers]
    totals=[zeros(Int64,nb) for _ in 1:workers]
    origins=exhaustive ? collect(1:n) : selected
    done=Threads.Atomic{Int}(0);started=time()
    Threads.@threads :static for worker in 1:workers
        a=OSMModeShares.Workspace(g);b=OSMModeShares.Workspace(nonowner)
        for oi in worker:workers:length(origins)
            o=origins[oi];OSMModeShares.shortest_modes!(a,g,o);OSMModeShares.shortest_modes!(b,nonowner,o)
            si=index[o]
            if si>0
                for (j,d) in enumerate(selected)
                    labels[si,j,1]=a.labels[d];labels[si,j,2]=b.labels[d]
                    dist[si,j,1]=a.distance[d];dist[si,j,2]=b.distance[d]
                end
            end
            if exhaustive
                @inbounds for d in 1:n
                    o==d && continue
                    band=OSMDemand.band(1.3*OSMDemand.distance_km(coords[o],coords[d]),profile["upper_bounds_km"])
                    hist[worker][Int(a.labels[d])+1,band,1]+=1
                    hist[worker][Int(b.labels[d])+1,band,2]+=1
                    totals[worker][band]+=1
                end
            end
            q=Threads.atomic_add!(done,1)+1
            if q%2000==0 || q==length(origins)
                println("    ",q,"/",length(origins)," origins ",round(time()-started,digits=1)," s");flush(stdout)
            end
        end
    end
    labels,dist,reduce(+,hist),reduce(+,totals)
end

function selected_hist(labels,coords,profile)
    k=length(coords);nb=length(profile["upper_bounds_km"])+1
    h=zeros(Int64,65536,nb,2);t=zeros(Int64,nb)
    for i in 1:k,j in 1:k
        i==j && continue
        band=OSMDemand.band(1.3*OSMDemand.distance_km(coords[i],coords[j]),profile["upper_bounds_km"])
        t[band]+=1
        for c in 1:2;h[Int(labels[i,j,c])+1,band,c]+=1;end
    end
    h,t
end
const AGG_LOOKUP=demand_label_metrics()
function share_rows(h,t,profile,p,city,bike,population,k)
    @assert sum(h[:,:,1])==sum(h[:,:,2])==sum(t)==k*(k-1)
    @assert all(iszero(h[bits+1,b,2]) for bits in 0:65535 for b in eachindex(t) if bits&0xcccc!=0)
    units=BigInt.(profile["mass_units"]);mass=units.//sum(units)
    missing=findall(b->t[b]==0 && mass[b]>0,eachindex(t))
    metrics=zeros(Int64,9,length(t),2)
    for b in eachindex(t),c in 1:2,bits in 1:65536
        z=h[bits,b,c];z==0 && continue
        for metric in 1:9;metrics[metric,b,c]+=z*AGG_LOOKUP[metric,bits];end
    end
    rows=NamedTuple[]
    for weighting in ("uniform_od","distance_weighted")
        valid=weighting=="uniform_od" || isempty(missing)
        exact=zeros(Rational{BigInt},9)
        if valid
            for b in eachindex(t),metric in 1:9
                t[b]==0 && continue
                q=weighting=="uniform_od" ? BigInt(t[b])//BigInt(sum(t)) : mass[b]
                exact[metric]+=q*(p*metrics[metric,b,1]+(1-p)*metrics[metric,b,2])/t[b]
            end
            @assert sum(exact[1:7])==1
        end
        for metric in 1:9
            push!(rows,(;city,bike,population,physical_nodes=k,layered_nodes=4k,ordered_pairs=k*(k-1),weighting,
                metric=DEMAND_METRICS[metric],supported=valid,share=valid ? Float64(exact[metric]) : missing_value(),
                numerator=valid ? string(numerator(exact[metric])) : "",denominator=valid ? string(denominator(exact[metric])) : "",
                unsupported_bands=join(missing,";"),unsupported_target_mass=Float64(sum(mass[missing];init=0//BigInt(1)))))
        end
    end
    rows
end
missing_value()=NaN

