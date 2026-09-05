using CSV
using DataFrames
using Graphs
using JuBiC

include(joinpath(@__DIR__, "aachen_osm_gtfs_generation.jl"))

# Kept here rather than in the retired station-pair diagnostic script because
# the balanced station selector needs these distances as a core dependency.
function _balanced_heap_push!(heap, item)
    push!(heap, item)
    i = length(heap)
    while i > 1
        p = fld(i, 2)
        heap[p][1] <= heap[i][1] && break
        heap[p], heap[i] = heap[i], heap[p]
        i = p
    end
end

function _balanced_heap_pop!(heap)
    item = heap[1]
    last_item = pop!(heap)
    isempty(heap) && return item
    heap[1] = last_item
    i = 1
    while true
        left, right, smallest = 2i, 2i + 1, i
        left <= length(heap) && heap[left][1] < heap[smallest][1] && (smallest = left)
        right <= length(heap) && heap[right][1] < heap[smallest][1] && (smallest = right)
        smallest == i && break
        heap[i], heap[smallest] = heap[smallest], heap[i]
        i = smallest
    end
    item
end

function station_distances_to_targets(graph, costs, sources, targets, allowed)
    inf = typemax(Int)
    result = fill(inf, length(sources), length(targets))
    target_index = Dict(target => i for (i, target) in enumerate(targets))
    for (si, source) in enumerate(sources)
        distances = fill(inf, nv(graph))
        distances[source] = 0
        heap = [(0, source)]
        while !isempty(heap)
            current, node = _balanced_heap_pop!(heap)
            current == distances[node] || continue
            for next_node in outneighbors(graph, node)
                allowed(node, next_node) || continue
                weight = costs[node, next_node]
                weight == inf && continue
                candidate = current + weight
                candidate < distances[next_node] || continue
                distances[next_node] = candidate
                _balanced_heap_push!(heap, (candidate, next_node))
            end
        end
        for (target, ti) in target_index
            result[si, ti] = distances[target]
        end
    end
    result
end

const CITIES = [
    (name="cologne", osm="data/osm/comparative_cities/cologne_buffered.json", boundary="data/osm/comparative_cities/cologne_boundary.json", gtfs="data/osm/comparative_cities/cologne_gtfs"),
    (name="heidelberg", osm="data/osm/comparative_cities/heidelberg.json", boundary="data/osm/comparative_cities/heidelberg_boundary.json", gtfs="data/osm/comparative_cities/heidelberg_gtfs"),
    (name="karlsruhe", osm="data/osm/comparative_cities/karlsruhe_buffered.json", boundary="data/osm/comparative_cities/karlsruhe_boundary.json", gtfs="data/osm/comparative_cities/karlsruhe_gtfs"),
]
const TARGETS=[50,100,150,200]
const OUT="tmp_compare/reduced_station_balanced_spatial_sweep.csv"

function balanced_medoids(metric, target)
    n=size(metric,1); k=min(target,n); inf=typemax(Int)
    seeds=Int[1]
    while length(seeds)<k
        best=0; bestsep=-1.0
        for candidate in 1:n
            candidate in seeds && continue
            sep=minimum(metric[candidate,s] == inf ? 0.0 : Float64(metric[candidate,s]) for s in seeds)
            sep>bestsep && (bestsep=sep; best=candidate)
        end
        best==0 && break
        push!(seeds,best)
    end
    k=length(seeds); capacity=fill(div(n,k),k); capacity[1:rem(n,k)] .+= 1
    clusters=[Int[] for _ in 1:k]
    for (cluster,seed) in enumerate(seeds); push!(clusters[cluster],seed); end
    remaining=[x for x in 1:n if !(x in seeds)]
    sort!(remaining, by=x->minimum(metric[x,s] == inf ? typemax(Int) : metric[x,s] for s in seeds), rev=true)
    # Fill every cluster to its exact target size in rounds. This guarantees
    # cluster sizes differing by at most one while still assigning each point
    # to a nearby seed whenever possible.
    unassigned=Set(remaining)
    while !isempty(unassigned)
        for c in sortperm(capacity, by=c->length(clusters[c]))
            length(clusters[c]) >= capacity[c] && continue
            isempty(unassigned) && break
            x=argmin(y->(metric[y,seeds[c]] == inf ? typemax(Int) : metric[y,seeds[c]]), unassigned)
            push!(clusters[c],x); delete!(unassigned,x)
        end
    end
    medoids=Int[]; sizes=Int[]
    for cluster in clusters
        isempty(cluster) && continue
        medoid=argmin(x->sum(metric[x,y] == inf ? 1.0e18 : Float64(metric[x,y]) for y in cluster),cluster)
        push!(medoids,medoid); push!(sizes,length(cluster))
    end
    medoids, sizes
end

function geo_metric(coords)
    n=length(coords); result=zeros(Int,n,n)
    for i in 1:n, j in 1:n
        result[i,j]=round(Int,1_000_000*sqrt((coords[i][1]-coords[j][1])^2 + ((coords[i][2]-coords[j][2])*cosd((coords[i][1]+coords[j][1])/2))^2))
    end
    result
end

function stats(car,transit,bike,selected)
    inf=typemax(Int); total=0; c=0; t=0; b=0; ties=0
    cr=0; tr=0; br=0
    for i in selected,j in selected
        i==j && continue; total+=1; v=(car[i,j],transit[i,j],bike[i,j])
        cr+=v[1]!=inf; tr+=v[2]!=inf; br+=v[3]!=inf
        finite=[x for x in v if x!=inf]; isempty(finite) && continue
        w=findall(x->x==minimum(finite),v)
        length(w)==1 ? ((:car,:transit,:bike)[first(w)]==:car ? (c+=1) : (:car,:transit,:bike)[first(w)]==:transit ? (t+=1) : (b+=1)) : (ties+=1)
    end
    total,(c,t,b,ties),(cr,tr,br)
end

function main()
    rows=NamedTuple[]
    for city in CITIES
        println("Generating ",city.name)
        g=generate_aachen_osm_gtfs(osm_path=city.osm,boundary_path=city.boundary,use_boundary_compression=true,
            gtfs_dir=city.gtfs,nusers=1,seed=20260829,bike_station_count=10000,all_nodes_bike_stations=false,
            generalized_cost=true,bike_entry_cost_cents=0,bike_cost_per_km_cents=0,
            bike_daily_subscription_cost_cents=nothing,transit_monthly_cost_cents=6300,transit_rides_per_month=60)
        h=g.instance; n=g.metadata["nnodes"]÷4; stations=g.metadata["stations"]; costs=h.users[1].mcost; coords=g.metadata["coords"]
        car_allowed=(u,v)->(u<=n && n<v<=2n)||(n<u<=2n && v<=n)||(n<u<=2n && n<v<=2n)
        bike_allowed=(u,v)->(u<=n && 3n<v<=4n)||(3n<u<=4n && v<=n)||(3n<u<=4n && 3n<v<=4n)
        transit_allowed=(u,v)->(u<=n && v<=n)||(2n<u<=3n && 2n<v<=3n)||(u<=n && 2n<v<=3n)||(2n<u<=3n && v<=n)
        car=station_distances_to_targets(h.mygraph,costs,stations,stations,car_allowed)
        transit=station_distances_to_targets(h.mygraph,costs,stations,stations,transit_allowed)
        bike=station_distances_to_targets(h.mygraph,costs,stations,stations,bike_allowed)
        network_metric=copy(car); geographic_metric=geo_metric([coords[s] for s in stations])
        for i in axes(network_metric,1),j in axes(network_metric,2)
            network_metric[i,j]=network_metric[i,j]==typemax(Int) ? typemax(Int) : min(network_metric[i,j],network_metric[j,i])
        end
        for target in TARGETS, (strategy,metric) in (("network_balanced",network_metric),("geographic_balanced",geographic_metric))
            selected,cluster_sizes=balanced_medoids(metric,target); total,(c,t,b,ti),(cr,tr,br)=stats(car,transit,bike,selected); k=length(selected)
            arcs=sum(count(i!=j && matrix[selected[i],selected[j]]!=typemax(Int) for i in 1:k,j in 1:k) for matrix in (car,transit,bike))
            push!(rows,(city=city.name,strategy=strategy,target_stations=target,retained_stations=k,original_stations=length(stations),ordered_pairs=total,
                aggregate_nodes=3k,aggregate_modal_arcs=arcs,cluster_minimum=minimum(cluster_sizes),cluster_maximum=maximum(cluster_sizes),
                car_best=c,transit_walk_best=t,bike_best=b,ties=ti,car_best_share=c/total,transit_walk_best_share=t/total,bike_best_share=b/total,tie_share=ti/total,
                car_reachable=cr,transit_walk_reachable=tr,bike_reachable=br,cost_setting="generalized cost with explicit entry/exit; bike entry/km fees zero; Deutschlandticket=6300 cents/month over 60 rides"))
            CSV.write(OUT,DataFrame(rows)); println(city.name," ",strategy," k=",target," car=",c," transit=",t," bike=",b)
        end
    end
end
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
