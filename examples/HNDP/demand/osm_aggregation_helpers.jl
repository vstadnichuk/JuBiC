# Reusable aggregation primitives for the current OSM pipeline.
# Moved without algorithm changes from the retired benchmark drivers.
# Selection, closures and redundancy pruning share the same implementation.
module OSMAggregation
using Graphs, CSV

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
    Threads.@threads for si in eachindex(sources)
        source=sources[si]
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


"""Shortest complete journeys on every supplied layer/transfer, then cheaper-direction symmetry.
The caller supplies the full original graph with all EXISTING sharing sites open.
"""
function multimodal_site_metric(graph, costs, sites)
    directed=station_distances_to_targets(graph,costs,sites,sites,(u,v)->true)
    metric=min.(directed,permutedims(directed))
    all(x->0<=x<typemax(Int),metric) || throw(ArgumentError("Selection sites must be reachable on the full multimodal graph."))
    metric
end

"""Network-medoid first, then farthest-point coverage. Return an ordered nested sequence."""
function potential_station_sequence(metric, count; eligible=collect(axes(metric,1)))
    pool=sort!(unique(collect(eligible)))
    0<=count<=length(pool) || throw(ArgumentError("Too few existing eligible sites."))
    size(metric,1)==size(metric,2) || throw(ArgumentError("Square selection metric required."))
    count==0 && return Int[]
    all(i->1<=i<=size(metric,1),pool) || throw(ArgumentError("Invalid eligible ID."))
    all(0<=metric[i,j]<typemax(Int) for i in pool,j in pool) || throw(ArgumentError("Unreachable candidate reference."))
    chosen=[argmin(i->(sum(Int128(metric[i,j]) for j in pool),i),pool)]
    while length(chosen)<count
        push!(chosen,argmax(i->(i in chosen ? Int128(-1) : Int128(minimum(metric[i,j] for j in chosen)),-i),pool))
    end
    chosen
end

function save_selection_metric(folder,k,metric)
    size(metric)==(k,k) || throw(ArgumentError("Selection metric size mismatch."))
    CSV.write(joinpath(folder,"selection_metric_$(k).csv"),
        ((source=i,target=j,cost_cents=metric[i,j]) for j in 1:k for i in 1:k))
end

function load_selection_metric(folder,k)
    file=joinpath(folder,"selection_metric_$(k).csv")
    isfile(file) || throw(ArgumentError("Regenerate aggregation with full-network selection costs: $file"))
    metric=fill(typemax(Int),k,k);seen=falses(k,k)
    for r in CSV.File(file)
        i,j=Int(r.source),Int(r.target)
        1<=i<=k && 1<=j<=k && !seen[i,j] || throw(ArgumentError("Invalid/repeated selection OD."))
        metric[i,j]=Int(r.cost_cents);seen[i,j]=true
    end
    all(seen) && metric==permutedims(metric) && all(i->metric[i,i]==0,1:k) || throw(ArgumentError("Incomplete/asymmetric selection metric."))
    metric
end

function balanced_medoids(metric, target)
    target>0 || throw(ArgumentError("Positive aggregation target required."))
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
        for c in sort!(collect(1:k), by=c->(length(clusters[c]),c))
            length(clusters[c]) >= capacity[c] && continue
            isempty(unassigned) && break
            x=argmin(y->(metric[y,seeds[c]] == inf ? typemax(Int) : metric[y,seeds[c]],y), unassigned)
            push!(clusters[c],x); delete!(unassigned,x)
        end
    end
    medoids=Int[]; sizes=Int[]
    for cluster in clusters
        isempty(cluster) && continue
        medoid=argmin(x->(sum(metric[x,y] == inf ? 1.0e18 : Float64(metric[x,y]) for y in cluster),x),cluster)
        push!(medoids,medoid); push!(sizes,length(cluster))
    end
    medoids, sizes
end


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
    Threads.@threads for si in eachindex(sources)
        source=sources[si]
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
    Threads.@threads for si in eachindex(sources)
        source=sources[si]
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

end # module OSMAggregation
