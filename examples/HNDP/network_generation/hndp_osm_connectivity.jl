module OSMConnectivity
using Graphs

"""Retain the largest strongly connected walking component, with stable ties.

Walking arcs are normally bidirectional. SCCs also guarantee reachability if
future inputs contain directed pedestrian arcs. Ties retain the component
with the smallest original node ID. The returned node IDs are sorted.
"""
function walking_backbone_nodes(n::Int, walk_edges)
    n>0 || throw(ArgumentError("cannot clean an empty network"))
    graph=DiGraph(n)
    for (u,v) in walk_edges
        add_edge!(graph,u,v)
    end
    components=strongly_connected_components(graph)
    sort!(components,by=c->(-length(c),minimum(c)))
    selected=sort!(first(components))
    length(selected)>=2 || throw(ArgumentError("walking graph has no component with two reachable locations"))
    selected
end

function restrict_road_layers(coords,car,bike,walk)
    keep=walking_backbone_nodes(length(coords),keys(walk))
    index=Dict(v=>i for (i,v) in enumerate(keep))
    restrict(edges)=Dict((index[u],index[v])=>a for ((u,v),a) in edges if haskey(index,u)&&haskey(index,v))
    coords[keep],restrict(car),restrict(bike),restrict(walk),keep
end
end
