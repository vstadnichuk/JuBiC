using CSV, JSON
if !isdefined(@__MODULE__, :_read_osm_ways)
    include("hndp_osm_network_generation.jl")
end
"""Recover exact OSM lengths for the retained saved road arcs, without resnapping.

The source bundle contains original road geometry and municipal boundaries.
The target bundle supplies the saved cleaned node coordinates and arc set.
"""
function recover_osm_road_lengths(city;raw_root,network_root,output_root)
    folder=joinpath(raw_root,"sources",city)
    meta=JSON.parsefile(joinpath(network_root,"graphs",city,"metadata.json"));n=Int(meta["nnodes"])÷4
    coords,car,bike,walk=_read_osm_ways(joinpath(folder,"roads_buffered.json"))
    coords,car,bike,walk=_compress_external_car_edges(coords,car,bike,walk,joinpath(folder,"boundary.json"))
    ids=Dict((Float64(c[1]),Float64(c[2]))=>i for (i,c) in enumerate(meta["coords"]))
    lengths=Dict{Tuple{Int,Int},Tuple{Float64,String}}()
    for (layer,roads) in [(0,walk),(n,car),(3n,bike)]
        for ((u,v),attr) in roads
            a=get(ids,coords[u],0);b=get(ids,coords[v],0)
            a>0 && b>0 || continue
            lengths[(a+layer,b+layer)]=(Float64(attr.distance),attr.kind=="external_shortcut" ? "osm_external_shortcut" : "osm_geometry")
        end
    end
    rows=NamedTuple[]
    for edge in CSV.File(joinpath(network_root,"graphs",city,"edges.csv"))
        u=Int(edge.source);v=Int(edge.target)
        cld(u,n)==cld(v,n) || continue
        2n<u<=3n && continue
        haskey(lengths,(u,v)) || throw(ArgumentError("Missing road geometry for $city arc $u → $v."))
        length,source=lengths[(u,v)]
        push!(rows,(source=u,target=v,length_km=length,length_source=source))
    end
    mkpath(output_root);path=joinpath(output_root,city*".csv");CSV.write(path,rows)
    path
end
