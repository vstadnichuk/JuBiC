# Full-OD OSM benchmark instances, independent of solver selection.
# Historical benchmark outputs retain their original two-sided cost convention.
using CSV, JSON, Graphs, SparseArrays, SHA
if !isdefined(@__MODULE__, :HNDPGeneratedNetwork)
    include("../network_generation/hndp_osm_network_io.jl")
end
if !isdefined(@__MODULE__, :OSMAggregation)
    include("../demand/osm_aggregation_helpers.jl")
end

"""
    build_osm_full_od_benchmark(city, k; aggregation_root, network_root, demand_root,
                                station_construction_cost=0.0, candidate_count=20, station_budget=10)

Reuse saved aggregate topology/costs, select the requested sharing candidates by full-network cost, allow the budgeted subset
open, and load full-network demand already aggregated separately by class.
Diagonal mass is retained in metadata, without creating zero-length followers.
Charge K once on station entry; return remains controllable with zero cost.
"""
function build_osm_full_od_benchmark(city, k;
    aggregation_root::AbstractString, network_root::AbstractString, demand_root::AbstractString,
    station_construction_cost=0.0, candidate_count=20, station_budget=10)
    repo_root = normpath(joinpath(@__DIR__, "..", "..", ".."))
    directory=joinpath(aggregation_root,city)
    sites=filter(r->r.target==k,collect(CSV.File(joinpath(directory,"selected_nodes.csv"))))
    sort!(sites,by=r->r.aggregate_node)
    @assert 0<=station_budget<=candidate_count<=k
    metadata=JSON.parsefile(joinpath(network_root,"graphs",city,"metadata.json"))
    allowed=Set(Int.(metadata["all_candidate_stations"]))
    eligible=[i for i in 1:k if Int(sites[i].original_node) in allowed]
    metric=OSMAggregation.load_selection_metric(directory,k)
    chosen=sort!(OSMAggregation.potential_station_sequence(metric,candidate_count;eligible))
    station_set=Set(chosen)
    graph=DiGraph(4k); costs=spzeros(Float64,4k,4k);risk=spzeros(Float64,4k,4k)
    for row in CSV.File(joinpath(directory,"graph_$(k)_edges.csv"))
        u=Int(row.source);v=Int(row.target)
        bike_transfer=(u<=k && v==u+3k)||(v<=k && u==v+3k)
        bike_transfer && !(min(u,v) in station_set) && continue
        add_edge!(graph,u,v);costs[u,v]=row.cost_cents
    end
    metadata=JSON.parsefile(joinpath(network_root,"graphs",city,"metadata.json"))
    fee=Float64(metadata["bike_entry_cost_cents"])
    groups=[[(i,i+3k),(i+3k,i)] for i in chosen]
    decisions, prices = hndp_sharing_station_decisions(groups; construction_cost=station_construction_cost)
    for i in chosen
        @assert has_edge(graph,i,i+3k) && has_edge(graph,i+3k,i)
        risk[i,i+3k]=-fee
    end
    noncar=copy(graph)
    for edge in collect(edges(noncar))
        u=src(edge);v=dst(edge)
        (k<u<=2k || k<v<=2k) && rem_edge!(noncar,u,v)
    end
    profiles=map((graph,noncar)) do allowed
        fixed=copy(allowed)
        for arc in decisions;rem_edge!(fixed,arc...);end
        HNDPProfileCosts(costs,allowed,Float64.(floyd_warshall_shortest_paths(fixed,costs).dists))
    end
    own=only(filter(c->c["city"]==city,JSON.parsefile(joinpath(repo_root,"examples/HNDP/sources/osm_car_ownership_sources.json"))["cities"]))
    p=BigInt(own["numerator"])//BigInt(own["denominator"])
    demand_file=joinpath(demand_root,city,"demand_$(k).csv")
    isfile(demand_file) || throw(ArgumentError("Generate full-network fixed-path demand with prepare_osm_fixed_demand before building this instance: $demand_file"))
    owner=zeros(Rational{BigInt},k,k);nonowner=zeros(Rational{BigInt},k,k)
    seen=Set{Tuple{Int,Int}}()
    for row in CSV.File(demand_file;types=String)
        i=parse(Int,row.origin);j=parse(Int,row.destination)
        1<=i<=k && 1<=j<=k && !((i,j) in seen) || throw(ArgumentError("Invalid or repeated aggregated OD."))
        push!(seen,(i,j))
        owner[i,j]=parse(BigInt,row.car_numerator)//parse(BigInt,row.car_denominator)
        nonowner[i,j]=parse(BigInt,row.noncar_numerator)//parse(BigInt,row.noncar_denominator)
        owner[i,j]>=0 && nonowner[i,j]>=0 || throw(ArgumentError("Negative aggregated demand."))
    end
    @assert length(seen)==k*k && sum(owner)==p && sum(nonowner)==1-p
    diagonal_owner=sum(owner[i,i] for i in 1:k);diagonal_nonowner=sum(nonowner[i,i] for i in 1:k)
    users=User[];mass=zero(p);owner_mass=zero(p)
    for i=1:k,j=1:k
        i==j && continue
        for (profile,demand) in enumerate((owner[i,j],nonowner[i,j]))
            mass+=demand;profile==1 && (owner_mass+=demand)
            push!(users,User("$(i)_$(j)_$(profile)",i,j,HNDPWeightedRisk(risk,Float64(demand)),profiles[profile],nothing,nothing))
        end
    end
    @assert mass+diagonal_owner+diagonal_nonowner==1 && owner_mass+diagonal_owner==p && length(users)==2k*(k-1)
    @assert all(!(k<src(e)<=2k || k<dst(e)<=2k) for e in edges(noncar))
    @assert all(isfinite(profiles[profile].fixed_distances[i,j]) for profile=1:2,i=1:k,j=1:k)
    name="$(city)_n$(k)_s$(candidate_count)_b$(station_budget)_fullod_fixed_demand"
    meta=Dict{String,Any}("name"=>name,"nnodes"=>4k,"narcs"=>ne(graph),"nusers"=>length(users),"aggregation_target_nodes"=>k,"bike_possible_stations"=>candidate_count,"decision_budget"=>station_budget,"candidate_ids"=>chosen,"demand_total"=>"1","owner_mass"=>string(p),"station_construction_costs"=>[prices[first(group)] for group in groups],"station_construction_cost_convention"=>"entry_only","bike_entry_revenue_cents"=>fee,
        "demand_method"=>"class_specific_full_fixed_network_before_aggregation","modeled_inter_site_mass"=>string(mass),"modeled_owner_mass"=>string(owner_mass),
        "intrazonal_owner_mass"=>string(diagonal_owner),"intrazonal_nonowner_mass"=>string(diagonal_nonowner),"demand_sha256"=>bytes2hex(sha256(read(demand_file))))
    meta["candidate_original_ids"]=[Int(sites[i].original_node) for i in chosen]
    meta["eligible_sharing_site_count"]=length(eligible)
    meta["selection_metric"]="full_multimodal_all_existing_sharing_open"
    instance=HNDPwC(graph,users,decisions,prices,nothing,groups,station_budget)
    HNDPGeneratedNetwork(name,instance,meta)
end
