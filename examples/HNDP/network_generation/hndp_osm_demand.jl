module OSMDemand
using Random, JSON

"""Published distance-band masses; integer units avoid rounding the input shares."""
function profile(city; path=joinpath(@__DIR__,"..","sources","osm_distance_demand_sources.json"))
    records=JSON.parsefile(path)["cities"]
    matches=filter(r->r["city"]==city,records)
    length(matches)==1 || throw(ArgumentError("No unique demand profile for $city"))
    p=only(matches);validate(p);return p
end
function validate(p)
    edges=Float64.(p["upper_bounds_km"]);units=Int.(p["mass_units"])
    length(units)==length(edges)+1 || throw(ArgumentError("one mass required per distance band"))
    all(isfinite,edges)&&all(>(0),edges)&&all(>(0),diff(edges)) || throw(ArgumentError("bounds must increase and be positive"))
    all(>=(0),units)&&sum(units)>0 || throw(ArgumentError("invalid demand masses"))
    return nothing
end
@inline band(km,edges)=searchsortedfirst(edges,km)
@inline distance_km(a,b)=2*6371.0088*asin(min(1.0,sqrt(sind((b[1]-a[1])/2)^2+cosd(a[1])*cosd(b[1])*sind((b[2]-a[2])/2)^2)))

"""Finite equal-user demand construction, with largest-remainder band quotas.

Uniform rejection sampling within each band (with replacement). Eligibility is
explicit: the importer uses its existing directed-car-reachability condition.
No silent fallback when a requested band cannot be sampled. A detour multiplier
is an assumption, not an observed route length. Exact all-OD demand weighting is
provided by station_location_paper/scripts/prepare_instances.jl on full fixed-network paths.
"""
function sample_pairs(coords,count,rng,p; detour_factor=1.3,
    eligible=(o,d)->true,max_attempts=max(10000,10000count))
    validate(p)
    count>=0 || throw(ArgumentError("negative user count"))
    isfinite(detour_factor)&&detour_factor>0 || throw(ArgumentError("invalid detour factor"))
    length(coords)>=2 || throw(ArgumentError("two physical nodes required"))
    units=Int.(p["mass_units"]);total=sum(units)
    quotas=fld.(count .* units,total)
    order=sortperm(eachindex(units),by=i->(-mod(count*units[i],total),i))
    for i in order[1:count-sum(quotas)];quotas[i]+=1;end
    remaining=copy(quotas);pairs=Tuple{Int,Int}[];edges=Float64.(p["upper_bounds_km"])
    attempts=0
    while length(pairs)<count && attempts<max_attempts
        attempts+=1;o=rand(rng,eachindex(coords));d=rand(rng,eachindex(coords))
        o==d && continue
        b=band(detour_factor*distance_km(coords[o],coords[d]),edges)
        remaining[b]>0 || continue
        eligible(o,d) || continue
        push!(pairs,(o,d));remaining[b]-=1
    end
    all(iszero,remaining) || error("Cannot realize distance demand; unfilled band quotas $remaining after $attempts attempts. Check geographic support/eligibility; no fallback was used.")
    shuffle!(rng,pairs)
    return pairs,Dict("profile"=>p,"distance_proxy"=>"detour_factor * haversine distance",
        "detour_factor"=>detour_factor,"band_user_counts"=>quotas,
        "sampling"=>"largest-remainder quotas; uniform eligible ordered pairs within each band, with replacement",
        "eligibility"=>"caller predicate","attempts"=>attempts)
end
end
