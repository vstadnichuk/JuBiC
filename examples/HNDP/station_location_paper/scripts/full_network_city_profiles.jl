#= All full-network ordered nonself ODs; all EXISTING sharing candidates open.

Usage: julia --project=<CSV/JSON environment> --threads=8 full_network_city_profiles.jl [cache-root] [city ...]
Missing fixed-network caches are rebuilt from bundled geometry/graphs. Caches
are optional acceleration, not extra demand inputs. No optimization is solved.
=#
using CSV, JSON, SHA, Mmap, Dates, Base.Threads
include(joinpath(@__DIR__, "paths.jl"))
include(joinpath(SRC, "demand", "osm_mode_share_analysis.jl"))
include(joinpath(SRC, "demand", "osm_fixed_network_demand.jl"))
const MS=OSMModeShares
const FD=OSMFixedNetworkDemand
const NAMES=["walking","car","transit","bike_sharing","multi_powered","tie","unreachable"]

function category(bits::Int)
    bits==0 && return 7
    families=0
    for mask in 0:15
        bits&(1<<mask)==0 && continue
        p=mask>>1
        f=p==0 ? 1 : p==1 ? 2 : p==2 ? 3 : p==4 ? 4 : 5
        families |= 1<<(f-1)
    end
    return count_ones(families)==1 ? trailing_zeros(families)+1 : 6
end

function smoke()
    # Walking access/egress does not change the single non-walking family.
    @assert category(1<<1)==1
    @assert category(1<<2)==category(1<<3)==2
    @assert category((1<<4)|(1<<5))==3
    @assert category(1<<9)==4
    @assert category(1<<7)==category(1<<13)==5
    @assert category((1<<3)|(1<<5))==6
    # A car/transit optimum tie propagated through zero-cost exits.
    g=MS.mode_graph(3,[(1,2,10),(1,4,2),(4,5,2),(5,2,1),
        (1,7,2),(7,8,3),(8,2,0)])
    w=MS.Workspace(g);MS.shortest_modes!(w,g,1)
    @assert w.distance[2]==5 && category(Int(w.labels[2]))==6
end

function count_hist(g,lengths,bounds,progress)
    n=g.physical_nodes;nb=length(bounds)+1;workers=min(n,Threads.nthreads())
    counters=[zeros(Int64,65536,nb) for _ in 1:workers]
    done=Threads.Atomic{Int}(0);started=time()
    Threads.@threads :static for worker in 1:workers
        w=MS.Workspace(g);h=counters[worker]
        for origin in worker:workers:n
            MS.shortest_modes!(w,g,origin)
            @inbounds for dest in 1:n
                origin==dest && continue
                km=lengths[dest,origin]
                isfinite(km) || error("Unreachable fixed reference $origin->$dest")
                band=searchsortedfirst(bounds,km)
                h[Int(w.labels[dest])+1,band]+=1
            end
            finished=Threads.atomic_add!(done,1)+1
            if finished%1000==0 || finished==n
                progress(finished,n,round(time()-started,digits=2))
            end
        end
    end
    reduce(+,counters)
end

function run_city(city,cache_root)
    start=time();folder=joinpath(B,"raw","full_network_city_profiles",city);mkpath(folder)
    root=joinpath(B,"inputs","city_networks_walk_connected")
    profiles=JSON.parsefile(joinpath(PROVENANCE,"osm_distance_demand_sources.json"))["cities"]
    profile=only(filter(x->x["city"]==city,profiles))
    owner=only(filter(x->x["city"]==city,JSON.parsefile(joinpath(PROVENANCE,"osm_car_ownership_sources.json"))["cities"]))
    p=BigInt(owner["numerator"])//BigInt(owner["denominator"])
    units=BigInt.(profile["mass_units"]);mass=units .// sum(units)
    bounds=Float64.(profile["upper_bounds_km"])
    g,info=MS.load_city(root,city;station_scope="gtfs_candidates")
    g=MS.with_car_costs(g,joinpath(root,"speed_overrides",city*".csv"))
    n=g.physical_nodes;total=n*(n-1)
    lookup=[category(bits) for bits in 0:65535]
    combined=zeros(Rational{BigInt},7);uniform=zeros(Rational{BigInt},7)
    classrows=NamedTuple[];histrows=NamedTuple[];bandrows=NamedTuple[];cachemeta=Dict()
    for (class,available,weight) in [("car",true,p),("noncar",false,1-p)]
        fixed,provenance=FD.load_fixed_graph(joinpath(root,"graphs",city),
            [joinpath(B,"inputs","geometry","road_lengths",city*".csv"),joinpath(B,"inputs","geometry","transit_lengths",city*".csv")];
            car_available=available,car_override=joinpath(root,"speed_overrides",city*".csv"))
        hash=FD.fingerprint(fixed);cache=joinpath(cache_root,city,"paths_"*class)
        if !isfile(joinpath(cache,"complete.json"))
            println(city," building optional reference cache ",class);flush(stdout)
            FD.cache_all_pairs(fixed,cache;progress=(d,n)->begin
                if d%1000==0 || d==n
                    println(city," fixed reference ",class," ",d,"/",n);flush(stdout)
                end
            end)
        end
        marker=JSON.parsefile(joinpath(cache,"complete.json"))
        @assert marker["graph_sha256"]==hash && marker["physical_nodes"]==n
        lengths=FD.matrix_file(joinpath(cache,"lengths_km.f64"),Float64,n)
        graph=available ? g : MS.without_car(g)
        progress=function(done,total,seconds)
            println(city," ",class," origins ",done,"/",total,"; ",seconds,"s");flush(stdout)
            open(joinpath(folder,"status.json"),"w") do io
                JSON.print(io,Dict("class"=>class,"completed_origins"=>done,"origins"=>total,"seconds"=>seconds,"updated"=>string(now())))
            end
        end
        println(city," routing ",class,"; ",n," locations; ",info["available_bike_stations"]," existing sharing stations");flush(stdout)
        timer=time();h=count_hist(graph,lengths,bounds,progress);routing=time()-timer
        counts=vec(sum(h,dims=1));@assert sum(counts)==total
        vals=zeros(Rational{BigInt},7);raw=zeros(Int64,7)
        for band in eachindex(mass)
            counts[band]>0 || iszero(mass[band]) || error("Unsupported positive-mass distance band")
            push!(bandrows,(;city,class,band,original_pairs=counts[band],target_mass=string(mass[band]),class_mass=string(weight)))
            for j in 1:65536
                c=h[j,band];c==0 && continue
                !available && ((j-1)&0xcccc)!=0 && error("Non-car route uses car")
                f=lookup[j];raw[f]+=c
                vals[f]+=mass[band]*(BigInt(c)//BigInt(counts[band]))
                push!(histrows,(;city,class,band,optimal_mode_set_bits=j-1,family=NAMES[f],pairs=c))
            end
        end
        @assert sum(vals)==1 && sum(raw)==total && raw[7]==0
        combined .+= weight.*vals;uniform .+= weight.*(BigInt.(raw).//BigInt(total))
        for k in eachindex(NAMES)
            push!(classrows,(;city,class,car_available=available,physical_nodes=n,ordered_pairs=total,
                stations=info["available_bike_stations"],family=NAMES[k],conditional_share=Float64(vals[k]),
                conditional_share_exact=string(vals[k]),class_mass_exact=string(weight),
                original_total_share=Float64(weight*vals[k]),uniform_conditional_share=raw[k]/total,routing_seconds=routing))
        end
        cachemeta[class]=merge(marker,Dict("length_sources"=>provenance))
        lengths=nothing;fixed=nothing;graph=nothing;h=nothing;GC.gc()
    end
    @assert sum(combined)==sum(uniform)==1 && combined[7]==0
    CSV.write(joinpath(folder,"class_shares.csv"),classrows)
    CSV.write(joinpath(folder,"band_mode_counts.csv"),histrows)
    CSV.write(joinpath(folder,"distance_bands.csv"),bandrows)
    fingerprints=Dict()
    for file in [joinpath(root,"graphs",city,"metadata.json"),joinpath(root,"graphs",city,"edges.csv"),
        joinpath(root,"speed_overrides",city*".csv"),joinpath(@__FILE__),joinpath(SRC, "demand", "osm_mode_share_analysis.jl"),
        joinpath(SRC, "demand", "osm_fixed_network_demand.jl")]
        fingerprints[relpath(file,B)]=bytes2hex(sha256(read(file)))
    end
    result=merge(info,Dict("completed"=>true,"ordered_pairs_per_class"=>total,"class_od_cases"=>2total,
        "threads"=>Threads.nthreads(),"seconds"=>time()-start,"car_probability_exact"=>string(p),
        "demand_profile"=>profile,"fixed_reference_cache_metadata"=>cachemeta,"input_sha256"=>fingerprints,
        "classification"=>"Walking only if sole movement mode; one non-walking mode includes walking access/egress; >=2 non-walking modes is multi_powered; cross-family optimum ties separate",
        "station_scenario"=>"All existing retained GTFS-derived candidates open, no opening budget, no extra locations",
        "distance_weighting"=>"Full-network class-specific fixed-only minimum-cost reference path length, then survey mass uniform within bands; no geometric proxy",
        "shares"=>Dict(NAMES[k]=>Float64(combined[k]) for k in eachindex(NAMES)),
        "exact_shares"=>Dict(NAMES[k]=>string(combined[k]) for k in eachindex(NAMES)),
        "uniform_od_shares"=>Dict(NAMES[k]=>Float64(uniform[k]) for k in eachindex(NAMES))))
    open(joinpath(folder,"result.json"),"w") do io;JSON.print(io,result,2);end
    println("FINISHED ",city," ",JSON.json(result["shares"]));flush(stdout)
    result
end

function main()
    smoke()
    cache_root=length(ARGS)>0 ? abspath(ARGS[1]) : joinpath(B,"work","fixed_reference_cache")   # shared with prepare_instances.jl
    cities=length(ARGS)>1 ? ARGS[2:end] : ["aachen","cologne","heidelberg","karlsruhe","nuremberg","dresden"]
    for city in cities
        run_city(city,cache_root)
        rows=NamedTuple[];classrows=NamedTuple[]
        for key in ["aachen","cologne","heidelberg","karlsruhe","nuremberg","dresden"]
            f=joinpath(B,"raw","full_network_city_profiles",key,"result.json");isfile(f)||continue
            r=JSON.parsefile(f);shares=r["shares"]
            push!(rows,(;city=key,physical_nodes=r["physical_nodes"],ordered_pairs_per_class=r["ordered_pairs_per_class"],
                stations=r["available_bike_stations"],car_probability_exact=r["car_probability_exact"],seconds=r["seconds"],
                (Symbol(name*"_pct")=>100shares[name] for name in NAMES)...))
            append!(classrows,[NamedTuple{Tuple(propertynames(r))}(Tuple(r)) for r in CSV.File(joinpath(dirname(f),"class_shares.csv"))])
        end
        CSV.write(joinpath(B,"data","full_network_city_profiles.csv"),rows)
        CSV.write(joinpath(B,"data","full_network_city_profiles_by_class.csv"),classrows)
        GC.gc()
    end
end
main()
