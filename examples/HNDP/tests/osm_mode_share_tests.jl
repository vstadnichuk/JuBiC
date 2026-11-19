using Test, Random
include("../demand/osm_mode_share_analysis.jl")
using .OSMModeShares
const MS=OSMModeShares

# Independent Bellman-Ford distances and fixed-point propagation over tight
# edges; deliberately no heap or settlement-order assumptions.
function reference_modes(g,source)
    nv=4g.physical_nodes;inf=typemax(Int64);d=fill(inf,nv);d[source]=0
    for _ in 1:nv-1
        changed=false
        for u in 1:nv,e in g.offsets[u]:g.offsets[u+1]-1
            v=g.targets[e]
            if d[u]!=inf && d[u]+g.costs[e]<d[v]
                d[v]=d[u]+g.costs[e];changed=true
            end
        end
        changed || break
    end
    bits=zeros(UInt16,nv);bits[source]=1
    for _ in 1:16nv
        changed=false
        for u in 1:nv,e in g.offsets[u]:g.offsets[u+1]-1
            v=g.targets[e]
            d[u]==inf && continue
            d[u]+g.costs[e]==d[v] || continue
            propagated=UInt16(0)
            for m in 0:15
                bits[u]&(UInt16(1)<<m)==0 && continue
                propagated |= UInt16(1)<<(m|Int(g.flags[e]))
            end
            fresh=bits[v]|propagated
            fresh!=bits[v] && (bits[v]=fresh;changed=true)
        end
        changed || return d,bits
    end
    error("reference failed to converge")
end

@testset "OSM shortest-path mode shares" begin
    @testset "all bitset transitions" begin
        for bits in UInt16[0,1,2,0x5555,0xabcd,0xffff],flag in UInt8[0,1,2,4,8]
            expected=UInt16(0)
            for mask in 0:15
                bits&(UInt16(1)<<mask)!=0 && (expected |= UInt16(1)<<(mask|Int(flag)))
            end
            @test MS.add_mode(bits,flag)==expected
        end
    end
    @testset "ties through a zero-cost transit exit" begin
        g=MS.mode_graph(2,[(1,2,3),(1,5,1),(5,6,2),(6,2,0),(1,3,1),(3,4,1),(4,2,1)])
        w=MS.shortest_modes!(MS.Workspace(g),g,1)
        @test w.distance[2]==3
        @test w.labels[2]==0x0016 # walking, car, transit all equally good
        counts=MS.count_mode_shares(g;workers=1,progress=false)
        @test sum(counts)==2
        @test counts[1]==1 # reverse pair is unreachable
        info=MS.summarize(counts,Dict{String,Any}())
        @test info["family_counts"]["tie_between_families"]==1
        @test info["unreachable_pairs"]==1
        @test_throws ArgumentError MS.mode_graph(2,[(1,2,-1)])
        @test_throws ArgumentError MS.mode_graph(2,[(1,2,0)])
        @test_throws ArgumentError MS.mode_graph(2,[(1,2,1),(1,2,2)])
    end
    @testset "distinct exclusive and inclusive CSV columns" begin
        counts=zeros(Int64,65536)
        counts[(1<<2)+(1<<6)+1]=1 # tied car and car+transit optima
        info=Dict{String,Any}("city"=>"fixture","physical_nodes"=>2,
            "available_bike_stations"=>2,"analysis_seconds"=>0)
        row=MS.summary_row(MS.summarize(counts,info))
        @test row.multiple_powered_modes_share==0
        @test row.any_multiple_powered_modes_share==1
        @test row.tie_between_families_share==1
    end
    @testset "random graphs against independent reference" begin
        rng=MersenneTwister(27)
        for trial in 1:40
            n=3;arcs=Tuple{Int,Int,Int}[]
            for layer in 0:3,i in 1:n,j in 1:n
                i!=j && rand(rng)<0.7 && push!(arcs,(i+layer*n,j+layer*n,rand(rng,1:5)))
            end
            for layer in 1:3,i in 1:n
                push!(arcs,(i,i+layer*n,rand(rng,1:3)))
                push!(arcs,(i+layer*n,i,layer==2 ? 0 : rand(rng,1:3)))
            end
            g=MS.mode_graph(n,arcs);w=MS.Workspace(g)
            for source in 1:n
                distances,bits=reference_modes(g,source)
                MS.shortest_modes!(w,g,source)
                @test w.distance==distances
                @test w.labels==bits
            end
            serial=MS.count_mode_shares(g;workers=1,progress=false)
            threaded=MS.count_mode_shares(g;progress=false)
            @test serial==threaded
            @test sum(serial)==n*(n-1)
        end
    end
    @testset "exact ownership weights and car exclusion" begin
        owner=zeros(Int64,65536);nonowner=copy(owner)
        owner[(1<<2)+1]=2;nonowner[(1<<1)+1]=2
        mix=MS.ownership_mix(owner,nonowner,632//1000)
        @test mix.owner_units//(mix.owner_units+mix.nonowner_units)==632//1000
        @test mix.weight_denominator==125
        @test mix.weight_units[(1<<2)+1]==158
        @test mix.weight_units[(1<<1)+1]==92
        @test MS.ownership_mix(owner,nonowner,0//1).weight_units==nonowner
        @test MS.ownership_mix(owner,nonowner,1//1).weight_units==owner
        @test_throws ArgumentError MS.ownership_mix(owner,owner,1//2)
        @test_throws ArgumentError MS.ownership_mix(owner,nonowner,101//100)
        g=MS.mode_graph(2,[(1,2,10),(2,1,10),(1,3,1),(3,4,1),(4,2,1)])
        no_car=MS.without_car(g)
        @test all(f->f!=2,no_car.flags)
        w=MS.shortest_modes!(MS.Workspace(no_car),no_car,1)
        @test w.distance[2]==10
        @test w.labels[2]==1<<1
        @test w.distance[3]==typemax(Int64)
        @test MS.shortest_modes!(MS.Workspace(g),g,1).distance[2]==3
        counts=MS.count_mode_shares(no_car;progress=false)
        @test all(counts[i]==0 || ((i-1)&0xcccc)==0 for i in eachindex(counts))
        mktempdir() do dir
            file=joinpath(dir,"costs.csv")
            write(file,"source,target,cost_cents\n3,4,20\n")
            updated=MS.with_car_costs(g,file)
            @test MS.shortest_modes!(MS.Workspace(updated),updated,1).distance[2]==10
            write(file,"source,target,cost_cents\n")
            @test_throws ArgumentError MS.with_car_costs(g,file)
            GC.gc()
        end
    end
    @testset "spatial station matching matches the importer" begin
        rng=MersenneTwister(31)
        coords=[[50+0.03rand(rng),6+0.03rand(rng)] for _ in 1:100]
        stops=[[50+0.04rand(rng),6+0.04rand(rng)] for _ in 1:80]
        expected=Set{Int}()
        for p in stops
            distances=[MS.haversine(p[1],p[2],q[1],q[2]) for q in coords]
            minimum(distances)<=0.25 && push!(expected,argmin(distances))
        end
        mktempdir() do dir
            path=joinpath(dir,"stops.txt")
            open(path,"w") do io
                println(io,"stop_id,stop_lat,stop_lon")
                for (i,p) in enumerate(stops);println(io,"$i,$(p[1]),$(p[2])");end
            end
            @test MS.candidate_stations(coords,path)==sort!(collect(expected))
            GC.gc()
        end
    end
end
