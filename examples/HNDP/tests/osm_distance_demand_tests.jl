using Test,Random,JSON
include("../demand/osm_distance_demand_analysis.jl")

@testset "Distance demand construction" begin
    for city in ["aachen","cologne","heidelberg","karlsruhe","nuremberg","dresden"]
        p=OSMDemand.profile(city)
        @test OSMDemand.validate(p)===nothing
        @test sum(BigInt.(p["mass_units"]) .// sum(p["mass_units"]))==1
    end
    @test OSMDemand.band(1.0,[1.,3.])==1
    @test OSMDemand.band(1.00001,[1.,3.])==2
    @test OSMDemand.band(10.,[1.,3.])==3
    coords=[[50.,6.],[50.001,6.],[50.03,6.]]
    p=Dict("upper_bounds_km"=>[1.0],"mass_units"=>[2,1])
    pairs,meta=OSMDemand.sample_pairs(coords,30,MersenneTwister(1),p;detour_factor=1.)
    @test count(x->OSMDemand.distance_km(coords[x[1]],coords[x[2]])<=1,pairs)==20
    @test meta["band_user_counts"]==[20,10]
    @test all(x->x[1]!=x[2],pairs)
    @test OSMDemand.sample_pairs(coords,30,MersenneTwister(1),p;detour_factor=1.)[1]==pairs
    @test_throws ErrorException OSMDemand.sample_pairs(coords,3,MersenneTwister(1),p;eligible=(o,d)->false,max_attempts=30)
    @test_throws ArgumentError OSMDemand.sample_pairs(coords,3,MersenneTwister(1),p;detour_factor=0.)
end

@testset "Mode-family ties and powered combinations" begin
    lookup=demand_label_metrics()
    for bits in [0,2,4,8,16,32,64,256,512,20,12,1024,0x8000]
        counts=zeros(Int64,65536);counts[Int(bits)+1]=1
        s=OSMModeShares.summarize(counts,Dict{String,Any}())
        for k in 1:6
            @test lookup[k,Int(bits)+1]==s["family_counts"][DEMAND_METRICS[k]]
        end
        @test lookup[7,Int(bits)+1]==s["unreachable_pairs"]
        @test lookup[8,Int(bits)+1]==s["any_multimodal_optimum_including_walking"]
        @test lookup[9,Int(bits)+1]==s["any_multiple_powered_modes_optimum"]
    end
end

@testset "Band routing matches existing all-OD analysis" begin
    # One disconnected physical node; reachable paths use walk and car layers.
    g=OSMModeShares.mode_graph(3,[(1,2,20),(2,1,20),(1,4,1),(4,1,1),(2,5,1),(5,2,1),(4,5,1),(5,4,1)])
    no=OSMModeShares.without_car(g)
    coords=[[50.,6.],[50.001,6.],[50.03,6.]]
    p=Dict("upper_bounds_km"=>[1.0],"mass_units"=>[2,1])
    h,t,ds=count_distance_demand(g,no,coords,p;factors=[1.0])
    @test sum(t)==6
    @test sum(h[7,:,:,1])==4
    lookup=demand_label_metrics()
    for (c,graph) in enumerate((g,no))
        counts=OSMModeShares.count_mode_shares(graph;progress=false)
        for k in 1:9
            @test sum(h[k,:,:,c])==sum(counts .* lookup[k,:])
        end
    end
    @test sum(h[1:7,:,1,1])==6
    owner=79//125;mass=[2//3,1//3]
    @test sum(mass .* owner)==owner
    @test sum(mass .* (1-owner))==1-owner
    @test sum(mass[b]*(owner*sum(h[1:7,b,1,1])+(1-owner)*sum(h[1:7,b,1,2]))/t[b,1] for b in 1:2)==1
end
