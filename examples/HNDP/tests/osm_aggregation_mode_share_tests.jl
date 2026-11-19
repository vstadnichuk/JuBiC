include("../demand/osm_aggregation_mode_share_analysis.jl")
using Test
@testset "Existing aggregation and exact OD weighting" begin
    # Every layer has a bidirectional 1--2--3 chain. Only node 2 can board
    # transit; discarding it preserves pure-mode distances but loses access.
    n=3;arcs=Tuple{Int,Int,Int}[]
    for layer in 0:3,(a,b) in ((1,2),(2,1),(2,3),(3,2))
        push!(arcs,(a+layer*n,b+layer*n,layer==2 ? 1 : 20))
    end
    append!(arcs,[(2,8,1),(8,2,0)])
    g=OSMModeShares.mode_graph(n,arcs);graphs,costs=split_layers(g)
    ag=aggregate_existing(g,[1,3],graphs,costs)
    @test ag.physical_nodes==2
    @test length(ag.offsets)==9
    a=OSMModeShares.Workspace(g);b=OSMModeShares.Workspace(ag)
    OSMModeShares.shortest_modes!(a,g,1);OSMModeShares.shortest_modes!(b,ag,1)
    @test a.distance[3]==40 # no exit at node 3; boarding at 2 cannot help
    @test b.distance[2]==40
    # Add a transit exit at node 3: full can board at 2 and exit at 3.
    push!(arcs,(9,3,0));g=OSMModeShares.mode_graph(n,arcs)
    graphs,costs=split_layers(g);ag=aggregate_existing(g,[1,3],graphs,costs)
    a=OSMModeShares.Workspace(g);b=OSMModeShares.Workspace(ag)
    OSMModeShares.shortest_modes!(a,g,1);OSMModeShares.shortest_modes!(b,ag,1)
    @test a.distance[3]==22
    @test b.distance[2]==40
    @test b.distance[2]>a.distance[3]
    profile=Dict("upper_bounds_km"=>[1.0,3.0],"mass_units"=>[30,40,30])
    coords=[[50.,6.],[50.02,6.]]
    labels=fill(UInt16(2),2,2,2) # walking mask=1, encoded as bit 2
    h,t=selected_hist(labels,coords,profile)
    rows=share_rows(h,t,profile,BigInt(7)//BigInt(10),"fixture","sharing","aggregated",2)
    @test sum(t)==2
    uniform=filter(r->r.weighting=="uniform_od",rows)
    @test only(filter(r->r.metric=="walking",uniform)).share==1.0
    @test all(r.supported for r in uniform)
    @test all(!r.supported && isnan(r.share) for r in rows if r.weighting=="distance_weighted")
    medoids,sizes=ExistingAggregation.balanced_medoids([abs(i-j) for i in 1:9,j in 1:9],4)
    @test length(unique(medoids))==4
    @test sum(sizes)==9
    @test maximum(sizes)-minimum(sizes)<=1
end

@testset "Sharing-only active setting and cache compatibility" begin
    @test AGG_BIKE_SETTING=="sharing"
    @test AGG_PRIMARY_WEIGHTING=="distance_weighted"
    @test !sharing_marker_compatible(Dict("city"=>"fixture"))
    @test !sharing_marker_compatible(Dict("analysis_schema"=>2,"bike_setting"=>"private_bike","primary_weighting"=>"distance_weighted"))
    @test sharing_marker_compatible(Dict("analysis_schema"=>5,"cost_configuration"=>AGG_COST_CONFIGURATION,"selection_metric"=>"full_multimodal_all_existing_sharing_open","bike_setting"=>"sharing","primary_weighting"=>"distance_weighted","targets"=>AGG_TARGETS,"walking_connectivity_policy"=>"largest_walking_scc"))
end

@testset "Walking connectivity cleanup" begin
    edges=[(1,2),(2,1),(2,3),(3,2),(4,5),(5,4)]
    @test OSMConnectivity.walking_backbone_nodes(6,edges)==[1,2,3]
    @test OSMConnectivity.walking_backbone_nodes(4,[(1,2),(2,1),(3,4),(4,3)])==[1,2]
    @test OSMConnectivity.walking_backbone_nodes(3,[(1,2),(2,1),(2,3)])==[1,2]
    @test_throws ArgumentError OSMConnectivity.walking_backbone_nodes(2,[(1,2)])
    coords=collect(1:5);road=Dict((u,v)=>10 for (u,v) in edges)
    co,car,bike,walk,keep=OSMConnectivity.restrict_road_layers(coords,road,road,road)
    @test co==keep==[1,2,3]
    @test all(1<=u<=3 && 1<=v<=3 for (u,v) in keys(car))
    @test length(walk)==4
end
