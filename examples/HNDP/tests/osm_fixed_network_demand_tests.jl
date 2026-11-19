using Test
include("../demand/osm_fixed_network_demand.jl")
const FD=OSMFixedNetworkDemand
@testset "Fixed cost paths retain physical path length" begin
    edges=[(1,2,1,4.0),(1,4,1,0.0),(4,2,0,1.0),(2,1,2,2.0),
           (2,3,1,2.0),(3,2,1,3.0),(3,1,3,5.0),(1,3,9,0.2)]
    g=FD.FixedGraph(3,edges)
    w=FD.Workspace(g);FD.shortest_paths!(w,g,1)
    @test w.cost[2]==1 && w.length[2]==1 # Equal cost → shorter physical length.
    @test w.cost[3]==2 && w.length[3]==3 # Cost-optimal ≠ physically shortest.
    # Windows can retain mapped-file handles until process exit.
    let folder=mktempdir(;cleanup=false)
        costs,lengths=FD.cache_all_pairs(g,folder)
        @test costs[2,1]==1 && costs[1,2]==2
        @test lengths[3,1]==3
        access,egress=FD.nearest_sites(costs,lengths,[1,3])
        @test access==[1,2,2]
        @test egress==[1,1,2] # Destination mapping has opposite direction.
        result=FD.aggregate_demand(costs,lengths,[1,3],[1.,3.],[2,5,3];car_probability=71//100)
        @test result.band_counts==[1,4,1]
        @test result.total==[1//5 1//8;11//20 1//8]
        @test result.diagonal_mass==13//40
        @test result.inter_site_mass==27//40
        @test sum(result.owner)==71//100 && sum(result.nonowner)==29//100
        @test all(result.owner.==result.total.*(71//100))
        @test sum(result.band_aggregate_pair_counts)==6 # Original self-pairs excluded.
        loaded_costs,loaded_lengths=FD.load_cache(folder;expected_graph_hash=FD.fingerprint(g))
        @test loaded_costs==costs && loaded_lengths==lengths
        @test_throws ArgumentError FD.load_cache(folder;expected_graph_hash="wrong")
        @test_throws ArgumentError FD.aggregate_demand(costs,lengths,[1,3],[0.,1.,3.],[1,1,1,1])
        @test_throws ArgumentError FD.nearest_sites(costs,lengths,Int[])
        # Release Windows mappings before mktempdir cleanup.
        costs=nothing;lengths=nothing;loaded_costs=nothing;loaded_lengths=nothing;GC.gc()
    end
    @test_throws ArgumentError FD.FixedGraph(2,[(1,2,-1,1.0)])
    @test_throws ArgumentError FD.FixedGraph(2,[(1,2,1,-1.0)])
end
