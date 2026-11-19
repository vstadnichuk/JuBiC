using Test, Graphs, SparseArrays
include("../demand/osm_aggregation_helpers.jl")
using .OSMAggregation: multimodal_site_metric, potential_station_sequence,
    balanced_medoids, save_selection_metric, load_selection_metric
@testset "Complete multimodal selection and nested candidates" begin
    # Physical1->physical2 uses transit (nodes3,4), cheaper than direct travel.
    g=DiGraph(4);c=spzeros(Int,4,4)
    for (u,v,w) in [(1,2,100),(2,1,100),(1,3,3),(3,4,2),(4,2,0),(2,4,3),(4,3,2),(3,1,0)]
        add_edge!(g,u,v);c[u,v]=w
    end
    metric=multimodal_site_metric(g,c,[1,2])
    @test metric==[0 5;5 0]
    # Network medoid2 differs from the old first-site/geographic center choices.
    m=[0 2 7 9;2 0 5 7;7 5 0 2;9 7 2 0]
    order=potential_station_sequence(m,4)
    @test order==[2,4,1,3]
    @test potential_station_sequence(m,2)==order[1:2]
    @test potential_station_sequence(m,3)==order[1:3]
    @test Set(potential_station_sequence(m,2;eligible=[1,3,4])) ⊆ Set([1,3,4])
    @test isempty(potential_station_sequence(m,0))
    @test_throws ArgumentError potential_station_sequence(m,5)
    # The former capacity-index confusion failed when capacity>cluster count.
    line=[abs(i-j) for i in 1:6,j in 1:6]
    medoids,sizes=balanced_medoids(line,2)
    @test sizes==[3,3] && sort(medoids)==[2,5]
    @test balanced_medoids(line,1)==([3],[6])
    mktempdir() do folder
        save_selection_metric(folder,4,m)
        @test load_selection_metric(folder,4)==m
    end
end
