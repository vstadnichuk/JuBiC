using Graphs, JuBiC


struct TestLabel
    node::Any
    costs::Any
    distance::Any
end


function neighboring_function(graph, label::TestLabel, costm, distm, maxDist)
    neighbors_list = []
    for n in neighbors(graph, label.node)
        new_label = TestLabel(
            n,
            label.costs + costm[label.node, n],
            label.distance + distm[label.node, n],
        )
        if new_label.distance <= maxDist
            push!(neighbors_list, new_label)
        end
    end
    return neighbors_list
end

function test_labeling()
    graph = DiGraph(4)
    add_edge!(graph, 1, 2)
    add_edge!(graph, 2, 4)
    add_edge!(graph, 1, 3)
    add_edge!(graph, 3, 4)
    add_edge!(graph, 2, 3)

    edge_costs = [0 3 7 0; 0 0 2 8; 0 0 0 4; 0 0 0 0]
    distances = floyd_warshall_shortest_paths(graph, edge_costs).dists

    starts = TestLabel(1, 0, 0)
    ends = TestLabel(4, 0, 0)

    always_false(label::TestLabel) = false
    isgoal(label::TestLabel) = label.node == ends.node
    dominance(lv, lw) =
        (lv.me.node == lw.me.node) &&
        (lv.cost + lv.heu >= lw.cost + lw.heu) &&
        (lv.me.distance > lw.me.distance)
    heuristicf_none(lv::TestLabel, lw::TestLabel) = 0
    heuristicf_optimal(lv::TestLabel, lw::TestLabel) = distances[lv.node, lw.node]
    costf(lv::TestLabel, lw::TestLabel) = edge_costs[lv.node, lw.node]
    neighboring(label::TestLabel) =
        neighboring_function(graph, label, edge_costs, zeros(Float32, 4, 4), 10)


    @testset "shortest_path_labeling Tests" begin
        solution_1 = JuBiC.shortest_path_labeling(
            neighboring,
            starts,
            ends,
            heuristicf_optimal,
            costf,
            isgoal,
            dominance,
            12,
            Inf,
        )
        @test solution_1.status == JuBiC.LS_OPTIMAL
        @test solution_1.goal.node == ends.node
        @test solution_1.costs == 9

        solution_2 = JuBiC.shortest_path_labeling(
            neighboring,
            starts,
            ends,
            heuristicf_none,
            costf,
            isgoal,
            dominance,
            8,
            Inf,
        )
        @test solution_2.status == JuBiC.LS_NO_SOLUTION

        solution_3 = JuBiC.shortest_path_labeling(
            neighboring,
            starts,
            ends,
            heuristicf_none,
            costf,
            always_false,
            dominance,
            12,
            Inf,
        )
        @test solution_3.status == JuBiC.LS_NO_SOLUTION

        solution_4 = JuBiC.shortest_path_labeling(
            neighboring,
            starts,
            ends,
            heuristicf_none,
            costf,
            isgoal,
            dominance,
            12,
            -1,
        )
        @test solution_4.status == JuBiC.LS_TIMEOUT
    end

    @testset "enumerate_all_paths Tests" begin
        solution_1 = JuBiC.enumerate_all_paths(
            neighboring,
            starts,
            ends,
            costf,
            isgoal,
            dominance,
            12,
            Inf,
        )
        @test solution_1.status == JuBiC.LS_OPTIMAL
        @test length(solution_1.goal) == 3

        solution_2 = JuBiC.enumerate_all_paths(
            neighboring,
            starts,
            ends,
            costf,
            isgoal,
            dominance,
            8,
            Inf,
        )
        @test solution_2.status == JuBiC.LS_NO_SOLUTION

        solution_3 = JuBiC.enumerate_all_paths(
            neighboring,
            starts,
            ends,
            costf,
            isgoal,
            dominance,
            12,
            -1,
        )
        @test solution_3.status == JuBiC.LS_TIMEOUT
    end
end


test_labeling()
