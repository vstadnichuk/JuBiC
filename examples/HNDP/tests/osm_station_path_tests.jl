using Test
if !isdefined(@__MODULE__, :build_hndp_station_path_instance)
    include("../models/hndp_model_generation.jl")
end

"""
    station_path_counterexample()

Station s can be entered from the walking layer (node 2) and from a transit stop
(node 3); both sharing routes return the bike at station t. Both routes require
{s, t}, but they use different entry arcs. Decision-arc dominance keeps both,
station dominance removes the more expensive one. Fixed path 1 -> 2 costs 10.

Nodes: 1 origin, 2 walking node at s, 3 transit stop, 4 bike node s,
5 bike node t, 6 walking node at t, 7 destination.
"""
function station_path_counterexample()
    graph = DiGraph(7); cost = spzeros(Float64, 7, 7); risk = spzeros(Float64, 7, 7)
    for (u, v, c) in [(1, 7, 10.0), (1, 2, 2.0), (2, 4, 1.0), (1, 3, 1.0), (3, 4, 3.0),
                      (4, 5, 1.0), (5, 6, 1.0), (6, 7, 1.0), (6, 5, 1.0)]
        add_edge!(graph, u, v); cost[u, v] = c
    end
    risk[2, 4] = -1.0; risk[3, 4] = -1.0
    user = User("u", 1, 7, risk, cost, nothing, nothing)
    groups = [[(2, 4), (3, 4)], [(6, 5), (5, 6)]]
    arcs = reduce(vcat, groups)
    return HNDPwC(graph, [user], arcs, Dict(a => 0.0 for a in arcs), nothing, groups, 2), user
end

function hndp_station_path_tests(; solver_factory=() -> GurobiSolver())
    @testset "Station data" begin
        h, _ = station_path_counterexample()
        data = hndp_station_data(h)
        @test length(data.stations) == 2 && all(data.budgeted)
        @test data.arc_station[(2, 4)] == data.arc_station[(3, 4)] == 1
        @test data.arc_station[(6, 5)] == data.arc_station[(5, 6)] == 2
        @test _hndp_station_anchor(data, 1) == (2, 4)
        ungrouped = HNDPwC(h.mygraph, h.users, vcat(h.edgeA, [(1, 2)]), merge(h.edge_price, Dict((1, 2) => 0.0)), nothing, h.decision_groups, 2)
        d2 = hndp_station_data(ungrouped)
        @test length(d2.stations) == 3 && d2.budgeted == [true, true, false]
        @test_throws ArgumentError hndp_station_data(HNDPwC(h.mygraph, h.users, h.edgeA, h.edge_price, nothing, [[(2, 4)], [(2, 4)]], 2))
    end

    @testset "Station dominance removes the counterexample path" begin
        h, u = station_path_counterexample()
        raw, _, t1 = _enumerate_user_paths_station_dominance(u, h, 10.0, time() + 60; use_station_dominance=false, check_decision_budget=false)
        legacy, _, t2 = _enumerate_user_paths_v2(u, h, 10.0, time() + 60; use_decision_arc_dominance=true)
        station, _, t3 = _enumerate_user_paths_station_dominance(u, h, 10.0, time() + 60)
        @test !t1 && !t2 && !t3
        @test sort(getfield.(raw, :cost)) == [6.0, 7.0, 10.0]
        @test sort(getfield.(legacy, :cost)) == [6.0, 7.0, 10.0]
        @test sort(getfield.(station, :cost)) == [6.0, 10.0]
        # Without dominance, the new enumerator returns exactly the legacy raw path set.
        legacy_raw, _, _ = _enumerate_user_paths_v2(u, h, 10.0, time() + 60; use_decision_arc_dominance=false, check_decision_budget=false)
        @test Set(Tuple(p.arcs) for p in raw) == Set(Tuple(p.arcs) for p in legacy_raw)
        # Budget 1 prunes both sharing routes (each requires two stations).
        h1 = HNDPwC(h.mygraph, h.users, h.edgeA, h.edge_price, nothing, h.decision_groups, 1)
        pruned, _, _ = _enumerate_user_paths_station_dominance(u, h1, 10.0, time() + 60)
        @test getfield.(pruned, :cost) == [10.0]
    end

    @testset "Station-based cleanup of complete paths" begin
        h, _ = station_path_counterexample()
        data = hndp_station_data(h)
        make(arcs, cost; risk=0.0) = HNDPPath(arcs, cost, risk, 0.0)
        cheap = make([(1, 2), (2, 4), (4, 5), (5, 6), (6, 7)], 6.0)
        other_entry = make([(1, 3), (3, 4), (4, 5), (5, 6), (6, 7)], 7.0)
        tie = make([(1, 3), (3, 4), (4, 5), (5, 6), (6, 7)], 6.0; risk=-5.0)
        fixed = make([(1, 7)], 10.0)
        @test _cleanup_station_dominated_hndp_paths([other_entry, cheap, fixed], data) == [cheap, fixed]
        @test _cleanup_station_dominated_hndp_paths([cheap, tie, fixed], data) == [cheap, tie, fixed]
        @test _cleanup_station_dominated_hndp_paths([fixed], data) == [fixed]
        # The decision-arc cleanup keeps the path with the other entry arc.
        @test length(_cleanup_strictly_dominated_hndp_paths([other_entry, cheap, fixed], h.edgeA)) == 3
    end

    @testset "Shared completion bounds" begin
        h, u = station_path_counterexample()
        cache = _station_completion_bounds(h)
        @test length(cache) == 1
        cached = _station_cached_completion_bounds(cache, u, h)
        @test cached == _hndp_reverse_dijkstra_distances(h.mygraph, u.mcost, u.destination)
        # A user with another cost matrix is not served from the cache.
        other = User("v", 1, 7, u.mrisk, copy(u.mcost), nothing, nothing)
        @test isnothing(_station_cached_completion_bounds(cache, other, h))
        @test isnothing(_station_cached_completion_bounds(nothing, u, h))
        data = hndp_station_data(h)
        a, _, _ = _enumerate_user_paths_station_dominance(u, h, 10.0, time() + 60; station_data=data)
        b, _, _ = _enumerate_user_paths_station_dominance(u, h, 10.0, time() + 60; station_data=data, heuristic_distances=cached)
        @test Set(Tuple(p.arcs) for p in a) == Set(Tuple(p.arcs) for p in b)
        # Batched completion bounds give the same model as the shared cache.
        for batch in (1, 2)
            i1, _, c1 = build_hndp_station_path_instance(h, solver_factory(); enumeration_time_limit=60.0)
            i2, _, c2 = build_hndp_station_path_instance(h, solver_factory(); enumeration_time_limit=60.0,
                completion_bound_batch=batch)
            @test c1 == c2
            m1 = _hndp_master_model(i1); m2 = _hndp_master_model(i2)
            @test num_variables(m1) == num_variables(m2)
        end
    end

    @testset "Station-based model: objective, size, and solution export" begin
        h, _ = station_path_counterexample()
        for budget in (0, 1, 2), dominance in (false, true)
            hb = HNDPwC(h.mygraph, h.users, h.edgeA, h.edge_price, nothing, h.decision_groups, budget)
            inst, _, counts = build_hndp_station_path_instance(hb, solver_factory();
                enumeration_time_limit=60.0, use_station_dominance=dominance)
            m = _hndp_master_model(inst); set_silent(m); optimize!(m)
            @test termination_status(m) == MOI.OPTIMAL
            @test isapprox(objective_value(m), budget == 2 ? -1.0 : 0.0; atol=1e-8)
            @test counts["u"] == (budget < 2 ? 1 : dominance ? 2 : 3)
            # One binary opening variable per station; none for fixed arcs.
            @test count(is_binary, all_variables(m)) == 2
            @test haskey(object_dictionary(m), :x) && length(m[:x]) == 2
            @test all(g -> m[:x][first(g)] isa VariableRef, hb.decision_groups)
            weights = hndp_path_solution_weights(m)
            @test isapprox(sum(weights["u"]), 1.0; atol=1e-8)
            # Same optimum as the legacy path builder.
            legacy, _, _ = build_hndp_path_instance(hb, solver_factory(); enumeration_time_limit=60.0,
                parallelize=true, use_decision_arc_dominance=dominance)
            lm = _hndp_master_model(legacy); set_silent(lm); optimize!(lm)
            @test isapprox(objective_value(lm), objective_value(m); atol=1e-8)
        end
    end

    @testset "Station linking aggregates all paths of a user per station" begin
        # Two paths share the entry arc of station 1; a second user needs station 2.
        graph = DiGraph(9); cost = spzeros(Float64, 9, 9); risk = spzeros(Float64, 9, 9)
        for (u, v, c) in [(1, 4, 10.0), (1, 2, 3.0), (2, 3, 1.0), (2, 9, 1.0), (3, 4, 1.0), (9, 4, 2.0),
                          (5, 7, 10.0), (5, 6, 3.0), (6, 7, 2.0)]
            add_edge!(graph, u, v); cost[u, v] = c
        end
        risk[1, 2] = -1.0; risk[5, 6] = -1.0
        users = [User("u1", 1, 4, risk, cost, nothing, nothing), User("u2", 5, 7, risk, cost, nothing, nothing)]
        groups = [[(1, 2)], [(5, 6)]]
        arcs = reduce(vcat, groups)
        h = HNDPwC(graph, users, arcs, Dict(a => 0.0 for a in arcs), nothing, groups, 1)
        inst, _, counts = build_hndp_station_path_instance(h, solver_factory(); enumeration_time_limit=60.0,
            use_station_dominance=false)
        m = _hndp_master_model(inst)
        @test counts["u1"] == 3
        links = [c for c in all_constraints(m; include_variable_in_set_constraints=false) if startswith(name(c), "station_link_u1")]
        @test length(links) == 1
        relax_integrality(m); set_silent(m); optimize!(m)
        # Aggregated linking: LP bound equals the integer optimum -1 (per-path linking would give -1.5).
        @test isapprox(objective_value(m), -1.0; atol=1e-8)
    end
end
