using JuBiC, JuMP, Gurobi, Test
import MathOptInterface as MOI

@testset "v0.6 branch integration" begin
    wrapper = GurobiSolver()
    output = JuBiC.repo_local_tempdir("tests", "integration"; prefix="integration")
    @testset "Positional constructor compatibility" begin
        blc = BLCparam(wrapper, false, output, "lp", JuBiC.RunStats(), 60, 42, 1, 1, false)
        @test blc.parallel_subsolvers == 1
        @test !blc.integer_obj
        integer_blc = BLCparam(wrapper, false, output, "lp", JuBiC.RunStats(), 60, 42, 1, 1, false, true)
        @test integer_blc.integer_obj
        lag = BlCLagparam(wrapper, false, output, "lp", JuBiC.RunStats(), 60, 42, 1, 1,
            false, PARETO_OPTIMALITY_ONLY, true, 1e9, 2e-4)
        @test lag.blc_pareto_band_tolerance == 2e-4
        @test lag.parallel_subsolvers == 1
    end

    @testset "Cost-aware callbacks and integer cuts" begin
        instance = generate_blc_simple_bilevel_instance()
        master = instance.master
        calls = Ref(0)
        callback = (a, sub, cost::Real, x::AbstractDict) -> begin
            calls[] += 1
            @test isfinite(cost)
            @test haskey(x, a)
            master.big_m(a, sub)
        end
        dynamic_master = JuBiC.BlCMaster(master.hpr, master.A, master.link_vars,
            callback, master.sub_names, master.sub_objectives)
        params = BLCparam(wrapper, false, output, "lp", 60, true)
        results = solve_instance!(Instance(dynamic_master, instance.subproblems), params)
        @test calls[] > 0
        @test isapprox(results.data["Opt"], 1; atol=2e-4)
        @test JuBiC._evaluate_blc_big_m((a, sub) -> 3, 1, "follower", 2.0, Dict(1 => 1.0)) == 3
        @test_throws DomainError JuBiC._evaluate_blc_big_m((a, sub) -> -1, 1, "follower", 2.0, Dict())
        @test JuBiC._round_integer_objective_if_close(0.2, params, "follower", "test") == 0.2
        @test JuBiC._round_integer_objective_if_close(1.000001, params, "follower", "test") == 1.0
    end

    @testset "Linking-variable branching priorities" begin
        model = Model(() -> JuBiC.get_next_optimizer(wrapper))
        @variable(model, linking[1:2], Bin)
        @variable(model, other, Bin)
        links = Dict(i => linking[i] for i in 1:2)
        params = GBCparam(wrapper, false, output, "lp")
        params.stats.data["branching_rule"] = "linking_first"
        JuBiC._apply_linking_branching_rule!(model, links, params)
        optimize!(model)
        attr = Gurobi.VariableAttribute("BranchPriority")
        priority(v) = MOI.get(unsafe_backend(model), attr, JuMP.optimizer_index(v))
        @test all(priority(v) > priority(other) for v in linking)
        params.stats.data["branching_rule"] = "fixed_linking_order"
        JuBiC._apply_linking_branching_rule!(model, links, params)
        optimize!(model)
        @test priority(linking[1]) > priority(linking[2]) > priority(other)
        params.stats.data["branching_rule"] = "misspelled_rule"
        @test_throws ArgumentError JuBiC._apply_linking_branching_rule!(model, links, params)
    end

    if Threads.nthreads() >= 2
        @testset "Two concurrent follower workers" begin
            results = test_gbc_two_follower(; parallel_subsolvers=2)
            @test results.data["parallel_connector_workers_used"] == 2
            @test isapprox(results.data["Opt"], 0; atol=2e-4)
        end
    end
end
