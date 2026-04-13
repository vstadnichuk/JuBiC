using Test
using JuBiC
using JuMP
using Gurobi
using Graphs

include("hndp_model_generation_v2.jl")

"""
    build_toy_hndp_two_users()

Build a small deterministic HNDP instance with two independent user OD-pairs.
Each user has a fixed-network path and a cheaper decision-dependent path. This
lets us test both big-M rules while keeping the instance tiny and reproducible.
"""
function build_toy_hndp_two_users()
    graph = DiGraph(8)
    fixed_arcs = [(1, 2), (2, 4), (5, 6), (6, 8)]
    decision_arcs = [(1, 3), (3, 4), (5, 7), (7, 8)]
    for arc in vcat(fixed_arcs, decision_arcs)
        add_edge!(graph, arc...)
    end

    n = nv(graph)
    risk = zeros(Int, n, n)
    cost = zeros(Int, n, n)
    weight = zeros(Int, n, n)

    cost[1, 2] = 4
    cost[2, 4] = 5
    cost[1, 3] = 2
    cost[3, 4] = 2

    cost[5, 6] = 3
    cost[6, 8] = 4
    cost[5, 7] = 1
    cost[7, 8] = 2

    risk[1, 3] = -5
    risk[3, 4] = -5
    risk[5, 7] = -4
    risk[7, 8] = -4

    users = User[
        User("U1", 1, 4, risk, cost, weight, nothing),
        User("U2", 5, 8, risk, cost, weight, nothing),
    ]
    edge_price = Dict(
        (1, 3) => 3,
        (3, 4) => 3,
        (5, 7) => 2,
        (7, 8) => 3,
    )
    return HNDPwC(graph, users, decision_arcs, edge_price, nothing)
end

"""
    build_toy_hndp_two_users_all_decision_arcs()

Variant of the two-user HNDP test instance without any fixed-arc competitor
network. All arcs are decision arcs, so the path-based model must fall back to
the `n-1` heuristic bound before enumerating feasible paths.
"""
function build_toy_hndp_two_users_all_decision_arcs()
    graph = DiGraph(8)
    all_decision_arcs = [(1, 2), (2, 4), (5, 6), (6, 8), (1, 3), (3, 4), (5, 7), (7, 8)]
    for arc in all_decision_arcs
        add_edge!(graph, arc...)
    end

    n = nv(graph)
    risk = zeros(Int, n, n)
    cost = zeros(Int, n, n)
    weight = zeros(Int, n, n)

    cost[1, 2] = 4
    cost[2, 4] = 5
    cost[1, 3] = 2
    cost[3, 4] = 2

    cost[5, 6] = 3
    cost[6, 8] = 4
    cost[5, 7] = 1
    cost[7, 8] = 2

    risk[1, 3] = -5
    risk[3, 4] = -5
    risk[5, 7] = -4
    risk[7, 8] = -4

    users = User[
        User("U1", 1, 4, risk, cost, weight, nothing),
        User("U2", 5, 8, risk, cost, weight, nothing),
    ]
    edge_price = Dict(
        (1, 2) => 6,
        (2, 4) => 6,
        (5, 6) => 5,
        (6, 8) => 5,
        (1, 3) => 3,
        (3, 4) => 3,
        (5, 7) => 2,
        (7, 8) => 3,
    )
    return HNDPwC(graph, users, all_decision_arcs, edge_price, nothing)
end

"""
    build_toy_hndp_weighted_user()

Build a small weighted HNDP instance with one user. Besides the feasible fixed
and decision paths, the graph also contains a cheaper but overweight path so we
exercise the constrained-shortest-path logic in both the MIP and AStar-based
subsolvers.
"""
function build_toy_hndp_weighted_user()
    graph = DiGraph(5)
    fixed_arcs = [(1, 2), (2, 4)]
    decision_arcs = [(1, 3), (3, 4), (1, 5), (5, 4)]
    for arc in vcat(fixed_arcs, decision_arcs)
        add_edge!(graph, arc...)
    end

    n = nv(graph)
    risk = zeros(Int, n, n)
    cost = zeros(Int, n, n)
    weight = zeros(Int, n, n)

    cost[1, 2] = 4
    cost[2, 4] = 4
    weight[1, 2] = 2
    weight[2, 4] = 2

    cost[1, 3] = 2
    cost[3, 4] = 2
    weight[1, 3] = 2
    weight[3, 4] = 2
    risk[1, 3] = -4
    risk[3, 4] = -4

    cost[1, 5] = 1
    cost[5, 4] = 1
    weight[1, 5] = 4
    weight[5, 4] = 4
    risk[1, 5] = -6
    risk[5, 4] = -6

    users = User[
        User("UW", 1, 4, risk, cost, weight, 5),
    ]
    edge_price = Dict(
        (1, 3) => 2,
        (3, 4) => 2,
        (1, 5) => 2,
        (5, 4) => 2,
    )
    return HNDPwC(graph, users, decision_arcs, edge_price, nothing)
end

function _run_hndp_model_pair_test(big_m_mode::Symbol)
    hndp = build_toy_hndp_two_users()
    solver = GurobiSolver(Gurobi.Env())

    blc_instance = build_hndp_blc_instance(hndp, solver; big_m_mode=big_m_mode)
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=big_m_mode)

    blc_dir = mktempdir()
    sd_dir = mktempdir()

    blc_param = BLCparam(solver, false, blc_dir, "lp", 60)
    sd_param = MIPparam(solver, false, sd_dir, "lp", 60)

    blc_stats = solve_instance!(blc_instance, blc_param)
    sd_stats = solve_instance!(sd_instance, sd_param)

    @test haskey(blc_stats.data, "Opt")
    @test haskey(sd_stats.data, "Opt")
    @test blc_stats.data["Opt"] ≈ sd_stats.data["Opt"] atol = 1e-6
    @test blc_stats.data["Opt"] ≈ -7 atol = 1e-6
end

function _run_hndp_blc_astar_pair_test(hndp::HNDPwC, big_m_mode::Symbol, expected_opt::Union{Nothing,Float64}=nothing)
    solver = GurobiSolver(Gurobi.Env())

    blc_mip_instance = build_hndp_blc_instance(
        hndp,
        solver;
        big_m_mode=big_m_mode,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )
    blc_astar_instance = build_hndp_blc_instance(
        hndp,
        solver;
        big_m_mode=big_m_mode,
        subproblem_method=HNDP_SUBPROBLEM_ASTAR,
    )

    blc_mip_stats = solve_instance!(blc_mip_instance, BLCparam(solver, false, mktempdir(), "lp", 60))
    blc_astar_stats = solve_instance!(blc_astar_instance, BLCparam(solver, false, mktempdir(), "lp", 60))

    @test haskey(blc_mip_stats.data, "Opt")
    @test haskey(blc_astar_stats.data, "Opt")
    @test blc_mip_stats.data["Opt"] ≈ blc_astar_stats.data["Opt"] atol = 1e-6
    if !isnothing(expected_opt)
        @test blc_astar_stats.data["Opt"] ≈ expected_opt atol = 1e-6
    end
end

function _run_hndp_path_model_test()
    hndp = build_toy_hndp_two_users()
    solver = GurobiSolver(Gurobi.Env())

    blc_instance = build_hndp_blc_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH)
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH)
    path_instance, enum_runtime, path_counts = build_hndp_path_instance(hndp, solver; enumeration_time_limit=60.0)

    blc_param = BLCparam(solver, false, mktempdir(), "lp", 60)
    sd_param = MIPparam(solver, false, mktempdir(), "lp", 60)
    path_param = MIPparam(solver, false, mktempdir(), "lp", 60)

    blc_stats = solve_instance!(blc_instance, blc_param)
    sd_stats = solve_instance!(sd_instance, sd_param)
    path_stats = solve_instance!(path_instance, path_param)

    @test enum_runtime >= 0
    @test all(v >= 1 for v in values(path_counts))
    @test blc_stats.data["Opt"] ≈ sd_stats.data["Opt"] atol = 1e-6
    @test blc_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    @test path_stats.data["Opt"] ≈ -7 atol = 1e-6
end

function _run_hndp_path_model_all_decision_arcs_test()
    hndp = build_toy_hndp_two_users_all_decision_arcs()
    solver = GurobiSolver(Gurobi.Env())

    blc_instance = build_hndp_blc_instance(hndp, solver; big_m_mode=HNDP_BIGM_N_MINUS_ONE)
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=HNDP_BIGM_N_MINUS_ONE)
    path_instance, enum_runtime, path_counts = build_hndp_path_instance(hndp, solver; enumeration_time_limit=60.0)

    blc_param = BLCparam(solver, false, mktempdir(), "lp", 60)
    sd_param = MIPparam(solver, false, mktempdir(), "lp", 60)
    path_param = MIPparam(solver, false, mktempdir(), "lp", 60)

    blc_stats = solve_instance!(blc_instance, blc_param)
    sd_stats = solve_instance!(sd_instance, sd_param)
    path_stats = solve_instance!(path_instance, path_param)

    @test enum_runtime >= 0
    @test all(v >= 1 for v in values(path_counts))
    @test haskey(blc_stats.data, "Opt")
    @test haskey(sd_stats.data, "Opt")
    @test haskey(path_stats.data, "Opt")
    @test blc_stats.data["Opt"] ≈ sd_stats.data["Opt"] atol = 1e-6
    @test blc_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    @test sd_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    # With no fixed-arc network available, all three formulations fall back to
    # the same all-decision-arc design problem.
    @test path_stats.data["Opt"] ≈ 15 atol = 1e-6
end

function _run_hndp_path_model_weighted_test()
    hndp = build_toy_hndp_weighted_user()
    solver = GurobiSolver(Gurobi.Env())

    blc_instance = build_hndp_blc_instance(
        hndp,
        solver;
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )
    path_instance, enum_runtime, path_counts = build_hndp_path_instance(hndp, solver; enumeration_time_limit=60.0)

    blc_stats = solve_instance!(blc_instance, BLCparam(solver, false, mktempdir(), "lp", 60))
    path_stats = solve_instance!(path_instance, MIPparam(solver, false, mktempdir(), "lp", 60))

    @test enum_runtime >= 0
    @test path_counts == Dict("UW" => 2)
    @test haskey(blc_stats.data, "Opt")
    @test haskey(path_stats.data, "Opt")
    @test blc_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    @test path_stats.data["Opt"] ≈ 0 atol = 1e-6
end

@testset "HNDP Model Generation V2 Tests" begin
    _run_hndp_model_pair_test(HNDP_BIGM_FIXED_NETWORK_PATH)
    _run_hndp_model_pair_test(HNDP_BIGM_N_MINUS_ONE)
    _run_hndp_path_model_test()
    _run_hndp_path_model_all_decision_arcs_test()
    _run_hndp_path_model_weighted_test()
    _run_hndp_blc_astar_pair_test(build_toy_hndp_two_users(), HNDP_BIGM_FIXED_NETWORK_PATH, -7.0)
    _run_hndp_blc_astar_pair_test(build_toy_hndp_two_users_all_decision_arcs(), HNDP_BIGM_N_MINUS_ONE, 15.0)
    _run_hndp_blc_astar_pair_test(build_toy_hndp_weighted_user(), HNDP_BIGM_FIXED_NETWORK_PATH)
end
