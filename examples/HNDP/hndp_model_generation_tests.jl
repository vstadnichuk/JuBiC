using Test
using JuBiC
using JuMP
using Gurobi
using Graphs

include("hndp_model_generation.jl")

const HNDP_TEST_TIME_LIMIT = 60
const TOY_TWO_USER_OPT = -7.0
const TOY_ALL_DECISION_OPT = -7.0
const TOY_WEIGHTED_OPT = -4.0
# These are the optima after the GBC regression tests convert the toy risks
# to nonnegative values. They were checked by enumerating every feasible
# leader decision pattern and each follower path.
const TOY_TWO_USER_GBC_OPT = 0.0
const TOY_ALL_DECISION_GBC_OPT = 22.0
const TOY_WEIGHTED_GBC_OPT = 0.0

_hndp_test_tempdir() = JuBiC.repo_local_tempdir("tests", "hndp_model_generation"; prefix="hndp_test")

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

function _copy_hndp_with_nonnegative_risk(hndp::HNDPwC)
    users = User[]
    for user in hndp.users
        risk = copy(user.mrisk)
        for i in axes(risk, 1), j in axes(risk, 2)
            risk[i, j] = abs(risk[i, j])
        end
        push!(users, User(user.uname, user.origin, user.destination, risk, user.mcost, user.mweight, user.weighlimit))
    end
    return HNDPwC(hndp.mygraph, users, hndp.edgeA, hndp.edge_price, hndp.minweights)
end

function _run_hndp_model_pair_test(big_m_mode::Symbol)
    hndp = build_toy_hndp_two_users()
    solver = GurobiSolver()

    blc_instance = build_hndp_blc_instance(hndp, solver; big_m_mode=big_m_mode)
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=big_m_mode)

    blc_dir = _hndp_test_tempdir()
    sd_dir = _hndp_test_tempdir()

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
    solver = GurobiSolver()

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

    blc_mip_stats = solve_instance!(blc_mip_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    blc_astar_stats = solve_instance!(blc_astar_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test haskey(blc_mip_stats.data, "Opt")
    @test haskey(blc_astar_stats.data, "Opt")
    @test blc_mip_stats.data["Opt"] ≈ blc_astar_stats.data["Opt"] atol = 1e-6
    if !isnothing(expected_opt)
        @test blc_astar_stats.data["Opt"] ≈ expected_opt atol = 1e-6
    end
end

function _run_hndp_path_model_test()
    hndp = build_toy_hndp_two_users()
    solver = GurobiSolver()

    blc_instance = build_hndp_blc_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH)
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH)
    path_instance, enum_runtime, path_counts = build_hndp_path_instance(hndp, solver; enumeration_time_limit=60.0)

    blc_param = BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60)
    sd_param = MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60)
    path_param = MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60)

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
    solver = GurobiSolver()

    blc_instance = build_hndp_blc_instance(hndp, solver; big_m_mode=HNDP_BIGM_N_MINUS_ONE)
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=HNDP_BIGM_N_MINUS_ONE)
    path_instance, enum_runtime, path_counts = build_hndp_path_instance(hndp, solver; enumeration_time_limit=60.0)

    blc_param = BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60)
    sd_param = MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60)
    path_param = MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60)

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
    @test path_stats.data["Opt"] ≈ TOY_ALL_DECISION_OPT atol = 1e-6
end

function _run_hndp_path_model_weighted_test()
    hndp = build_toy_hndp_weighted_user()
    solver = GurobiSolver()

    blc_instance = build_hndp_blc_instance(
        hndp,
        solver;
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )
    path_instance, enum_runtime, path_counts = build_hndp_path_instance(hndp, solver; enumeration_time_limit=60.0)

    blc_stats = solve_instance!(blc_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    path_stats = solve_instance!(path_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test enum_runtime >= 0
    @test path_counts == Dict("UW" => 2)
    @test haskey(blc_stats.data, "Opt")
    @test haskey(path_stats.data, "Opt")
    @test blc_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    @test path_stats.data["Opt"] ≈ TOY_WEIGHTED_OPT atol = 1e-6
end

function _run_hndp_path_parallel_equivalence_test(hndp::HNDPwC, expected_opt::Float64)
    solver = GurobiSolver()

    seq_instance, seq_runtime, seq_counts = build_hndp_path_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        parallelize=false,
    )
    par_instance, par_runtime, par_counts = build_hndp_path_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        parallelize=true,
    )

    seq_stats = solve_instance!(seq_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    par_stats = solve_instance!(par_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test seq_runtime >= 0
    @test par_runtime >= 0
    @test seq_counts == par_counts
    @test haskey(seq_stats.data, "Opt")
    @test haskey(par_stats.data, "Opt")
    @test seq_stats.data["Opt"] ≈ par_stats.data["Opt"] atol = 1e-6
    @test par_stats.data["Opt"] ≈ expected_opt atol = 1e-6
end

function _run_hndp_hybrid_path_equivalence_test(hndp::HNDPwC, expected_opt::Float64)
    solver = GurobiSolver()

    path_instance, path_runtime, path_counts = build_hndp_path_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        parallelize=true,
    )
    hybrid_instance, hybrid_runtime, hybrid_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        fallback_mode=HNDP_HYBRID_FALLBACK_SD,
    )

    path_stats = solve_instance!(path_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    hybrid_stats = solve_instance!(hybrid_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test path_runtime >= 0
    @test hybrid_runtime >= 0
    @test hybrid_counts == path_counts
    @test isempty(fallback_users)
    @test isempty(path_count_fallback_users)
    @test haskey(path_stats.data, "Opt")
    @test haskey(hybrid_stats.data, "Opt")
    @test hybrid_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    @test hybrid_stats.data["Opt"] ≈ expected_opt atol = 1e-6
end

function _run_hndp_hybrid_all_fallback_test(hndp::HNDPwC, big_m_mode::Symbol, expected_opt::Float64)
    solver = GurobiSolver()

    hybrid_instance, hybrid_runtime, hybrid_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_instance(
        hndp,
        solver;
        enumeration_time_limit=0.0,
        fallback_mode=HNDP_HYBRID_FALLBACK_SD,
        big_m_mode=big_m_mode,
    )
    sd_instance = build_hndp_sd_instance(hndp, solver; big_m_mode=big_m_mode)

    hybrid_stats = solve_instance!(hybrid_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    sd_stats = solve_instance!(sd_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test hybrid_runtime ≈ 0.0 atol = 1e-6
    @test isempty(hybrid_counts)
    @test length(fallback_users) == length(hndp.users)
    @test isempty(path_count_fallback_users)
    @test Set(fallback_users) == Set(string(user.uname) for user in hndp.users)
    @test haskey(hybrid_stats.data, "Opt")
    @test haskey(sd_stats.data, "Opt")
    @test hybrid_stats.data["Opt"] ≈ sd_stats.data["Opt"] atol = 1e-6
    @test hybrid_stats.data["Opt"] ≈ expected_opt atol = 1e-6
end

function _run_hndp_hybrid_blc_path_equivalence_test(hndp::HNDPwC, expected_opt::Float64)
    solver = GurobiSolver()

    path_instance, path_runtime, path_counts = build_hndp_path_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        parallelize=true,
    )
    hybrid_blc_instance, hybrid_runtime, hybrid_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_blc_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )

    path_stats = solve_instance!(path_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    hybrid_stats = solve_instance!(hybrid_blc_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test path_runtime >= 0
    @test hybrid_runtime >= 0
    @test hybrid_counts == path_counts
    @test isempty(fallback_users)
    @test isempty(path_count_fallback_users)
    @test isempty(hybrid_blc_instance.subproblems)
    @test haskey(path_stats.data, "Opt")
    @test haskey(hybrid_stats.data, "Opt")
    @test hybrid_stats.data["Opt"] ≈ path_stats.data["Opt"] atol = 1e-6
    @test hybrid_stats.data["Opt"] ≈ expected_opt atol = 1e-6
end

function _run_hndp_hybrid_blc_all_fallback_test(hndp::HNDPwC, big_m_mode::Symbol, expected_opt::Float64)
    solver = GurobiSolver()

    hybrid_blc_instance, hybrid_runtime, hybrid_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_blc_instance(
        hndp,
        solver;
        enumeration_time_limit=0.0,
        big_m_mode=big_m_mode,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )
    blc_instance = build_hndp_blc_instance(
        hndp,
        solver;
        big_m_mode=big_m_mode,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )

    hybrid_stats = solve_instance!(hybrid_blc_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    blc_stats = solve_instance!(blc_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test hybrid_runtime ≈ 0.0 atol = 1e-6
    @test isempty(hybrid_counts)
    @test length(fallback_users) == length(hndp.users)
    @test isempty(path_count_fallback_users)
    @test Set(fallback_users) == Set(string(user.uname) for user in hndp.users)
    @test length(hybrid_blc_instance.subproblems) == length(hndp.users)
    @test haskey(hybrid_stats.data, "Opt")
    @test haskey(blc_stats.data, "Opt")
    @test hybrid_stats.data["Opt"] ≈ blc_stats.data["Opt"] atol = 1e-6
    @test hybrid_stats.data["Opt"] ≈ expected_opt atol = 1e-6
end

function _run_hndp_hybrid_path_count_fallback_test()
    solver = GurobiSolver()
    hndp = build_toy_hndp_many_paths_user()
    all_arc_count = length(_hndp_all_arcs(hndp))

    hybrid_instance, hybrid_runtime, hybrid_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        fallback_mode=HNDP_HYBRID_FALLBACK_SD,
        fallback_if_paths_exceed_flow_vars=true,
    )
    sd_instance = build_hndp_sd_instance(hndp, solver)

    hybrid_blc_instance, hybrid_blc_runtime, hybrid_blc_counts, hybrid_blc_fallback_users, hybrid_blc_path_count_fallback_users = build_hndp_hybrid_blc_instance(
        hndp,
        solver;
        enumeration_time_limit=60.0,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        fallback_if_paths_exceed_flow_vars=true,
    )
    blc_instance = build_hndp_blc_instance(hndp, solver; subproblem_method=HNDP_SUBPROBLEM_MIP)

    @test hybrid_runtime >= 0
    @test hybrid_blc_runtime >= 0
    @test hybrid_counts["UM"] > all_arc_count
    @test hybrid_blc_counts["UM"] > all_arc_count
    @test fallback_users == ["UM"]
    @test hybrid_blc_fallback_users == ["UM"]
    @test path_count_fallback_users == ["UM"]
    @test hybrid_blc_path_count_fallback_users == ["UM"]

    hybrid_stats = solve_instance!(hybrid_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    sd_stats = solve_instance!(sd_instance, MIPparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    hybrid_blc_stats = solve_instance!(hybrid_blc_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))
    blc_stats = solve_instance!(blc_instance, BLCparam(solver, false, _hndp_test_tempdir(), "lp", 60))

    @test isapprox(hybrid_stats.data["Opt"], sd_stats.data["Opt"]; atol=1e-6)
    @test isapprox(hybrid_blc_stats.data["Opt"], blc_stats.data["Opt"]; atol=1e-6)
end

function _run_hndp_gbc_subsolver_triplet_test(hndp::HNDPwC, big_m_mode::Symbol, expected_opt::Float64)
    solver = GurobiSolver()

    gbc_mip = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        big_m_mode=big_m_mode,
    )
    gbc_blc = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_BLC_JUMP,
        big_m_mode=big_m_mode,
    )
    gbc_astar = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_ASTAR,
        big_m_mode=big_m_mode,
    )

    stats_mip = solve_instance!(gbc_mip, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, 60, true))
    stats_blc = solve_instance!(gbc_blc, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, 60, true))
    stats_astar = solve_instance!(gbc_astar, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, 60, true))

    @test haskey(stats_mip.data, "Opt")
    @test haskey(stats_blc.data, "Opt")
    @test haskey(stats_astar.data, "Opt")
    @test stats_mip.data["Opt"] ≈ stats_blc.data["Opt"] atol = 1e-6
    @test stats_mip.data["Opt"] ≈ stats_astar.data["Opt"] atol = 1e-6
    @test stats_astar.data["Opt"] ≈ expected_opt atol = 1e-6
end

function _run_hndp_gbc_weighted_subsolver_pair_test(hndp::HNDPwC, expected_opt::Float64)
    solver = GurobiSolver()

    gbc_mip = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
    )
    gbc_astar = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_ASTAR,
    )

    stats_mip = solve_instance!(gbc_mip, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, 60, true))
    stats_astar = solve_instance!(gbc_astar, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, 60, true))

    @test haskey(stats_mip.data, "Opt")
    @test haskey(stats_astar.data, "Opt")
    @test stats_mip.data["Opt"] ≈ stats_astar.data["Opt"] atol = 1e-6
    @test stats_astar.data["Opt"] ≈ expected_opt atol = 1e-6
end




"""
    _run_hndp_mibs_toy_solve_test()

Validate the single-follower MiBS transformation on the original toy HNDP
instance against the corresponding GBC model.
"""
function _run_hndp_mibs_toy_solve_test()
    solver = GurobiSolver()
    hndp = build_toy_HNDPwC()

    gbc_instance = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        big_m_mode=HNDP_BIGM_N_MINUS_ONE,
    )
    mibs_instance = build_hndp_mibs_instance(hndp, solver)

    gbc_stats = solve_instance!(gbc_instance, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, HNDP_TEST_TIME_LIMIT, true))
    mibs_stats = solve_instance!(mibs_instance, MibSparam(false, _hndp_test_tempdir(), HNDP_TEST_TIME_LIMIT))

    @test haskey(gbc_stats.data, "Opt")
    @test haskey(mibs_stats.data, "Opt")
    @test gbc_stats.data["Opt"] ≈ mibs_stats.data["Opt"] atol = 1e-6
end

"""
    _run_hndp_mibs_multi_follower_toy_solve_test()

Validate the multi-follower MiBS path on the two-user toy HNDP instance. This
covers the `merge_subproblems` transformation before the model is exported to
MiBS and checks that the resulting single-follower reformulation matches GBC.
"""
function _run_hndp_mibs_multi_follower_toy_solve_test()
    solver = GurobiSolver()
    hndp = build_toy_hndp_two_users()

    gbc_instance = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    )
    mibs_instance = build_hndp_mibs_instance(hndp, solver)

    gbc_stats = solve_instance!(gbc_instance, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, HNDP_TEST_TIME_LIMIT, true))
    mibs_stats = solve_instance!(mibs_instance, MibSparam(false, _hndp_test_tempdir(), HNDP_TEST_TIME_LIMIT))

    @test haskey(gbc_stats.data, "Opt")
    @test haskey(mibs_stats.data, "Opt")
    @test gbc_stats.data["Opt"] ≈ mibs_stats.data["Opt"] atol = 1e-6
    @test gbc_stats.data["Opt"] ≈ TOY_TWO_USER_OPT atol = 1e-6
end

"""
    _run_hndp_mibs_runtime_toy_solve_test()

Check that the JuBiC-side MiBS runner stores the requested time limit in the
generated parameter file path and still returns the same objective as GBC.
"""
function _run_hndp_mibs_runtime_toy_solve_test()
    solver = GurobiSolver()
    hndp = build_toy_hndp_two_users()

    gbc_instance = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    )
    mibs_instance = build_hndp_mibs_instance(hndp, solver)

    gbc_stats = solve_instance!(gbc_instance, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, HNDP_TEST_TIME_LIMIT, true))
    mibs_stats = solve_instance!(mibs_instance, MibSparam(false, _hndp_test_tempdir(), HNDP_TEST_TIME_LIMIT))

    @test haskey(gbc_stats.data, "Opt")
    @test haskey(mibs_stats.data, "Opt")
    @test get(mibs_stats.data, "MibSExecution", nothing) == "JuBiCParamFileRunner"
    @test get(mibs_stats.data, "time_limit", nothing) == HNDP_TEST_TIME_LIMIT
    @test gbc_stats.data["Opt"] ≈ mibs_stats.data["Opt"] atol = 1e-6
end

function _run_hndp_mibs_runtime_toy_solve_test()
    solver = GurobiSolver()
    hndp = build_toy_hndp_two_users()
    mibs_outdir = _hndp_test_tempdir()

    gbc_instance = build_hndp_gbc_instance(
        hndp,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH,
    )
    mibs_instance = build_hndp_mibs_instance(hndp, solver)

    gbc_stats = solve_instance!(gbc_instance, GBCparam(solver, false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, HNDP_TEST_TIME_LIMIT, true))
    mibs_stats = solve_instance!(mibs_instance, MibSparam(false, mibs_outdir, HNDP_TEST_TIME_LIMIT))

    @test haskey(gbc_stats.data, "Opt")
    @test haskey(mibs_stats.data, "Opt")
    @test get(mibs_stats.data, "MibSExecution", nothing) == "JuBiCParamFileRunner"
    @test get(mibs_stats.data, "time_limit", nothing) == HNDP_TEST_TIME_LIMIT
    @test isfile(joinpath(mibs_outdir, "mibs.par"))
    @test occursin("Alps_timeLimit $(Float64(HNDP_TEST_TIME_LIMIT))", read(joinpath(mibs_outdir, "mibs.par"), String))
    @test isapprox(gbc_stats.data["Opt"], mibs_stats.data["Opt"]; atol=1e-6)
end

function _run_hndp_astar_negative_master_cost_guard_test()
    hndp = build_toy_hndp_two_users()
    subsolver = build_hndp_astar_user(hndp.users[1], hndp, hndp.edgeA)
    params = GBCparam(GurobiSolver(), false, _hndp_test_tempdir(), "lp", PARETO_OPTIMALITY_ONLY, 60, true)

    @test_throws ArgumentError JuBiC.compute_lower_bound_master_contribution(subsolver, params, HNDP_TEST_TIME_LIMIT)
end

function build_toy_hndp_many_paths_user()
    graph = DiGraph(6)
    all_arcs = Tuple{Int,Int}[]
    for i in 1:5
        for j in (i + 1):6
            add_edge!(graph, i, j)
            push!(all_arcs, (i, j))
        end
    end

    n = nv(graph)
    risk = zeros(Int, n, n)
    cost = zeros(Int, n, n)
    weight = zeros(Int, n, n)
    for (i, j) in all_arcs
        cost[i, j] = 1
        risk[i, j] = -1
    end

    users = User[
        User("UM", 1, 6, risk, cost, weight, nothing),
    ]
    edge_price = Dict(arc => 1 for arc in all_arcs)
    return HNDPwC(graph, users, all_arcs, edge_price, nothing)
end


@testset "HNDP Model Generation Tests" begin
    @testset "Arc-Based Models" begin
        _run_hndp_model_pair_test(HNDP_BIGM_FIXED_NETWORK_PATH)
        _run_hndp_model_pair_test(HNDP_BIGM_N_MINUS_ONE)
        _run_hndp_blc_astar_pair_test(build_toy_hndp_two_users(), HNDP_BIGM_FIXED_NETWORK_PATH, TOY_TWO_USER_OPT)
        _run_hndp_blc_astar_pair_test(build_toy_hndp_two_users_all_decision_arcs(), HNDP_BIGM_N_MINUS_ONE, TOY_ALL_DECISION_OPT)
        _run_hndp_blc_astar_pair_test(build_toy_hndp_weighted_user(), HNDP_BIGM_FIXED_NETWORK_PATH)
    end

    @testset "GBC Models" begin
        _run_hndp_gbc_subsolver_triplet_test(_copy_hndp_with_nonnegative_risk(build_toy_hndp_two_users()), HNDP_BIGM_FIXED_NETWORK_PATH, TOY_TWO_USER_GBC_OPT)
        _run_hndp_gbc_subsolver_triplet_test(_copy_hndp_with_nonnegative_risk(build_toy_hndp_two_users_all_decision_arcs()), HNDP_BIGM_N_MINUS_ONE, TOY_ALL_DECISION_GBC_OPT)
        _run_hndp_gbc_weighted_subsolver_pair_test(_copy_hndp_with_nonnegative_risk(build_toy_hndp_weighted_user()), TOY_WEIGHTED_GBC_OPT)
        _run_hndp_astar_negative_master_cost_guard_test()
    end

    @testset "Path Model" begin
        _run_hndp_path_model_test()
        _run_hndp_path_model_all_decision_arcs_test()
        _run_hndp_path_model_weighted_test()
        _run_hndp_path_parallel_equivalence_test(build_toy_hndp_two_users(), TOY_TWO_USER_OPT)
        _run_hndp_path_parallel_equivalence_test(build_toy_hndp_two_users_all_decision_arcs(), TOY_ALL_DECISION_OPT)
        _run_hndp_path_parallel_equivalence_test(build_toy_hndp_weighted_user(), TOY_WEIGHTED_OPT)
    end

    @testset "Hybrid SD Fallback" begin
        _run_hndp_hybrid_path_equivalence_test(build_toy_hndp_two_users(), TOY_TWO_USER_OPT)
        _run_hndp_hybrid_path_equivalence_test(build_toy_hndp_two_users_all_decision_arcs(), TOY_ALL_DECISION_OPT)
        _run_hndp_hybrid_all_fallback_test(build_toy_hndp_two_users(), HNDP_BIGM_FIXED_NETWORK_PATH, TOY_TWO_USER_OPT)
        _run_hndp_hybrid_all_fallback_test(build_toy_hndp_two_users_all_decision_arcs(), HNDP_BIGM_N_MINUS_ONE, TOY_ALL_DECISION_OPT)
        _run_hndp_hybrid_path_count_fallback_test()
    end

    @testset "Hybrid BlC Fallback" begin
        _run_hndp_hybrid_blc_path_equivalence_test(build_toy_hndp_two_users(), TOY_TWO_USER_OPT)
        _run_hndp_hybrid_blc_path_equivalence_test(build_toy_hndp_two_users_all_decision_arcs(), TOY_ALL_DECISION_OPT)
        _run_hndp_hybrid_blc_all_fallback_test(build_toy_hndp_two_users(), HNDP_BIGM_FIXED_NETWORK_PATH, TOY_TWO_USER_OPT)
        _run_hndp_hybrid_blc_all_fallback_test(build_toy_hndp_two_users_all_decision_arcs(), HNDP_BIGM_N_MINUS_ONE, TOY_ALL_DECISION_OPT)
    end

    @testset "MibS toy fixtures" begin
        _run_hndp_mibs_toy_solve_test()
        _run_hndp_mibs_multi_follower_toy_solve_test()
        _run_hndp_mibs_runtime_toy_solve_test()
    end

end
