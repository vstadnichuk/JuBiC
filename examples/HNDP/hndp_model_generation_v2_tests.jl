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

@testset "HNDP Model Generation V2 Tests" begin
    _run_hndp_model_pair_test(HNDP_BIGM_FIXED_NETWORK_PATH)
    _run_hndp_model_pair_test(HNDP_BIGM_N_MINUS_ONE)
end
