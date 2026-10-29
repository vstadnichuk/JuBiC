using JuBiC
include(joinpath(@__DIR__, "..", "..", "..", "examples", "HNDP", "hndp_model_generation.jl"))

# A tiny in-memory routing fixture keeps this solver illustration offline.
hndp = build_toy_HNDPwC()
solver = GurobiSolver()
instance = build_hndp_gbc_instance(hndp, solver;
    subproblem_method=HNDP_SUBPROBLEM_ASTAR)
output = JuBiC.repo_local_tempdir("docs", "routing_oracle"; prefix="gbc")
params = GBCparam(solver, false, output, "lp", PARETO_OPTIMALITY_ONLY, 60.0)
params.stats.data["enable_output_logs"] = false
result = solve_instance!(instance, params)
println(result.data["Opt_status"])
