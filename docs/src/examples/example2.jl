using JuBiC

solver = GurobiSolver()
optimizer = () -> get_next_optimizer(solver)
instance = get_GBC_instance(joinpath(@__DIR__, "example2.mps"),
    joinpath(@__DIR__, "example2.aux"), optimizer)
output = JuBiC.repo_local_tempdir("docs", "mps_aux"; prefix="gbc")
params = GBCparam(solver, false, output, "lp", PARETO_OPTIMALITY_ONLY, 60.0)
params.stats.data["enable_output_logs"] = false
result = solve_instance!(instance, params)
println(result.data["Opt_status"])
