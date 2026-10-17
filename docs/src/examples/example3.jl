using JuBiC
include(joinpath(@__DIR__, "..", "..", "..", "examples", "HNDP", "hndp_model_generation.jl"))

config = Dict{String,Any}("parameter_seeds" => [17], "instances" => [
    Dict{String,Any}("name"=>"routing_oracle", "instance_type"=>"constrained_shortest_path",
        "topologies"=>["sioux_falls"], "nusers"=>[2],
        "length_constrained"=>[false])
])
generated = only(generate_hndp_networks(config))
solver = GurobiSolver()
instance = build_hndp_gbc_instance(generated.instance, solver;
    subproblem_method=HNDP_SUBPROBLEM_ASTAR)
output = JuBiC.repo_local_tempdir("docs", "routing_oracle"; prefix="gbc")
params = GBCparam(solver, false, output, "lp", PARETO_OPTIMALITY_ONLY, 60.0)
params.stats.data["enable_output_logs"] = false
result = solve_instance!(instance, params)
println(result.data["Opt_status"])
