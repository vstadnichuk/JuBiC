using JSON
using Dates

ENV["GUROBI_HOME"] = get(ENV, "GUROBI_HOME", "C:\\gurobi1301\\win64")
ENV["GUROBI_JL_USE_GUROBI_JLL"] = "false"
include(joinpath(@__DIR__, "hndp_experiment_runner.jl"))

const output_root = get(ENV, "JUBIC_REPRO_OUTPUT", joinpath(pwd(), "benchmark_run", "runs", "fixed_order_big_runs_reproduced_" * Dates.format(now(), "yyyymmdd_HHMMSS")))

function write_json(path, value)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, JSON.json(value, 2))
    end
end

instance_config = Dict{String,Any}(
    "parameter_seeds" => collect(2:4),
    "instances" => [Dict{String,Any}(
        "name" => "sioux_fixed_order_big_runs",
        "instance_type" => "competition",
        "topologies" => ["layered_sioux_falls"],
        "nusers" => [100, 125, 150, 175, 200],
        "length_constrained" => [true, false],
        "alpha" => [0.01, 0.05, 0.1],
        "beta" => [1.0],
        "user_parameter_mode" => ["shared"],
        "construction_cost_min" => 0,
        "construction_cost_max" => 0,
        "availability_budget_fraction" => [1.0],
        "od_pair_mode" => "sampled",
    )],
)

base_models = [
    Dict{String,Any}("name" => "blc_n_minus_one_fixed_order", "model_type" => "blc", "subproblem_method" => "mip", "big_m_mode" => "n_minus_one_most_expensive", "write_solver_instance_files" => false),
    Dict{String,Any}("name" => "blc_fixed_path_fixed_order", "model_type" => "blc", "subproblem_method" => "mip", "big_m_mode" => "fixed_network_path", "write_solver_instance_files" => false),
    Dict{String,Any}("name" => "blc_fixed_path_current_cost_fixed_order", "model_type" => "blc", "subproblem_method" => "mip", "big_m_mode" => "fixed_network_path_current_cost", "write_solver_instance_files" => false),
]
blclag_models = [Dict{String,Any}("name" => "blclag_n_minus_one_fixed_order", "model_type" => "blclag", "subproblem_method" => "blc_jump", "big_m_mode" => "n_minus_one_most_expensive", "write_solver_instance_files" => false)]

base_param = Dict{String,Any}(
    "mip_solver" => "Gurobi",
    "runtime" => 600.0,
    "threads_master" => 8,
    "threads_sub_con" => 8,
    "parallel_separation" => false,
    "seed" => 42,
    "branching_rule" => "fixed_linking_order",
    "print_collected_cuts" => false,
    "print_big_m_summary" => true,
)
blc_params = [merge(copy(base_param), Dict{String,Any}("name" => "gurobi_10min_fixed_order_seq8", "warmstart" => true))]
blclag_params = [
    merge(copy(base_param), Dict{String,Any}("name" => "gurobi_10min_fixed_order_seq8_warm", "warmstart" => true)),
    merge(copy(base_param), Dict{String,Any}("name" => "gurobi_10min_fixed_order_seq8_cold", "warmstart" => false)),
]

function run_suite!(name, models, params)
    suite_root = joinpath(output_root, name)
    input_root = joinpath(suite_root, "input_configs")
    instance_path = joinpath(input_root, "instances.json")
    model_path = joinpath(input_root, "models.json")
    param_path = joinpath(input_root, "params.json")
    write_json(instance_path, instance_config)
    write_json(model_path, Dict{String,Any}("models" => models))
    write_json(param_path, Dict{String,Any}("execution_mode" => "subprocess_per_experiment", "write_run_logs" => true, "params" => params))
    println("START_SUITE=$(name) OUTPUT_ROOT=$(suite_root)")
    println("EXPECTED_SUITE_JOBS=$(60 * length(models) * length(params))")
    run_hndp_experiments!(instance_path, model_path, param_path; output_root=suite_root, resume=true)
    println("DONE_SUITE=$(name) OUTPUT_ROOT=$(suite_root)")
end

mkpath(output_root)
println("START_OUTPUT_ROOT=$(output_root)")
run_suite!("blc", base_models, blc_params)
run_suite!("blclag", blclag_models, blclag_params)
println("DONE_OUTPUT_ROOT=$(output_root)")
