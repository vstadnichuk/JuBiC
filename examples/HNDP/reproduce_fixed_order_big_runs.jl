using JSON
using Dates

ENV["GUROBI_HOME"] = get(ENV, "GUROBI_HOME", "C:\\gurobi1301\\win64")
ENV["GUROBI_JL_USE_GUROBI_JLL"] = "false"
include(joinpath(@__DIR__, "hndp_experiment_runner.jl"))

const output_root = get(ENV, "JUBIC_REPRO_OUTPUT", joinpath(pwd(), "benchmark_run", "runs", "fixed_order_big_runs_reproduced_" * Dates.format(now(), "yyyymmdd_HHMMSS")))
const gurobi_threads = parse(Int, get(ENV, "JUBIC_THREADS", "8"))
const sweep_users = parse.(Int, split(get(ENV, "JUBIC_USERS", "100,125,150,175,200"), ","))
const branching_rule = get(ENV, "JUBIC_BRANCHING_RULE", "fixed_linking_order")
const gurobi_minimal = get(ENV, "JUBIC_GUROBI_MINIMAL", "0")
const topology = get(ENV, "JUBIC_TOPOLOGY", "layered_sioux_falls")
const decision_arc_values = parse.(Int, filter(!isempty, split(get(ENV, "JUBIC_DECISION_ARCS", ""), ",")))
const single_layer_modes = filter(!isempty, split(get(ENV, "JUBIC_SINGLE_LAYER_MODES", ""), ","))
const is_single_layer = !isempty(single_layer_modes)
const instance_count = length(sweep_users) * 3 * (3 + 1) * (is_single_layer ? length(decision_arc_values) * length(single_layer_modes) : 1)
const decision_arc_text = isempty(decision_arc_values) ? "none" : join(decision_arc_values, ',')
const single_layer_mode_text = isempty(single_layer_modes) ? "none" : join(single_layer_modes, ',')
const run_blclag = lowercase(get(ENV, "JUBIC_RUN_BLCLAG", "1")) in ("1", "true", "yes")
const run_blc = lowercase(get(ENV, "JUBIC_RUN_BLC", "1")) in ("1", "true", "yes")
const blclag_warm_only = lowercase(get(ENV, "JUBIC_BLCLAG_WARM_ONLY", "0")) in ("1", "true", "yes")

const console_log_path = get(ENV, "JUBIC_CONSOLE_LOG", joinpath(output_root, "console.log"))
mkpath(dirname(console_log_path))
const console_log_io = open(console_log_path, "w")
const original_stdout = stdout
const original_stderr = stderr
const output_pipe = Pipe()
Base.link_pipe!(output_pipe)
redirect_stdout(output_pipe.in)
redirect_stderr(output_pipe.in)
@async begin
    while !eof(output_pipe.out)
        data = readavailable(output_pipe.out)
        if !isempty(data)
            write(original_stdout, data)
            write(console_log_io, data)
            flush(original_stdout)
            flush(console_log_io)
        else
            yield()
        end
    end
end

function write_json(path, value)
    mkpath(dirname(path))
    open(path, "w") do io
        write(io, JSON.json(value, 2))
    end
end

instance_config = Dict{String,Any}(
    "parameter_seeds" => collect(2:4),
    "instances" => [Dict{String,Any}(
        "name" => "fixed_order_big_runs",
        "instance_type" => "competition",
        "topologies" => [topology],
        "nusers" => sweep_users,
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
if is_single_layer
    instance_config["instances"][1]["decision_arc_count"] = decision_arc_values
    instance_config["instances"][1]["single_layer_arc_mode"] = single_layer_modes
end

base_models = [
    Dict{String,Any}("name" => "blc_n_minus_one_fixed_order", "model_type" => "blc", "subproblem_method" => "mip", "big_m_mode" => "n_minus_one_most_expensive", "write_solver_instance_files" => false),
    Dict{String,Any}("name" => "blc_fixed_path_fixed_order", "model_type" => "blc", "subproblem_method" => "mip", "big_m_mode" => "fixed_network_path", "write_solver_instance_files" => false),
    Dict{String,Any}("name" => "blc_fixed_path_current_cost_fixed_order", "model_type" => "blc", "subproblem_method" => "mip", "big_m_mode" => "fixed_network_path_current_cost", "write_solver_instance_files" => false),
]
blclag_models = [Dict{String,Any}("name" => "blclag_n_minus_one_fixed_order", "model_type" => "blclag", "subproblem_method" => "blc_jump", "big_m_mode" => "n_minus_one_most_expensive", "write_solver_instance_files" => false)]

base_param = Dict{String,Any}(
    "mip_solver" => "Gurobi",
    "runtime" => 600.0,
    "threads_master" => gurobi_threads,
    "threads_sub_con" => gurobi_threads,
    "parallel_separation" => false,
    "seed" => 42,
    "branching_rule" => branching_rule,
    "print_collected_cuts" => false,
    "print_big_m_summary" => true,
)
blc_params = [merge(copy(base_param), Dict{String,Any}("name" => "gurobi_10min_fixed_order_seq8", "warmstart" => true))]
blclag_params = [
    merge(copy(base_param), Dict{String,Any}("name" => "gurobi_10min_fixed_order_seq8_warm", "warmstart" => true)),
    merge(copy(base_param), Dict{String,Any}("name" => "gurobi_10min_fixed_order_seq8_cold", "warmstart" => false)),
]
if blclag_warm_only
    blclag_params = blclag_params[1:1]
end

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
    println("EXPECTED_SUITE_JOBS=$(instance_count * length(models) * length(params))")
    run_hndp_experiments!(instance_path, model_path, param_path; output_root=suite_root, resume=true)
    println("DONE_SUITE=$(name) OUTPUT_ROOT=$(suite_root)")
end

mkpath(output_root)
println("START_OUTPUT_ROOT=$(output_root)")
println("GUROBI_THREADS=$(gurobi_threads)")
println("SWEEP_USERS=$(join(sweep_users, ','))")
println("BRANCHING_RULE=$(branching_rule)")
println("GUROBI_MINIMAL=$(gurobi_minimal)")
println("TOPOLOGY=$(topology)")
println("DECISION_ARCS=$(decision_arc_text)")
println("SINGLE_LAYER_MODES=$(single_layer_mode_text)")
if run_blc
    run_suite!("blc", base_models, blc_params)
else
    println("SKIP_SUITE=blc")
end
if run_blclag
    run_suite!("blclag", blclag_models, blclag_params)
else
    println("SKIP_SUITE=blclag")
end
println("DONE_OUTPUT_ROOT=$(output_root)")
