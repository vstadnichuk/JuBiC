# Reproduce the positive-risk alpha=0.2 and unrestricted benchmark rows from
# the committed instances in this folder. No random instances are generated.
using JuBiC
using JSON

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "examples", "HNDP", "hndp_experiment_runner.jl"))

const INSTANCE_ROOT = joinpath(@__DIR__, "instances")
const OUTPUT_ROOT = joinpath(@__DIR__, "reproduction_output")

function local_instance_config()
    paths = sort(joinpath(path, "instance.json") for path in
                 readdir(INSTANCE_ROOT; join=true)
                 if isfile(joinpath(path, "instance.json")))
    isempty(paths) && error("No saved HNDP instances found in $INSTANCE_ROOT")
    path = joinpath(OUTPUT_ROOT, "instances.json")
    mkpath(dirname(path))
    write(path, JSON.json(Dict("saved_instances" => paths)))
    return path
end

const MODE_SPECS = [
    ("master15_oneworker_sub1_1h", 15, 1, 1,
     ["BlD", "ND_WS_PD", "ND_WS_PD_ASTAR"]),
    ("master1_oneworker_sub15_1h", 1, 15, 1,
     ["BlD", "ND_WS_PD"]),
    ("master1_workers15_sub1_1h", 1, 1, 15,
     ["BlD", "ND_WS_PD", "ND_WS_PD_ASTAR"]),
]

function model_spec(name)
    name == "BlD" && return Dict("name"=>name, "model_type"=>"blc",
        "subproblem_method"=>"mip", "big_m_mode"=>"n_minus_one_most_expensive")
    name == "ND_WS_PD_ASTAR" && return Dict("name"=>name, "model_type"=>"gbc",
        "partial_decomposition"=>true, "subproblem_method"=>"astar")
    return Dict("name"=>name, "model_type"=>"gbc",
        "partial_decomposition"=>true, "subproblem_method"=>"mip")
end

function run_mode!(mode)
    name, master, sub, workers, model_names = mode
    models_path = joinpath(OUTPUT_ROOT, "$(name)_models.json")
    params_path = joinpath(OUTPUT_ROOT, "$(name)_params.json")
    write(models_path, JSON.json(Dict("models" => [model_spec(m) for m in model_names])))
    params = Dict{String,Any}(
        "name"=>name, "mip_solver"=>"Gurobi", "runtime"=>3600.0, "seed"=>42,
        "pareto"=>"OPT", "warmstart"=>true, "bigMwithLC"=>false,
        "trim_coeff"=>false, "integer_obj"=>true, "enable_output_logs"=>false,
        "threads_master"=>master, "threads_sub_con"=>sub,
        "parallel_subsolvers"=>workers, "parallel_separation"=>true,
    )
    write(params_path, JSON.json(Dict("write_run_logs"=>false, "params"=>[params])))
    println("Running $(name): master=$(master), workers=$(workers), subsolver=$(sub)")
    run_hndp_experiments!(local_instance_config(), models_path, params_path;
        output_root=OUTPUT_ROOT, resume=true)
end

selected = isempty(ARGS) ? "all" : ARGS[1]
selected == "all" && foreach(run_mode!, MODE_SPECS)
selected != "all" && begin
    mode = findfirst(m -> m[1] == selected, MODE_SPECS)
    isnothing(mode) && error("Unknown mode $selected. Use all or one of: $(join(first.(MODE_SPECS), ", "))")
    run_mode!(MODE_SPECS[mode])
end
