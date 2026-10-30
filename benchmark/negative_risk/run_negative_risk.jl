# Reproduce the negative-risk alpha=0.2 and unrestricted benchmark rows from
# the committed instances in this folder. No random instances are generated.
using JuBiC
using JSON
using CSV
using DataFrames

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
    mkpath(dirname(path)); write(path, JSON.json(Dict("saved_instances"=>paths)))
    return path
end

const MODELS = [
    Dict("name"=>"BlD", "model_type"=>"blc", "subproblem_method"=>"mip",
         "big_m_mode"=>"n_minus_one_most_expensive"),
    Dict("name"=>"ND", "model_type"=>"gbc", "partial_decomposition"=>false,
         "subproblem_method"=>"mip", "warmstart"=>false),
    Dict("name"=>"ND_WS", "model_type"=>"gbc", "partial_decomposition"=>false,
         "subproblem_method"=>"mip", "warmstart"=>true),
    Dict("name"=>"ND_WS_PD", "model_type"=>"gbc", "partial_decomposition"=>true,
         "subproblem_method"=>"mip", "warmstart"=>true),
    Dict("name"=>"ND_WS_PD_BS", "model_type"=>"gbc", "partial_decomposition"=>true,
         "subproblem_method"=>"blc_jump", "warmstart"=>true),
]

function run_model!(model)
    name = String(model["name"])
    model_path = joinpath(OUTPUT_ROOT, "$(lowercase(name))_models.json")
    params_path = joinpath(OUTPUT_ROOT, "$(lowercase(name))_params.json")
    write(model_path, JSON.json(Dict("models"=>[model])))
    params = Dict{String,Any}(
        "name"=>"pricing_negative_risk_k25_$(lowercase(name))_1h",
        "mip_solver"=>"Gurobi", "runtime"=>3600.0, "seed"=>42,
        "pareto"=>"OPT", "warmstart"=>true, "bigMwithLC"=>false,
        "trim_coeff"=>false, "integer_obj"=>true, "enable_output_logs"=>false,
        "threads_master"=>(name == "BlD" ? 15 : 1),
        "threads_sub_con"=>1, "parallel_subsolvers"=>(name == "BlD" ? 1 : 15),
        "parallel_separation"=>true,
    )
    write(params_path, JSON.json(Dict("write_run_logs"=>false, "params"=>[params])))
    run_hndp_experiments!(local_instance_config(), model_path, params_path;
        output_root=joinpath(OUTPUT_ROOT, lowercase(name)), resume=true)
end

foreach(run_model!, MODELS)
frames = DataFrame[]
for model in MODELS
    path = joinpath(OUTPUT_ROOT, lowercase(String(model["name"])), "results", "batch_summary.csv")
    isfile(path) && push!(frames, CSV.read(path, DataFrame))
end
if !isempty(frames)
    CSV.write(joinpath(@__DIR__, "results_network_design_negative_reproduced.csv"),
        reduce((a,b) -> vcat(a,b; cols=:union), frames))
end
