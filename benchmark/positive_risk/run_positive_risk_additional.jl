# Reproduce the supplementary positive-risk alpha=0.2 and unrestricted GBC
# rows from the committed instances in this folder.
using JuBiC
using JSON
using CSV
using DataFrames

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
include(joinpath(REPO_ROOT, "examples", "HNDP", "hndp_experiment_runner.jl"))

const INSTANCE_ROOT = joinpath(@__DIR__, "instances")
const OUTPUT_ROOT = joinpath(REPO_ROOT, "tmp", "positive_risk_additional_work")
const RESULT_CSV = joinpath(@__DIR__, "results_network_design_positive_additional.csv")

function local_instance_config()
    paths = sort(collect(joinpath(path, "instance.json") for path in
                        readdir(INSTANCE_ROOT; join=true)
                        if isfile(joinpath(path, "instance.json"))))
    isempty(paths) && error("No saved HNDP instances found in $INSTANCE_ROOT")
    path = joinpath(OUTPUT_ROOT, "instances.json")
    mkpath(dirname(path))
    write(path, JSON.json(Dict("saved_instances" => paths)))
    return path
end

const MODELS = [
    Dict("name"=>"ND", "model_type"=>"gbc", "partial_decomposition"=>false,
         "subproblem_method"=>"mip", "warmstart"=>false),
    Dict("name"=>"ND_WS", "model_type"=>"gbc", "partial_decomposition"=>false,
         "subproblem_method"=>"mip", "warmstart"=>true),
    # This is the bilevel-subsolver variant referred to as ND_WS_PD_BlD.
    Dict("name"=>"ND_WS_PD_BS", "model_type"=>"gbc", "partial_decomposition"=>true,
         "subproblem_method"=>"blc_jump", "warmstart"=>true),
]

function run_model!(model)
    name = String(model["name"])
    mkpath(OUTPUT_ROOT)
    model_path = joinpath(OUTPUT_ROOT, "$(lowercase(name))_models.json")
    params_path = joinpath(OUTPUT_ROOT, "$(lowercase(name))_params.json")
    write(model_path, JSON.json(Dict("models"=>[model])))
    params = Dict{String,Any}(
        "name"=>"positive_risk_additional_$(lowercase(name))_master1_workers15_1h",
        "mip_solver"=>"Gurobi", "runtime"=>3600.0, "seed"=>42,
        "pareto"=>"OPT", "warmstart"=>true, "bigMwithLC"=>false,
        "trim_coeff"=>false, "integer_obj"=>true, "enable_output_logs"=>false,
        "threads_master"=>1, "threads_sub_con"=>1,
        "parallel_subsolvers"=>15, "parallel_separation"=>true,
    )
    write(params_path, JSON.json(Dict("write_run_logs"=>false, "params"=>[params])))
    println("Running $(name): master=1, workers=15, subsolver=1")
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
    CSV.write(RESULT_CSV, reduce((a,b) -> vcat(a,b; cols=:union), frames))
end
