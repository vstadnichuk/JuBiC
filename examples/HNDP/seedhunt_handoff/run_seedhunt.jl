using JSON
using CSV
using DataFrames

include("../hndp_experiment_runner.jl")

"""
    writejson(path, obj)

Write a JSON object to `path`, creating parent folders as needed.
"""
function writejson(path, obj)
    mkpath(dirname(path))
    open(path, "w") do io
        JSON.print(io, obj, 4)
    end
end

"""
    seedhunt_configs()

Return the HNDP seed-hunt configs that were used for the original batch on the
current machine. The goal of this script is to be self-contained so it can be
copied to another PC and resumed without reconstructing old JSON files.
"""
function seedhunt_configs()
    instance_cfg = Dict(
        "parameter_seeds" => collect(1:8),
        "instances" => [
            Dict(
                "name" => "seedhunt",
                "instance_type" => "constrained_shortest_path",
                "topologies" => ["sioux_falls"],
                "nusers" => [1200, 1600, 2400],
                "length_constrained" => [true],
                "length_slack" => [0.0],
                "two_stage" => [false],
                "user_parameter_mode" => ["per_user"],
                "construction_cost" => 100,
                "max_cost" => 100,
                "max_risk" => 100,
                "max_weight" => 100,
            ),
        ],
    )

    model_cfg = Dict(
        "models" => [
            Dict(
                "name" => "gbc_mip_no_partial_ws",
                "model_type" => "gbc",
                "partial_decomposition" => false,
                "subproblem_method" => "mip",
                "big_m_mode" => "n_minus_one_most_expensive",
            ),
        ],
    )

    param_cfg = Dict(
        "write_run_logs" => false,
        "params" => [
            Dict(
                "name" => "seedhunt_params",
                "mip_solver" => "Gurobi",
                "runtime" => 60.0,
                "seed" => 42,
                "threads_master" => 2,
                "threads_sub_con" => 1,
                "file_format_output" => "lp",
                "debbug_out" => false,
                "pareto" => "OPT",
                "warmstart" => true,
            ),
        ],
    )

    return instance_cfg, model_cfg, param_cfg
end

"""
    write_seedhunt_manifests(root)

Write the generated config JSON files under `root/configs/` and return the
three file paths.
"""
function write_seedhunt_manifests(root::AbstractString)
    instance_cfg, model_cfg, param_cfg = seedhunt_configs()
    config_dir = joinpath(root, "configs")
    instance_cfg_path = joinpath(config_dir, "instances.json")
    model_cfg_path = joinpath(config_dir, "models.json")
    param_cfg_path = joinpath(config_dir, "params.json")

    writejson(instance_cfg_path, instance_cfg)
    writejson(model_cfg_path, model_cfg)
    writejson(param_cfg_path, param_cfg)
    return instance_cfg_path, model_cfg_path, param_cfg_path
end

"""
    print_seedhunt_summary(output_root)

Print a compact runtime/status summary if the batch summary CSV exists.
"""
function print_seedhunt_summary(output_root::AbstractString)
    summary_path = joinpath(output_root, "results", "batch_summary.csv")
    if !isfile(summary_path)
        println("No batch_summary.csv found yet at $(summary_path).")
        return nothing
    end

    df = CSV.read(summary_path, DataFrame)
    keep = [name for name in (
        :instance_name,
        :instance_parameter_seed,
        :instance_nusers,
        :runtime,
        :Opt_status,
        :GBCStatus,
        :gap,
    ) if name in names(df)]
    summary = select(df, keep)
    sort!(summary, [:instance_nusers, :instance_parameter_seed])
    println(summary)
    return nothing
end

"""
    run_seedhunt!(; output_root, resume=true)

Run or resume the HNDP seed-hunt batch.
"""
function run_seedhunt!(; output_root::AbstractString, resume::Bool=true, quiet::Bool=false)
    mkpath(output_root)
    instance_cfg_path, model_cfg_path, param_cfg_path = write_seedhunt_manifests(output_root)

    if !quiet
        println("Starting HNDP seed-hunt batch.")
        println("output_root = $(output_root)")
        println("resume = $(resume)")
    end

    if quiet
        log_path = joinpath(output_root, "seedhunt_console.log")
        open(log_path, "w") do io
            redirect_stdout(io) do
                redirect_stderr(io) do
                    run_hndp_experiments!(
                        instance_cfg_path,
                        model_cfg_path,
                        param_cfg_path;
                        output_root=output_root,
                        resume=resume,
                    )
                end
            end
        end
    else
        run_hndp_experiments!(
            instance_cfg_path,
            model_cfg_path,
            param_cfg_path;
            output_root=output_root,
            resume=resume,
        )
    end

    print_seedhunt_summary(output_root)
    return nothing
end

function _parse_bool_arg(value::AbstractString)
    lower = lowercase(String(value))
    if lower in ("1", "true", "yes", "y")
        return true
    elseif lower in ("0", "false", "no", "n")
        return false
    end
    throw(ArgumentError("Could not parse boolean argument '$value'. Use true/false."))
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_root = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "seedhunt_runs")
    resume = length(ARGS) >= 2 ? _parse_bool_arg(ARGS[2]) : true
    quiet = length(ARGS) >= 3 ? _parse_bool_arg(ARGS[3]) : false
    run_seedhunt!(; output_root=output_root, resume=resume, quiet=quiet)
end
