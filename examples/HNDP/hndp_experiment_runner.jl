using JuBiC
using JSON
using Dates

include("hndp_model_generation_v2.jl")

"""
    run_hndp_experiments!(
        instance_config_path,
        model_config_path,
        param_config_path;
        output_root,
        resume=true,
    )

Run a streamed HNDP experiment pipeline based on three JSON files:
- one for network instance generation
- one for HNDP model generation choices
- one for JuBiC solve parameters

The runner generates exactly one HNDP instance at a time, builds one JuBiC
model, solves it, appends a single CSV row to disk, and then continues with the
next combination. This keeps memory usage bounded and preserves completed runs
even if the machine crashes in the middle of the batch.
"""
function run_hndp_experiments!(
    instance_config_path::AbstractString,
    model_config_path::AbstractString,
    param_config_path::AbstractString;
    output_root::AbstractString,
    resume::Bool=true,
)
    instance_cfg = load_hndp_network_generation_config(String(instance_config_path))
    model_cfg = JSON.parsefile(String(model_config_path))
    param_cfg = JSON.parsefile(String(param_config_path))

    models = _load_hndp_model_specs(model_cfg)
    params = _load_hndp_param_specs(param_cfg)
    write_run_logs = get(param_cfg, "write_run_logs", true)

    _prepare_hndp_experiment_output!(
        output_root,
        instance_config_path,
        model_config_path,
        param_config_path,
    )

    summary_csv_path = joinpath(output_root, "results", "batch_summary.csv")
    completed_ids = resume ? JuBiC._read_completed_batch_ids(summary_csv_path) : Set{String}()
    results = Dict{String,JuBiC.RunStats}()

    visit_hndp_networks(instance_cfg, generated_network -> begin
        instance_metadata = deepcopy(generated_network.metadata)
        instance_name = String(instance_metadata["name"])

        for model_spec in models
            for param_spec in params
                experiment_id = _hndp_experiment_id(instance_name, model_spec, param_spec)
                if experiment_id in completed_ids
                    continue
                end

                run_output_path = write_run_logs ?
                    joinpath(output_root, "runs", experiment_id) :
                    mktempdir()
                if write_run_logs
                    mkpath(run_output_path)
                end

                aggregate_row = Dict{String,Any}()
                try
                    stats, aggregate_row = _run_hndp_experiment(
                        generated_network,
                        model_spec,
                        param_spec,
                        run_output_path,
                    )
                    results[experiment_id] = stats
                catch err
                    aggregate_row = _build_hndp_error_row(
                        generated_network,
                        model_spec,
                        param_spec,
                        experiment_id,
                        err,
                    )
                finally
                    aggregate_row["experiment_id"] = experiment_id
                    aggregate_row["instance_name"] = instance_name
                    JuBiC._append_batch_summary_csv!(summary_csv_path, aggregate_row)
                    push!(completed_ids, experiment_id)

                    if !write_run_logs
                        rm(run_output_path; recursive=true, force=true)
                    end
                end
            end
        end
    end)

    return results
end

function _run_hndp_experiment(
    generated_network::HNDPGeneratedNetwork,
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
    run_output_path::AbstractString,
)
    hndp = generated_network.instance
    solver_wrapper = _build_hndp_mip_wrapper(param_spec)

    model_instance, model_metadata = _build_hndp_model_from_spec(
        hndp,
        solver_wrapper,
        model_spec,
        param_spec,
    )

    solver_name = JuBiC._infer_solver_name(model_instance)
    resolved_param_config = _resolve_hndp_param_config(param_spec, solver_name, run_output_path, model_metadata)
    params = JuBiC._build_solver_params(solver_name, resolved_param_config)
    stats = solve_instance!(model_instance, params)

    aggregate_row = _build_hndp_aggregate_row(
        stats,
        generated_network.metadata,
        model_spec,
        param_spec,
        model_metadata,
    )
    return stats, aggregate_row
end

function _build_hndp_model_from_spec(
    hndp::HNDPwC,
    solver_wrapper::SolverWrapper,
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
)
    model_type = String(model_spec["model_type"])
    total_runtime = Float64(get(param_spec, "runtime", 3600))

    if model_type == "blc"
        instance = build_hndp_blc_instance(
            hndp,
            solver_wrapper;
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            subproblem_method=_parse_hndp_subproblem_method(get(model_spec, "subproblem_method", "mip")),
        )
        return instance, Dict{String,Any}()
    elseif model_type == "gbc"
        instance = build_hndp_gbc_instance(
            hndp,
            solver_wrapper;
            partial_decomposition=Bool(get(model_spec, "partial_decomposition", true)),
            include_objL2=Bool(get(model_spec, "include_objL2", false)),
            subproblem_method=_parse_hndp_subproblem_method(get(model_spec, "subproblem_method", "mip")),
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
        )
        return instance, Dict{String,Any}()
    elseif model_type == "mibs"
        instance = build_hndp_mibs_instance(
            hndp,
            solver_wrapper;
            partial_decomposition=Bool(get(model_spec, "partial_decomposition", false)),
        )
        return instance, Dict{String,Any}()
    elseif model_type == "sd"
        instance = build_hndp_sd_instance(
            hndp,
            solver_wrapper;
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            indicator_constraints=Bool(get(model_spec, "indicator_constraints", false)),
            bound_duals=Bool(get(model_spec, "bound_duals", true)),
        )
        return instance, Dict{String,Any}()
    elseif model_type == "path"
        enum_limit = _effective_hndp_enumeration_limit(model_spec, total_runtime)
        instance, enum_runtime, path_counts = build_hndp_path_instance(
            hndp,
            solver_wrapper;
            enumeration_time_limit=enum_limit,
            parallelize=Bool(get(model_spec, "parallelize", false)),
        )
        remaining_runtime = max(total_runtime - enum_runtime, 0.0)
        metadata = Dict{String,Any}(
            "enum_runtime" => enum_runtime,
            "path_counts" => path_counts,
            "runtime_override" => remaining_runtime,
        )
        return instance, metadata
    elseif model_type == "hybrid_sd"
        enum_limit = _effective_hndp_enumeration_limit(model_spec, total_runtime)
        instance, enum_runtime, path_counts, fallback_users = build_hndp_hybrid_instance(
            hndp,
            solver_wrapper;
            enumeration_time_limit=enum_limit,
            fallback_mode=HNDP_HYBRID_FALLBACK_SD,
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            indicator_constraints=Bool(get(model_spec, "indicator_constraints", false)),
            bound_duals=Bool(get(model_spec, "bound_duals", true)),
        )
        remaining_runtime = max(total_runtime - enum_runtime, 0.0)
        metadata = Dict{String,Any}(
            "enum_runtime" => enum_runtime,
            "path_counts" => path_counts,
            "fallback_users" => fallback_users,
            "runtime_override" => remaining_runtime,
        )
        return instance, metadata
    elseif model_type == "hybrid_blc"
        enum_limit = _effective_hndp_enumeration_limit(model_spec, total_runtime)
        instance, enum_runtime, path_counts, fallback_users = build_hndp_hybrid_blc_instance(
            hndp,
            solver_wrapper;
            enumeration_time_limit=enum_limit,
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            subproblem_method=_parse_hndp_subproblem_method(get(model_spec, "subproblem_method", "mip")),
        )
        remaining_runtime = max(total_runtime - enum_runtime, 0.0)
        metadata = Dict{String,Any}(
            "enum_runtime" => enum_runtime,
            "path_counts" => path_counts,
            "fallback_users" => fallback_users,
            "runtime_override" => remaining_runtime,
        )
        return instance, metadata
    end

    throw(ArgumentError("Unsupported HNDP model_type '$(model_type)'."))
end

function _effective_hndp_enumeration_limit(model_spec::Dict{String,Any}, total_runtime::Float64)
    requested = Float64(get(model_spec, "enumeration_time_limit", total_runtime))
    if requested > total_runtime
        @warn "The requested HNDP enumeration_time_limit $(requested) exceeds the total runtime $(total_runtime). Clamping it to the total runtime."
        return total_runtime
    end
    return requested
end

function _resolve_hndp_param_config(
    param_spec::Dict{String,Any},
    solver_name::AbstractString,
    run_output_path::AbstractString,
    model_metadata::Dict{String,Any},
)
    resolved = deepcopy(param_spec)
    resolved["solver"] = solver_name
    resolved["output_folder_path"] = String(run_output_path)
    if haskey(model_metadata, "runtime_override")
        resolved["runtime"] = model_metadata["runtime_override"]
    end
    return resolved
end

function _build_hndp_aggregate_row(
    stats::JuBiC.RunStats,
    instance_metadata::Dict{String,Any},
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
    model_metadata::Dict{String,Any},
)
    row = Dict{String,Any}()
    for (key, value) in instance_metadata
        row["instance_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    for (key, value) in model_spec
        row["model_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    for (key, value) in param_spec
        if key == "output_folder_path" || key == "write_run_logs"
            continue
        end
        row["param_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    for (key, value) in model_metadata
        row["model_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    for (key, value) in stats.data
        row[String(key)] = JuBiC._csv_safe_value(value)
    end
    return row
end

function _build_hndp_error_row(
    generated_network::HNDPGeneratedNetwork,
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
    experiment_id::String,
    err,
)
    row = Dict{String,Any}(
        "experiment_id" => experiment_id,
        "Error" => sprint(showerror, err),
        "RunStatus" => "Error",
    )
    for (key, value) in generated_network.metadata
        row["instance_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    for (key, value) in model_spec
        row["model_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    for (key, value) in param_spec
        if key == "output_folder_path" || key == "write_run_logs"
            continue
        end
        row["param_" * String(key)] = JuBiC._csv_safe_value(value)
    end
    return row
end

function _prepare_hndp_experiment_output!(
    output_root::AbstractString,
    instance_config_path::AbstractString,
    model_config_path::AbstractString,
    param_config_path::AbstractString,
)
    mkpath(joinpath(output_root, "results"))
    mkpath(joinpath(output_root, "manifests"))
    mkpath(joinpath(output_root, "runs"))

    cp(String(instance_config_path), joinpath(output_root, "manifests", "instances.json"); force=true)
    cp(String(model_config_path), joinpath(output_root, "manifests", "models.json"); force=true)
    cp(String(param_config_path), joinpath(output_root, "manifests", "params.json"); force=true)

    manifest = Dict(
        "created_at_utc" => string(Dates.now(Dates.UTC)),
        "instance_config" => basename(String(instance_config_path)),
        "model_config" => basename(String(model_config_path)),
        "param_config" => basename(String(param_config_path)),
    )
    open(joinpath(output_root, "manifests", "run_manifest.json"), "w") do io
        write(io, JSON.json(manifest, 2))
    end
    return nothing
end

function _load_hndp_model_specs(model_cfg::Dict{String,Any})
    specs = get(model_cfg, "models", nothing)
    specs === nothing && throw(ArgumentError("The HNDP model configuration requires a 'models' field."))
    return [Dict{String,Any}(spec) for spec in specs]
end

function _load_hndp_param_specs(param_cfg::Dict{String,Any})
    specs = get(param_cfg, "params", nothing)
    specs === nothing && throw(ArgumentError("The HNDP parameter configuration requires a 'params' field."))
    return [Dict{String,Any}(spec) for spec in specs]
end

function _hndp_experiment_id(instance_name::String, model_spec::Dict{String,Any}, param_spec::Dict{String,Any})
    model_name = String(get(model_spec, "name", get(model_spec, "model_type", "model")))
    param_name = String(get(param_spec, "name", "params"))
    return string(_slugify(instance_name), "__", _slugify(model_name), "__", _slugify(param_name))
end

function _slugify(text::AbstractString)
    return replace(lowercase(String(text)), r"[^a-z0-9]+" => "_")
end

function _parse_hndp_big_m_mode(value)
    value_str = String(value)
    if value_str == "fixed_network_path"
        return HNDP_BIGM_FIXED_NETWORK_PATH
    elseif value_str == "n_minus_one_most_expensive"
        return HNDP_BIGM_N_MINUS_ONE
    end
    throw(ArgumentError("Unsupported HNDP big_m_mode '$(value_str)'."))
end

function _parse_hndp_subproblem_method(value)
    value_str = String(value)
    if value_str == "mip"
        return HNDP_SUBPROBLEM_MIP
    elseif value_str == "astar"
        return HNDP_SUBPROBLEM_ASTAR
    elseif value_str == "blc_jump"
        return HNDP_SUBPROBLEM_BLC_JUMP
    end
    throw(ArgumentError("Unsupported HNDP subproblem_method '$(value_str)'."))
end

function _build_hndp_mip_wrapper(param_spec::Dict{String,Any})
    mip_solver = String(get(param_spec, "mip_solver", "Gurobi"))
    mip_solver == "Gurobi" || throw(ArgumentError("Unsupported HNDP mip_solver '$(mip_solver)'."))
    return GurobiSolver()
end
