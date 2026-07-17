using JuBiC
using JSON
using Dates
using SHA
using CSV
using DataFrames

include("hndp_model_generation.jl")

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

    execution_mode = String(get(param_cfg, "execution_mode", "subprocess_per_experiment"))
    if execution_mode == "subprocess_per_experiment"
        return _run_hndp_experiments_subprocess!(
            instance_cfg,
            model_cfg,
            param_cfg,
            String(instance_config_path),
            String(model_config_path),
            String(param_config_path);
            output_root=String(output_root),
            resume=resume,
        )
    elseif execution_mode != "in_process"
        throw(ArgumentError("Unsupported HNDP execution_mode '$(execution_mode)'. Supported values are 'subprocess_per_experiment' and 'in_process'."))
    end

    return _run_hndp_experiments_in_process!(
        instance_cfg,
        model_cfg,
        param_cfg,
        String(instance_config_path),
        String(model_config_path),
        String(param_config_path);
        output_root=String(output_root),
        resume=resume,
    )
end

function _run_hndp_experiments_in_process!(
    instance_cfg::Dict{String,Any},
    model_cfg::Dict{String,Any},
    param_cfg::Dict{String,Any},
    instance_config_path::String,
    model_config_path::String,
    param_config_path::String;
    output_root::String,
    resume::Bool=true,
)
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

                run_output_path = _allocate_hndp_run_output_path(
                    output_root,
                    experiment_id;
                    write_run_logs=write_run_logs,
                )

                aggregate_row = Dict{String,Any}()
                try
                    if Bool(get(model_spec, "write_solver_instance_files", false))
                        _maybe_export_hndp_solver_instance!(
                            generated_network,
                            model_spec,
                            param_spec,
                            output_root,
                        )
                    end
                    stats, aggregate_row = _run_hndp_experiment(
                        generated_network,
                        model_spec,
                        param_spec,
                        run_output_path,
                        write_run_logs;
                        force_post_cleanup=true,
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
                        try
                            rm(run_output_path; recursive=true, force=true, allow_delayed_delete=true)
                        catch err
                            @warn "Could not immediately remove transient HNDP run folder $(run_output_path). Continuing. Error: $(sprint(showerror, err))"
                        end
                    end
                end
            end
        end
    end)

    return results
end

function _run_hndp_experiments_subprocess!(
    instance_cfg::Dict{String,Any},
    model_cfg::Dict{String,Any},
    param_cfg::Dict{String,Any},
    instance_config_path::String,
    model_config_path::String,
    param_config_path::String;
    output_root::String,
    resume::Bool=true,
)
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
    subprocess_log_root = joinpath(output_root, "results", "subprocess_logs")
    mkpath(subprocess_log_root)

    visit_hndp_networks(instance_cfg, generated_network -> begin
        instance_metadata = deepcopy(generated_network.metadata)
        instance_name = String(instance_metadata["name"])

        for (model_idx, model_spec) in enumerate(models)
            for (param_idx, param_spec) in enumerate(params)
                experiment_id = _hndp_experiment_id(instance_name, model_spec, param_spec)
                if experiment_id in completed_ids
                    continue
                end

                if Bool(get(model_spec, "write_solver_instance_files", false))
                    _maybe_export_hndp_solver_instance!(
                        generated_network,
                        model_spec,
                        param_spec,
                        output_root,
                    )
                end

                ok = _run_hndp_experiment_subprocess!(
                    instance_name,
                    model_idx,
                    param_idx,
                    instance_config_path,
                    model_config_path,
                    param_config_path,
                    output_root,
                    experiment_id,
                    subprocess_log_root,
                )

                if !ok
                    aggregate_row = _build_hndp_error_row(
                        generated_network,
                        model_spec,
                        param_spec,
                        experiment_id,
                        ErrorException("HNDP child process terminated before producing a result row. See subprocess logs."),
                    )
                    aggregate_row["experiment_id"] = experiment_id
                    aggregate_row["instance_name"] = instance_name
                    aggregate_row["RunStatus"] = "SubprocessError"
                    JuBiC._append_batch_summary_csv!(summary_csv_path, aggregate_row)
                end

                push!(completed_ids, experiment_id)
            end
        end
    end)

    return Dict{String,JuBiC.RunStats}()
end

function _maybe_export_hndp_solver_instance!(
    generated_network::HNDPGeneratedNetwork,
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
    output_root::AbstractString,
)
    solver_wrapper = _build_hndp_mip_wrapper(param_spec)
    instance_name = String(generated_network.metadata["name"])
    export_root = joinpath(output_root, "mibs_instances")
    export_id = bytes2hex(sha1(instance_name))[1:16]
    export_dir = joinpath(export_root, export_id)
    export_basename = "instance"
    mps_path = joinpath(export_dir, export_basename * ".mps")
    aux_path = joinpath(export_dir, export_basename * ".aux")
    meta_path = joinpath(export_dir, "metadata.json")

    if isfile(mps_path) && isfile(aux_path) && isfile(meta_path)
        return nothing
    end

    instance = build_hndp_mibs_instance(
        generated_network.instance,
        solver_wrapper;
        partial_decomposition=false,
    )

    mkpath(export_dir)
    try
        JuBiC.output_MibS_instance(instance, export_basename, export_dir)
        open(meta_path, "w") do io
            write(
                io,
                JSON.json(
                    Dict(
                        "instance_name" => instance_name,
                        "export_id" => export_id,
                        "files" => Dict(
                            "mps" => export_basename * ".mps",
                            "aux" => export_basename * ".aux",
                        ),
                    ),
                    2,
                ),
            )
        end
    catch err
        @warn "Could not export the MiBS-format HNDP instance for $(instance_name) into $(export_dir): $(sprint(showerror, err))"
    end
    return nothing
end

function _legacy_maybe_export_hndp_solver_instance!(
    generated_network::HNDPGeneratedNetwork,
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
    run_output_path::AbstractString,
)
    solver_wrapper = _build_hndp_mip_wrapper(param_spec)
    model_instance, _ = _build_hndp_model_from_spec(
        generated_network.instance,
        generated_network.metadata,
        solver_wrapper,
        model_spec,
        param_spec,
    )

    solver_name = JuBiC._infer_solver_name(model_instance)
    if solver_name != "GBC"
        @warn "Skipping solver-instance export for HNDP experiment $(get(model_spec, "name", get(model_spec, "model_type", "model"))). The current export hook only supports GBC instances."
        return nothing
    end

    try
        output_GBC_solver_instance(model_instance, run_output_path)
    catch err
        @warn "Could not export the GBC solver instance for HNDP run folder $(run_output_path): $(sprint(showerror, err))"
    end
    return nothing
end

function _run_hndp_experiment(
    generated_network::HNDPGeneratedNetwork,
    model_spec::Dict{String,Any},
    param_spec::Dict{String,Any},
    run_output_path::AbstractString,
    write_run_logs::Bool,
    ;
    force_post_cleanup::Bool=true,
)
    hndp = generated_network.instance
    solver_wrapper = _build_hndp_mip_wrapper(param_spec)
    model_instance = nothing
    model_metadata = Dict{String,Any}()
    params = nothing
    stats = nothing
    aggregate_row = Dict{String,Any}()

    try
        model_instance, model_metadata = _build_hndp_model_from_spec(
            hndp,
            generated_network.metadata,
            solver_wrapper,
            model_spec,
            param_spec,
        )

        solver_name = JuBiC._infer_solver_name(model_instance)
        resolved_param_config = _resolve_hndp_param_config(
            param_spec,
            solver_name,
            run_output_path,
            model_metadata,
            write_run_logs,
        )
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
    finally
        model_instance = nothing
        params = nothing
        solver_wrapper = nothing
        force_post_cleanup && JuBiC._run_post_gurobi_cleanup!()
    end
end

function _run_hndp_experiment_subprocess!(
    instance_name::String,
    model_idx::Int,
    param_idx::Int,
    instance_config_path::String,
    model_config_path::String,
    param_config_path::String,
    output_root::String,
    experiment_id::String,
    subprocess_log_root::String,
)
    stdout_path = joinpath(subprocess_log_root, experiment_id * ".stdout.log")
    stderr_path = joinpath(subprocess_log_root, experiment_id * ".stderr.log")
    cmd = `$(Base.julia_cmd()) --project=$(Base.active_project()) $(abspath(@__FILE__))`
    env = Dict(
        "JUBIC_HNDP_CHILD_MODE" => "1",
        "JUBIC_HNDP_INSTANCE_CONFIG" => instance_config_path,
        "JUBIC_HNDP_MODEL_CONFIG" => model_config_path,
        "JUBIC_HNDP_PARAM_CONFIG" => param_config_path,
        "JUBIC_HNDP_OUTPUT_ROOT" => output_root,
        "JUBIC_HNDP_INSTANCE_NAME" => instance_name,
        "JUBIC_HNDP_MODEL_IDX" => string(model_idx),
        "JUBIC_HNDP_PARAM_IDX" => string(param_idx),
        "JUBIC_HNDP_EXPERIMENT_ID" => experiment_id,
    )

    proc = run(
        pipeline(
            ignorestatus(setenv(cmd, env));
            stdout=stdout_path,
            stderr=stderr_path,
        ),
    )
    summary_csv_path = joinpath(output_root, "results", "batch_summary.csv")
    row_written = _summary_contains_experiment_id(summary_csv_path, experiment_id)

    if success(proc)
        return row_written
    end

    if row_written
        @warn "HNDP child process for $(experiment_id) exited with a nonzero status after already writing its CSV row. Treating the experiment as completed and skipping the parent-side SubprocessError duplicate row."
        return true
    end
    return false
end

function _summary_contains_experiment_id(summary_csv_path::String, experiment_id::String)
    isfile(summary_csv_path) || return false
    df = CSV.read(summary_csv_path, DataFrame)
    if !("experiment_id" in names(df))
        return false
    end
    return any(isequal(experiment_id), df[!, "experiment_id"])
end

function _run_hndp_single_child_job_from_env!()
    instance_config_path = ENV["JUBIC_HNDP_INSTANCE_CONFIG"]
    model_config_path = ENV["JUBIC_HNDP_MODEL_CONFIG"]
    param_config_path = ENV["JUBIC_HNDP_PARAM_CONFIG"]
    output_root = ENV["JUBIC_HNDP_OUTPUT_ROOT"]
    instance_name = ENV["JUBIC_HNDP_INSTANCE_NAME"]
    model_idx = parse(Int, ENV["JUBIC_HNDP_MODEL_IDX"])
    param_idx = parse(Int, ENV["JUBIC_HNDP_PARAM_IDX"])
    experiment_id = ENV["JUBIC_HNDP_EXPERIMENT_ID"]

    instance_cfg = load_hndp_network_generation_config(instance_config_path)
    model_cfg = JSON.parsefile(model_config_path)
    param_cfg = JSON.parsefile(param_config_path)
    models = _load_hndp_model_specs(model_cfg)
    params = _load_hndp_param_specs(param_cfg)
    model_spec = models[model_idx]
    param_spec = params[param_idx]
    write_run_logs = get(param_cfg, "write_run_logs", true)
    summary_csv_path = joinpath(output_root, "results", "batch_summary.csv")

    found = Ref(false)
    visit_hndp_networks(instance_cfg, generated_network -> begin
        found[] && return nothing
        current_name = String(generated_network.metadata["name"])
        current_name == instance_name || return nothing
        found[] = true

        aggregate_row = Dict{String,Any}()
        run_output_path = _allocate_hndp_run_output_path(
            output_root,
            experiment_id;
            write_run_logs=write_run_logs,
        )

        try
            _, aggregate_row = _run_hndp_experiment(
                generated_network,
                model_spec,
                param_spec,
                run_output_path,
                write_run_logs;
                force_post_cleanup=false,
            )
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

            if !write_run_logs
                try
                    rm(run_output_path; recursive=true, force=true, allow_delayed_delete=true)
                catch err
                    @warn "Could not immediately remove transient HNDP run folder $(run_output_path). Continuing. Error: $(sprint(showerror, err))"
                end
            end
        end
        return nothing
    end)

    found[] || error("HNDP child job could not find generated instance named '$(instance_name)'.")
    return nothing
end

function _build_hndp_model_from_spec(
    hndp::HNDPwC,
    instance_metadata::Dict{String,Any},
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
        return instance, hndp_model_size_metadata(instance)
    elseif model_type == "blclag"
        instance = build_hndp_blclag_instance(
            hndp,
            solver_wrapper;
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            subproblem_method=_parse_hndp_subproblem_method(get(model_spec, "subproblem_method", "blc_jump")),
        )
        return instance, hndp_model_size_metadata(instance)
    elseif model_type == "gbc"
        instance = build_hndp_gbc_instance(
            hndp,
            solver_wrapper;
            partial_decomposition=Bool(get(model_spec, "partial_decomposition", true)),
            include_objL2=Bool(get(model_spec, "include_objL2", false)),
            subproblem_method=_parse_hndp_subproblem_method(get(model_spec, "subproblem_method", "mip")),
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
        )
        return instance, hndp_model_size_metadata(instance)
    elseif model_type == "mibs"
        instance = build_hndp_mibs_instance(
            hndp,
            solver_wrapper;
            partial_decomposition=Bool(get(model_spec, "partial_decomposition", false)),
        )
        return instance, hndp_model_size_metadata(instance)
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
            use_decision_arc_dominance=Bool(get(model_spec, "use_decision_arc_dominance", true)),
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
        instance, enum_runtime, path_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_instance(
            hndp,
            solver_wrapper;
            enumeration_time_limit=enum_limit,
            fallback_mode=HNDP_HYBRID_FALLBACK_SD,
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            indicator_constraints=Bool(get(model_spec, "indicator_constraints", false)),
            bound_duals=Bool(get(model_spec, "bound_duals", true)),
            use_decision_arc_dominance=Bool(get(model_spec, "use_decision_arc_dominance", true)),
            fallback_if_paths_exceed_flow_vars=Bool(get(model_spec, "fallback_if_paths_exceed_flow_vars", false)),
        )
        remaining_runtime = max(total_runtime - enum_runtime, 0.0)
        metadata = Dict{String,Any}(
            "enum_runtime" => enum_runtime,
            "path_counts" => path_counts,
            "fallback_users" => fallback_users,
            "path_count_fallback_users" => path_count_fallback_users,
            "path_count_fallback_user_count" => length(path_count_fallback_users),
            "runtime_override" => remaining_runtime,
        )
        return instance, metadata
    elseif model_type == "hybrid_blc"
        enum_limit = _effective_hndp_enumeration_limit(model_spec, total_runtime)
        instance, enum_runtime, path_counts, fallback_users, path_count_fallback_users = build_hndp_hybrid_blc_instance(
            hndp,
            solver_wrapper;
            enumeration_time_limit=enum_limit,
            big_m_mode=_parse_hndp_big_m_mode(get(model_spec, "big_m_mode", "fixed_network_path")),
            subproblem_method=_parse_hndp_subproblem_method(get(model_spec, "subproblem_method", "mip")),
            use_decision_arc_dominance=Bool(get(model_spec, "use_decision_arc_dominance", true)),
            fallback_if_paths_exceed_flow_vars=Bool(get(model_spec, "fallback_if_paths_exceed_flow_vars", false)),
            availability_budget_fraction=nothing,
            availability_budget_count=get(instance_metadata, "availability_budget_count", nothing),
        )
        remaining_runtime = max(total_runtime - enum_runtime, 0.0)
        metadata = Dict{String,Any}(
            "enum_runtime" => enum_runtime,
            "path_counts" => path_counts,
            "fallback_users" => fallback_users,
            "path_count_fallback_users" => path_count_fallback_users,
            "path_count_fallback_user_count" => length(path_count_fallback_users),
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
    write_run_logs::Bool,
)
    resolved = deepcopy(param_spec)
    resolved["solver"] = solver_name
    resolved["output_folder_path"] = String(run_output_path)
    resolved["enable_output_logs"] = write_run_logs
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
    instance_slug = _short_slug(instance_name, 32)
    model_slug = _short_slug(model_name, 20)
    param_slug = _short_slug(param_name, 20)
    digest = bytes2hex(sha1(string(instance_name, "::", model_name, "::", param_name)))[1:12]
    return string(instance_slug, "__", model_slug, "__", param_slug, "__", digest)
end

function _allocate_hndp_run_output_path(
    output_root::AbstractString,
    experiment_id::AbstractString;
    write_run_logs::Bool,
)
    runs_root = joinpath(String(output_root), "runs")
    target_path = write_run_logs ?
        joinpath(runs_root, String(experiment_id)) :
        joinpath(runs_root, "_transient", String(experiment_id))

    try
        mkpath(target_path)
        return target_path
    catch err
        @warn "Could not create HNDP run folder $(target_path). Error: $(sprint(showerror, err)). Falling back to $(runs_root)."
    end

    try
        mkpath(runs_root)
        return runs_root
    catch err
        @warn "Could not create fallback HNDP run folder $(runs_root). Error: $(sprint(showerror, err)). Falling back to $(output_root)."
    end

    if isdir(String(output_root))
        @warn "Could not create a dedicated HNDP run folder under $(runs_root). Logging will use the explicit output root $(output_root)."
        return String(output_root)
    end

    error("Could not allocate an HNDP run output folder under $(output_root).")
end

function _slugify(text::AbstractString)
    return replace(lowercase(String(text)), r"[^a-z0-9]+" => "_")
end

function _short_slug(text::AbstractString, max_len::Int)
    slug = strip(_slugify(text), '_')
    isempty(slug) && return "x"
    return first(slug, min(max_len, ncodeunits(slug)))
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
    elseif value_str == "mibs"
        return HNDP_SUBPROBLEM_MIBS
    end
    throw(ArgumentError("Unsupported HNDP subproblem_method '$(value_str)'."))
end

function _build_hndp_mip_wrapper(param_spec::Dict{String,Any})
    mip_solver = String(get(param_spec, "mip_solver", "Gurobi"))
    mip_solver == "Gurobi" || throw(ArgumentError("Unsupported HNDP mip_solver '$(mip_solver)'."))    
    return GurobiSolver()
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__) && get(ENV, "JUBIC_HNDP_CHILD_MODE", "0") == "1"
    _run_hndp_single_child_job_from_env!()
end
