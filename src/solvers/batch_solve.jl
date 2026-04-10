using JSON

"""
    solve_instance!(instance_factory::Function, config_path::AbstractString, experiment_name::AbstractString)

Load an experiment configuration from `config_path`, build a fresh instance using
`instance_factory`, and solve the named experiment.
"""
function solve_instance!(
    instance_factory::Function,
    config_path::AbstractString,
    experiment_name::AbstractString,
)
    config = _read_batch_config(config_path)
    experiment = _get_named_experiment(config, experiment_name)
    return _solve_experiment(instance_factory, config, experiment)
end

"""
    solve_batch!(instance_factory::Function, config_path::AbstractString)

Load the batch configuration from `config_path` and execute all configured
experiments sequentially. A fresh instance is generated for every experiment by
calling `instance_factory()`.
"""
function solve_batch!(instance_factory::Function, config_path::AbstractString)
    config = _read_batch_config(config_path)
    results = Dict{String,RunStats}()
    aggregate_rows = Dict{String,Any}[]
    for experiment in config["experiments"]
        name = experiment["name"]
        stats, aggregate_row = _solve_experiment(instance_factory, config, experiment)
        results[name] = stats
        push!(aggregate_rows, aggregate_row)
    end
    _write_batch_summary_csv(config["output_folder_path"], aggregate_rows)
    return results
end

"""
    solve_batch!(instance_factories::AbstractVector, config_path::AbstractString)

Load the batch configuration from `config_path` and execute all configured
experiments sequentially for every passed instance factory. The full cartesian
product of instances and experiments is evaluated. Each factory is called only
for the run that is currently executed, so only one instance is alive at a time.
"""
function solve_batch!(instance_factories::AbstractVector, config_path::AbstractString)
    config = _read_batch_config(config_path)
    normalized_factories = _normalize_instance_factories(instance_factories)
    results = Dict{String,Dict{String,RunStats}}()
    aggregate_rows = Dict{String,Any}[]
    for (instance_name, instance_factory) in normalized_factories
        instance_results = Dict{String,RunStats}()
        for experiment in config["experiments"]
            if !_experiment_matches_instance(experiment, instance_name)
                continue
            end
            experiment_name = experiment["name"]
            stats, aggregate_row = _solve_experiment(
                instance_factory,
                config,
                experiment;
                instance_name=instance_name,
            )
            instance_results[experiment_name] = stats
            push!(aggregate_rows, aggregate_row)
        end
        results[instance_name] = instance_results
    end
    _write_batch_summary_csv(config["output_folder_path"], aggregate_rows)
    return results
end

function _solve_experiment(
    instance_factory::Function,
    config::Dict{String,Any},
    experiment::Dict{String,Any};
    instance_name::Union{Nothing,String}=nothing,
)
    resolved = _merge_common_and_experiment(config, experiment)
    experiment_output_path = isnothing(instance_name) ?
        joinpath(config["output_folder_path"], experiment["name"]) :
        joinpath(config["output_folder_path"], experiment["name"], instance_name)
    mkpath(experiment_output_path)
    resolved["output_folder_path"] = experiment_output_path
    if !isnothing(instance_name)
        resolved["instance_name"] = instance_name
    end

    instance, instance_metadata = _materialize_instance(instance_factory)
    solver_name = _infer_solver_name(instance)
    resolved["solver"] = solver_name
    _write_resolved_experiment_config(experiment_output_path, resolved)

    params = _build_solver_params(solver_name, resolved)
    stats = solve_instance!(instance, params)
    aggregate_row = _build_aggregate_row(stats, resolved, instance_metadata)
    return stats, aggregate_row
end

function _read_batch_config(config_path::AbstractString)
    config = JSON.parsefile(config_path)
    isa(config, Dict{String,Any}) || error("The batch configuration at $(config_path) must be a JSON object.")

    if !haskey(config, "output_folder_path")
        error("The batch configuration at $(config_path) must define 'output_folder_path'.")
    end
    if !haskey(config, "experiments")
        error("The batch configuration at $(config_path) must define 'experiments'.")
    end

    common = get(config, "common", Dict{String,Any}())
    isa(common, Dict{String,Any}) || error("The 'common' entry in $(config_path) must be a JSON object.")
    if haskey(common, "solver")
        @warn "Ignoring the JuBiC solver setting in the common batch configuration. The JuBiC solver is inferred from the generated instance."
        delete!(common, "solver")
    end

    experiments = config["experiments"]
    isa(experiments, Vector) || error("The 'experiments' entry in $(config_path) must be a JSON array.")
    length(experiments) > 0 || error("The batch configuration at $(config_path) does not contain any experiments.")

    names = Set{String}()
    normalized_experiments = Dict{String,Any}[]
    for (idx, experiment_raw) in enumerate(experiments)
        isa(experiment_raw, Dict{String,Any}) || error("Experiment $(idx) in $(config_path) must be a JSON object.")
        experiment = experiment_raw
        if !haskey(experiment, "name")
            error("Experiment $(idx) in $(config_path) is missing the required field 'name'.")
        end
        name = String(experiment["name"])
        if name in names
            error("The batch configuration at $(config_path) defines the experiment name '$(name)' multiple times.")
        end
        if haskey(experiment, "solver")
            @warn "Ignoring the JuBiC solver setting for experiment $(name). The JuBiC solver is inferred from the generated instance."
            delete!(experiment, "solver")
        end
        push!(names, name)
        push!(normalized_experiments, experiment)
    end

    return Dict(
        "output_folder_path" => String(config["output_folder_path"]),
        "common" => common,
        "experiments" => normalized_experiments,
    )
end

function _experiment_matches_instance(experiment::Dict{String,Any}, instance_name::String)
    if !haskey(experiment, "instance_name")
        return true
    end
    return String(experiment["instance_name"]) == instance_name
end

function _normalize_instance_factories(instance_factories::AbstractVector)
    normalized = Pair{String,Function}[]
    for (idx, instance_factory) in enumerate(instance_factories)
        if instance_factory isa Pair{<:AbstractString,<:Function}
            push!(normalized, String(instance_factory.first) => instance_factory.second)
        elseif instance_factory isa Function
            push!(normalized, "instance_$(idx)" => instance_factory)
        else
            error("Instance factory entry $(idx) must be either a Function or a Pair{String, Function}.")
        end
    end

    names = Set{String}()
    for (name, _) in normalized
        if name in names
            error("The list of instance factories defines the instance name '$(name)' multiple times.")
        end
        push!(names, name)
    end
    return normalized
end

function _get_named_experiment(config::Dict{String,Any}, experiment_name::AbstractString)
    for experiment in config["experiments"]
        if experiment["name"] == experiment_name
            return experiment
        end
    end
    error("The batch configuration does not define an experiment named '$(experiment_name)'.")
end

function _merge_common_and_experiment(config::Dict{String,Any}, experiment::Dict{String,Any})
    resolved = deepcopy(config["common"])
    for (key, value) in experiment
        resolved[key] = value
    end
    return resolved
end

function _write_resolved_experiment_config(output_folder_path::AbstractString, resolved::Dict{String,Any})
    config_path = joinpath(output_folder_path, "resolved_experiment.json")
    open(config_path, "w") do io
        write(io, JSON.json(resolved, 2))
    end
end

function _materialize_instance(instance_factory::Function)
    instance_data = instance_factory()
    if instance_data isa Instance
        return instance_data, Dict{String,Any}()
    elseif instance_data isa Tuple && length(instance_data) == 2
        instance = instance_data[1]
        metadata = instance_data[2]
        instance isa Instance || error("An instance factory returning a tuple must return (instance::Instance, metadata::Dict-like).")
        metadata isa AbstractDict || error("An instance factory returning a tuple must return (instance::Instance, metadata::Dict-like).")
        return instance, _string_key_dict(metadata)
    elseif instance_data isa Pair
        instance = instance_data.first
        metadata = instance_data.second
        instance isa Instance || error("An instance factory returning a Pair must return instance => metadata.")
        metadata isa AbstractDict || error("An instance factory returning a Pair must return instance => metadata.")
        return instance, _string_key_dict(metadata)
    end
    error("An instance factory must return either an Instance or (Instance, metadata::Dict-like).")
end

function _build_aggregate_row(
    stats::RunStats,
    resolved::Dict{String,Any},
    instance_metadata::AbstractDict{String,<:Any},
)
    row = Dict{String,Any}()
    for (key, value) in _filtered_config_entries(resolved)
        row[key] = _csv_safe_value(value)
    end
    for (key, value) in instance_metadata
        row["instance_" * key] = _csv_safe_value(value)
    end
    for (key, value) in stats.data
        row[key] = _csv_safe_value(value)
    end
    return row
end

function _filtered_config_entries(resolved::Dict{String,Any})
    excluded = Set([
        "output_folder_path",
    ])
    entries = Pair{String,Any}[]
    for (key, value) in resolved
        if key in excluded
            continue
        end
        push!(entries, String(key) => value)
    end
    return sort(entries; by=first)
end

function _write_batch_summary_csv(output_folder_path::AbstractString, rows::Vector{Dict{String,Any}})
    isempty(rows) && return nothing

    all_keys = Set{String}()
    for row in rows
        union!(all_keys, keys(row))
    end
    sorted_keys = sort(collect(all_keys))

    df = DataFrames.DataFrame()
    for key in sorted_keys
        df[!, key] = [get(row, key, "") for row in rows]
    end

    CSV.write(joinpath(output_folder_path, "batch_summary.csv"), df)
    return nothing
end

function _string_key_dict(dict::AbstractDict)
    return Dict(String(key) => value for (key, value) in dict)
end

function _csv_safe_value(value)
    if value isa Number || value isa AbstractString || value isa Bool || isnothing(value)
        return value
    end
    return JSON.json(value)
end

function _infer_solver_name(instance::Instance)
    if instance.master isa Master
        return "GBC"
    elseif instance.master isa BlCMaster
        return "BLC"
    elseif instance.master isa BlCLagMaster
        return "BlCLag"
    elseif instance.master isa MIPMaster
        return "MIP"
    elseif instance.master isa MibSMaster
        return "MibS"
    end
    error("Could not infer a JuBiC solver for instance master type $(typeof(instance.master)).")
end

function _build_solver_params(solver_name::AbstractString, config::Dict{String,Any})
    debbug_out = _get_bool(config, "debbug_out", false)
    output_folder_path = _required_string(config, "output_folder_path")

    if solver_name == "GBC"
        wrapper = _build_mip_solver_wrapper(config)
        file_format_output = _get_string(config, "file_format_output", "lp")
        runtime = _get_number(config, "runtime", 3600)
        threads_master = _get_int(config, "threads_master", 8)
        threads_sub_con = _get_int(config, "threads_sub_con", 8)
        pareto = _parse_pareto_cut(get(config, "pareto", "OPT"))
        warmstart = _get_bool(config, "warmstart", true)
        bigMwithLC = _get_bool(config, "bigMwithLC", false)
        trim_coeff = _get_bool(config, "trim_coeff", false)
        infinity_num = _get_number(config, "infinity_num", 1e9)
        g_round_digit = _get_int(config, "g_round_digit", 0)
        return GBCparam(
            wrapper,
            debbug_out,
            output_folder_path,
            file_format_output,
            RunStats(),
            runtime,
            threads_master,
            threads_sub_con,
            pareto,
            warmstart,
            bigMwithLC,
            trim_coeff,
            infinity_num,
            g_round_digit,
        )
    elseif solver_name == "BLC"
        wrapper = _build_mip_solver_wrapper(config)
        file_format_output = _get_string(config, "file_format_output", "lp")
        runtime = _get_number(config, "runtime", 3600)
        threads_master = _get_int(config, "threads_master", 8)
        threads_sub_con = _get_int(config, "threads_sub_con", 8)
        return BLCparam(
            wrapper,
            debbug_out,
            output_folder_path,
            file_format_output,
            RunStats(),
            runtime,
            threads_master,
            threads_sub_con,
        )
    elseif solver_name == "BlCLag"
        wrapper = _build_mip_solver_wrapper(config)
        file_format_output = _get_string(config, "file_format_output", "lp")
        runtime = _get_number(config, "runtime", 3600)
        threads_master = _get_int(config, "threads_master", 8)
        threads_sub_con = _get_int(config, "threads_sub_con", 8)
        pareto = _parse_pareto_cut(get(config, "pareto", "OPT"))
        warmstart = _get_bool(config, "warmstart", true)
        infinity_num = _get_number(config, "infinity_num", 1e9)
        return BlCLagparam(
            wrapper,
            debbug_out,
            output_folder_path,
            file_format_output,
            RunStats(),
            runtime,
            threads_master,
            threads_sub_con,
            pareto,
            warmstart,
            infinity_num,
        )
    elseif solver_name == "MIP"
        wrapper = _build_mip_solver_wrapper(config)
        file_format_output = _get_string(config, "file_format_output", "lp")
        runtime = _get_number(config, "runtime", 3600)
        threads_master = _get_int(config, "threads_master", 8)
        return MIPparam(
            wrapper,
            debbug_out,
            output_folder_path,
            file_format_output,
            RunStats(),
            runtime,
            threads_master,
        )
    elseif solver_name == "MibS"
        return MibSparam(debbug_out, output_folder_path, RunStats())
    end

    error("Unsupported inferred JuBiC solver '$(solver_name)'.")
end

function _build_mip_solver_wrapper(config::Dict{String,Any})
    mip_solver = _required_string(config, "mip_solver")
    if mip_solver == "Gurobi"
        return GurobiSolver(Gurobi.Env())
    end
    error("Unsupported MIP solver '$(mip_solver)' in experiment configuration.")
end

function _parse_pareto_cut(value)
    value_str = String(value)
    if value_str == "None" || value_str == "PARETO_NONE"
        return PARETO_NONE
    elseif value_str == "OPT" || value_str == "PARETO_OPTIMALITY_ONLY"
        return PARETO_OPTIMALITY_ONLY
    elseif value_str == "Both" || value_str == "PARETO_OPTIMALITY_AND_FEASIBILITY"
        return PARETO_OPTIMALITY_AND_FEASIBILITY
    end
    error("Unsupported pareto setting '$(value_str)' in experiment configuration.")
end

function _required_string(config::Dict{String,Any}, key::String)
    haskey(config, key) || error("The experiment configuration is missing the required field '$(key)'.")
    return String(config[key])
end

function _get_string(config::Dict{String,Any}, key::String, default::AbstractString)
    return haskey(config, key) ? String(config[key]) : String(default)
end

function _get_bool(config::Dict{String,Any}, key::String, default::Bool)
    return haskey(config, key) ? Bool(config[key]) : default
end

function _get_int(config::Dict{String,Any}, key::String, default::Integer)
    return haskey(config, key) ? Int(config[key]) : Int(default)
end

function _get_number(config::Dict{String,Any}, key::String, default::Number)
    return haskey(config, key) ? config[key] : default
end
