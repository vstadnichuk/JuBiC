# A wrapper for the MibS solver
using BilevelJuMP, MibS_jll
import MathOptInterface as MOI

function _set_stat!(stats::RunStats, name::String, value::Any)
    stats.data[name] = value
end

function solve_with_MibS!(inst::Instance, param::MibSparam)
    model = inst.master.model

    @debug "starting setup of the MibS solver."
    new_stat!(param.stats, "Solver", "MibSSolver")
    _set_stat!(param.stats, "MibSStatus", "Started")
    _set_stat!(param.stats, "time_limit", param.runtime)
    _set_stat!(param.stats, "MibSExecution", "JuBiCParamFileRunner")

    try
        @info model
    catch
        @debug "Could not print BilevelJuMP model for current MiBS instance to debug file."
    end

    try
        @debug "Finished model construction. Now proceeding to optimization process with MibS solver."
        solution = _solve_with_MibS_param_file_runner(model, param)

        status = solution.status
        objective = solution.objective

        if status
            _set_stat!(param.stats, "MibSStatus", "Optimal")
            new_stat!(param.stats, "Opt", objective)
        else
            _set_stat!(param.stats, "MibSStatus", "NoOptimalSolution")
            @warn "The MibS solver did not find an optimal solution."
        end

        process_mibs_log_text!(param.stats, solution.log_output)
    catch e
        error_message = sprint(showerror, e)
        _set_stat!(param.stats, "MibSStatus", "Error")
        _set_stat!(param.stats, "Error", error_message)
        @error "An error occurred while solving with MibS: $error_message"
    end
end

function _solve_with_MibS_param_file_runner(model::BilevelJuMP.BilevelModel, param::MibSparam)
    orig_path = pwd()
    mktempdir() do path
        mps_filename = joinpath(path, "model.mps")
        aux_filename = joinpath(path, "model.aux")
        par_filename = joinpath(path, "mibs.par")
        output_filename = joinpath(path, "mibs_output.txt")
        error_filename = joinpath(path, "mibs_errors.txt")

        _validate_mibs_model_is_mip_mip(model)
        new_model, variables, objective, constraints, sense =
            BilevelJuMP._build_single_model(model, true)

        MOI.write_to_file(new_model, mps_filename)
        BilevelJuMP._write_auxillary_file(
            new_model,
            variables,
            objective,
            constraints,
            sense,
            aux_filename,
        )
        _write_mibs_parameter_file!(par_filename, mps_filename, aux_filename, param)

        output, err = _call_mibs_with_param_file(
            par_filename,
            output_filename,
            error_filename,
        )

        if should_debbug_print(param) || length(err) > 0
            _copy_mibs_debug_files!(
                path,
                output_file_path(param);
                prefix = "custom_",
            )
        end

        if length(err) > 0
            error(
                "MibS returned:\n\n" *
                "$err\n\n" *
                "JuBiC kept the debug files in $(output_file_path(param)) when debug mode is enabled or when errors occur.",
            )
        end
        if length(output) == 0
            error("MibS failed to return")
        end

        parsed = BilevelJuMP._parse_output(output, new_model, variables)
        return (
            status = parsed.status,
            objective = parsed.objective,
            nonzero_upper = parsed.nonzero_upper,
            nonzero_lower = parsed.nonzero_lower,
            all_upper = parsed.all_upper,
            all_lower = parsed.all_lower,
            all_var = parsed.all_var,
            log_output = output,
            working_directory = orig_path,
        )
    end
end

function _validate_mibs_model_is_mip_mip(model::BilevelJuMP.BilevelModel)
    non_mip_variables = String[]
    for var in all_variables(model)
        if !(JuMP.is_binary(var) || JuMP.is_integer(var))
            push!(non_mip_variables, JuMP.name(var))
        end
    end

    if !isempty(non_mip_variables)
        shown_vars = join(first(sort(non_mip_variables), min(8, length(non_mip_variables))), ", ")
        suffix = length(non_mip_variables) > 8 ? ", ..." : ""
        error(
            "MiBS currently supports only MIP-MIP models. JuBiC detected $(length(non_mip_variables)) non-integer variable(s) before exporting the model. " *
            "Examples: $(shown_vars)$(suffix). " *
            "Please ensure that all upper- and lower-level variables are binary or integer before calling MiBS.",
        )
    end

    return nothing
end

function _write_mibs_parameter_file!(
    filename::String,
    mps_filename::String,
    aux_filename::String,
    param::MibSparam,
)
    open(filename, "w") do io
        println(io, "Alps_instance ", mps_filename)
        println(io, "MibS_auxiliaryInfoFile ", aux_filename)
        println(io, "Alps_timeLimit ", param.runtime)
    end
    return nothing
end

function _call_mibs_with_param_file(
    par_filename::String,
    output_filename::String,
    error_filename::String,
)
    run_error = ""
    try
        MibS_jll.mibs() do exe
            run(
                pipeline(
                    `$(exe) -param $(par_filename)`;
                    stdout = output_filename,
                    stderr = error_filename,
                ),
            )
        end
    catch e
        run_error = sprint(showerror, e)
    end

    output = isfile(output_filename) ? read(output_filename, String) : ""
    err = isfile(error_filename) ? read(error_filename, String) : ""
    if length(err) == 0 && length(run_error) > 0
        err = run_error
    end
    return output, err
end

function _copy_mibs_debug_files!(source_dir::String, target_dir; prefix::String = "")
    mkpath(target_dir)
    for file_name in ("model.mps", "model.aux", "mibs.par", "mibs_output.txt", "mibs_errors.txt")
        source = joinpath(source_dir, file_name)
        if isfile(source)
            cp(source, joinpath(target_dir, prefix * file_name); force = true)
        end
    end
    return nothing
end

function process_mibs_log_text!(stats, logtext::String)
    match_nodes = match(r"Number of nodes processed:\s+(\d+)", logtext)
    match_runtime = match(r"Search wall-clock time:\s+([0-9.]+)", logtext)
    match_gap = match(r"Relative optimality gap is\s+([0-9.]+)%", logtext)

    nodes = isnothing(match_nodes) ? missing : parse(Int, match_nodes.captures[1])
    runtime = isnothing(match_runtime) ? missing : parse(Float64, match_runtime.captures[1])
    gap = isnothing(match_gap) ? missing : parse(Float64, match_gap.captures[1])

    _set_stat!(stats, "BNodes", nodes)
    _set_stat!(stats, "runtime", runtime)
    _set_stat!(stats, "gap", gap)
end
