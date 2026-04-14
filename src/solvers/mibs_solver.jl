# A wrapper for the MibS solver
using BilevelJuMP, MibS_jll

function _set_stat!(stats::RunStats, name::String, value::Any)
    stats.data[name] = value
end

function solve_with_MibS!(inst::Instance, param::MibSparam)
    model = inst.master.model

    # register solver
    @debug "starting setup of the MibS solver."
    new_stat!(param.stats, "Solver", "MibSSolver")
    _set_stat!(param.stats, "MibSStatus", "Started")

    try
        @info model
    catch
        @debug "Could not print BilevelJuMP model for current MiBS instance to debug file."
    end

    try
        @debug "Finished model construction. Now proceeding to optimization process with MibS solver."
        solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs)
        status = solution.status
        objective = solution.objective

        if status
            _set_stat!(param.stats, "MibSStatus", "Optimal")
            new_stat!(param.stats, "Opt", objective)
        else
            _set_stat!(param.stats, "MibSStatus", "NoOptimalSolution")
            @warn "The MibS solver did not find an optimal solution."
        end

        # process MibS output
        process_mibs_log!(param.stats)
    catch e
        error_message = sprint(showerror, e)
        _set_stat!(param.stats, "MibSStatus", "Error")
        _set_stat!(param.stats, "Error", error_message)
        @error "An error occurred while solving with MibS: $error_message"
    end
end

function process_mibs_log!(stats)
    mibs_output_file = joinpath(pwd(), "mibs_output.txt")
    if !isfile(mibs_output_file)
        @debug "MibS output file not found: $mibs_output_file"
        return
    end

    logtext = read(mibs_output_file, String)

    # extract relevant statistics
    match_nodes = match(r"Number of nodes processed:\s+(\d+)", logtext)
    match_runtime = match(r"Search wall-clock time:\s+([0-9.]+)", logtext)
    match_gap = match(r"Relative optimality gap is\s+([0-9.]+)%", logtext)

    nodes = isnothing(match_nodes) ? missing : parse(Int, match_nodes.captures[1])
    runtime = isnothing(match_runtime) ? missing : parse(Float64, match_runtime.captures[1])
    gap = isnothing(match_gap) ? missing : parse(Float64, match_gap.captures[1])

    new_stat!(stats, "BNodes", nodes)
    new_stat!(stats, "runtime", runtime)
    new_stat!(stats, "gap", gap)
end
