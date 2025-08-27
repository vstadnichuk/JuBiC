# Just a wrapper for the underlying MIP solver
using JSON

function solve_with_MIP!(inst::Instance, param::MIPparam)
    mipm::MIPMaster = inst.master

    # register solver
    @debug "Starting setup of the MIP solver."
    new_stat!(param.stats, "Solver", "MIPSolver")
    new_stat!(param.stats, "time_limit", param.runtime)

    # write model to lp file and set path to log file
    try
        write_to_file(mipm.mymip, param.output_folder_path * "/MIPinstance.$(param.file_format_output)")
        set_optimizer_attribute(mipm.mymip, "LogFile", param.output_folder_path *"/mip_log.txt")
    catch err 
        @error "Could not write model to file or set log file for folder $(param.output_folder_path). Error is $err"
    end

    # set runtime and number of threads
    set_time_limit_sec(mipm.mymip, param.runtime)
    set_attribute(mipm.mymip, MOI.NumberOfThreads(), param.threads_master)

    # add callback to master and solve 
    @debug "Finished model construction. Now proceeding to optimization process with MIP solver."
    optimize!(mipm.mymip)

    # print solution and collected data
    if termination_status(mipm.mymip) == MOI.OPTIMAL || termination_status(mipm.mymip) == MOI.LOCALLY_SOLVED || termination_status(mipm.mymip) == MOI.TIME_LIMIT
        if primal_status(mipm.mymip) == MOI.FEASIBLE_POINT
            # found a (suboptimal) solution (if time limit existence of solution not ensured)
            mobj = objective_value(mipm.mymip)
            @debug "The MIP models objective is $(mobj)."
            new_stat!(param.stats, "Opt", mobj)

            # print solution to file
            solution = Dict(JuMP.name(x) => JuMP.value(x) for x in all_variables(mipm.mymip))
            write(param.output_folder_path*"/solution.json", JSON.json(solution))
        end

        # set status 
        status = termination_status(mipm.mymip)
        set_optimization_status_stats(status, param)
    else
        @debug "The MIP is infeasible with status: $(termination_status(mipm.mymip))"
        set_optimization_status_stats(termination_status(mipm.mymip), param)
    end

    # save run data to statistics
    new_stat!(param.stats, "runtime", solve_time(mipm.mymip))
    new_stat!(param.stats, "gap", JuMP.relative_gap(mipm.mymip))
    new_stat!(param.stats, "BNodes", MOI.get(mipm.mymip, MOI.NodeCount()))
end