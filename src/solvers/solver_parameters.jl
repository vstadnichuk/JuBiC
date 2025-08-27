# Way to pass parameters within solution process within the implemented solvers
# TODO: I do not think that having an own file for all parameters is a good idea, or?
# TODO: parameters implemented not clean

abstract type SolverParam end

"""
    get_stats(param::SolverParam)

Return the 'RunStats' object where we should save the statistics.
"""
function get_stats(param::SolverParam)
    error("Please implement the function 'get_stats' for your own parameter settings.")
end

"""
    output_file_path(param::SolverParam)

Return the path to the output file (as String).
"""
function output_file_path(param::SolverParam)
    error(
        "Please implement the function 'output_file_path' for your own parameter settings.",
    )
end

"""
    should_debbug_print(param::SolverParam)

Return true if debbug text should be generated/printed.
"""
function should_debbug_print(param::SolverParam)
    error(
        "Please implement the function 'should_debbug_print' for your own parameter settings.",
    )
end

###### Parameters for JuMP based solver ######

@enum ParetoCut begin
    PARETO_NONE  # No Pareto cuts used
    PARETO_OPTIMALITY_ONLY  # Pareto cuts used only for optimality cuts
    PARETO_OPTIMALITY_AND_FEASIBILITY  # Pareto cuts used for both optimality and feasibility
end

struct GBCparam <: SolverParam
    solver::SolverWrapper
    debbug_out::Bool  # if true, print additional output to files (debug messages and files during run of model)
    output_folder_path::Any  # path to where to store output files 
    file_format_output::Any  # the ending for the file format to write models to file, e.g."lp" or "mps". 
    stats::RunStats  # Store statistics of the run here

    runtime::Number  # the maximal runtime of the master MIP (in seconds)
    threads_master::Integer # number of threads used in the master MIP problem
    threads_sub_con::Any  # number of threads used for LP solver in ConnectorLP. Suggested number of threads for sub_problem solver (but depends on solver if supported)

    pareto::ParetoCut  # setting for pareto optimality cuts

    infinity_num::Any  # Number used in subroblems to add sufisticated lower and upper bounds. Set it to some positiv value that can be considered infinity in your problem
end

GBCparam(solver, debbug_out, output_folder_path, file_format_output) = GBCparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), 3600, 8, 8, PARETO_OPTIMALITY_ONLY, 1e9)
GBCparam(solver, debbug_out, output_folder_path, file_format_output, pareto) = GBCparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), 3600, 8, 8, pareto, 1e9)
GBCparam(solver, debbug_out, output_folder_path, file_format_output, pareto, runtime) = GBCparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), runtime, 8, 8, pareto, 1e9)

function get_stats(param::GBCparam)
    return param.stats
end

function output_file_path(param::GBCparam)
    return param.output_folder_path
end

function should_debbug_print(param::GBCparam)
    return param.debbug_out
end


###### Parameter for Benders-like Cuts Solver (BlCSolver) ######
struct BLCparam <: SolverParam
    solver::SolverWrapper
    debbug_out::Bool  # if true, print additional output to files (debug messages and files during run of model)
    output_folder_path::Any  # path to where to store output files 
    file_format_output::Any  # the ending for the file format to write models to file, e.g."lp" or "mps". 
    stats::RunStats  # Store statistics of the run here

    runtime::Any  # the maximal runtime of the master MIP (in seconds)
    threads_master::Any  # number of threads used in the master MIP problem
    threads_sub_con::Any  # number of threads used for LP solver in ConnectorLP. Suggested number of threads for sub_problem solver (but depends on solver if supported)
end

BLCparam(solver, debbug_out, output_folder_path, file_format_output, runtime) = BLCparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), runtime, 8, 8)
BLCparam(solver, debbug_out, output_folder_path, file_format_output) = BLCparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), 3600, 8, 8)


function get_stats(param::BLCparam)
    return param.stats
end

function output_file_path(param::BLCparam)
    return param.output_folder_path
end

function should_debbug_print(param::BLCparam)
    return param.debbug_out
end


###### Parameter for the compact Solver (MIPSolver) that is a wrapper for underlying MIP solver ######
struct MIPparam <: SolverParam
    solver::SolverWrapper
    debbug_out::Bool  # if true, print additional output to files (debug messages and files during run of model)
    output_folder_path::Any  # path to where to store output files 
    file_format_output::Any  # the ending for the file format to write models to file, e.g."lp" or "mps". 
    stats::RunStats  # Store statistics of the run here

    runtime::Any  # the maximal runtime of the master MIP (in seconds)
    threads_master::Any  # number of threads used in the master MIP problem
end

MIPparam(solver, debbug_out, output_folder_path, file_format_output, runtime) = MIPparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), runtime, 8)
MIPparam(solver, debbug_out, output_folder_path, file_format_output) = MIPparam(solver, debbug_out, output_folder_path, file_format_output, RunStats(), 3600, 8)


function get_stats(param::MIPparam)
    return param.stats
end

function output_file_path(param::MIPparam)
    return param.output_folder_path
end

function should_debbug_print(param::MIPparam)
    return param.debbug_out
end


###### Parameter for MibSSolver that is a wrapper for the MibS solver ######
struct MibSparam <: SolverParam
    debbug_out::Bool  # if true, print additional output to files (debug messages and files during run of model)
    output_folder_path::Any  # path to where to store output files 
    stats::RunStats  # Store statistics of the run here
end

function get_stats(param::MibSparam)
    return param.stats
end

function output_file_path(param::MibSparam)
    return param.output_folder_path
end

function should_debbug_print(param::MibSparam)
    return param.debbug_out
end

