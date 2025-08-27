using Base.Filesystem, Logging, LoggingExtras


"""
    append_constraintlist_to_file(list, filepath)

Append to the passed file the list of JuMP model constraints.
It should also print other expression object nicely. 
"""
function append_constraintlist_to_file(list, filepath)
    open(filepath, "a") do f
        for c in list
            println(f, c)
        end
    end
end

function key_master_sol(msol::Dict, A)
    # generate a tuple from the current master solution that can be used as key for dict
    alist = []
    for a in A
        x = msol[a]
        nval = (x == -zero(x)) ? zero(x) : x  # normalize to avoid negative 0
        append!(alist, nval)
    end
    return Tuple(alist)
end

"""
    new_file_logger(filepath, print_debug)

Creates a new logger that writes to a file where the path is given by 'filepath'. If 'print_debug' is true, the logger also prints the @debug messages. 
"""
function new_file_logger(filepath, print_debug)
    if print_debug
        logger = MinLevelLogger(FileLogger(filepath), Logging.Debug)
        @info("Registered a new Logger for debuging. Please note that the debug messages are only printed to the logfile and not to the console. ")
    else
        logger = MinLevelLogger(FileLogger(filepath), Logging.Info)
        @info("Registered a new Logger but it will not print any debug messages as it was not set in the parameters of the solvers.")
    end
    return logger
end

"""
    print_collected_cuts(param::SolverParam, sol_to_lazy::Dict)

Print the passed list of JuMP constraints 'sol_to_lazy' to file.
"""
function print_collected_cuts(param::SolverParam, sol_to_lazy::Dict)
    filepath = param.output_folder_path * "/mastercuts_collection.txt"
    for lazylist in values(sol_to_lazy)
        append_constraintlist_to_file(lazylist, filepath)
    end
end

"""
    print_solution_to_file(mobj, xvars, params::SolverParam)

Print the passed solution value 'mobj' and a mapping of variables to value ('xvars') to file.
"""
function print_solution_to_file(mobj, xvars, params::SolverParam)
    outfilecut = params.output_folder_path * "/master_sol.txt"
    open(outfilecut, "w") do f
        println(f, "The objective value is $(mobj)")
        println(f, xvars)
    end
end


function round_master_solution(msol::Dict)
    # round values in the master solution mapping to 0 or 1. Creates new dict that is returned
    # TODO: nach Auxiliaries.jl verschieben?
    nsol = Dict{Any,Int}()

    for (a, sol) in msol
        if !(sol >= -0.000001 && sol <= 1.00001)
            error("The value of linking variable for $a is not binary but $sol")
        end
        nsol[a] = (sol > 0.5) ? 1 : 0
    end

    return nsol
end


function set_optimization_status_stats(status, param)
    if status == MOI.OPTIMAL
        new_stat!(param.stats, "Opt_status", "Optimal")
    elseif status == MOI.LOCALLY_SOLVED
        new_stat!(param.stats, "Opt_status", "Suboptimal")
    elseif status == MOI.TIME_LIMIT
        new_stat!(param.stats, "Opt_status", "Timelimit")
    else
        new_stat!(param.stats, "Opt_status", "Infeasible")
    end
end
