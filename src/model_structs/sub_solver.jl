using JuMP

abstract type SubSolver end

struct SubSolution
    vio::Bool  # vio=true iff we found a new violated constraint
    obj_first_level::Number  # The first level obj. value (of the solution)
    obj_second_level::Number  # The second level obj. value appearing in first level (of the solution)
    A_sub::Any  # The set of resources that are contained in found solution. 
end

############ Custom error types for subsolvers (and functions related to subsolvers)  ############
struct TimeoutException <: Exception
    message::String
end
Base.showerror(io::IO, err::TimeoutException) = print(io, err.message)

############ Functions you have to implement yourself for your subsolver ############
"""
    capacity_linking(sub_solver::SubSolver, a, params::SolverParam)

Return the capacity within the interdiction constraint assoziated with resource a, 
i.e., if constraint has form y_a <= C_a x_a, with y_a second and x_a first level variables, the C_a values. 
"""
function capacity_linking(sub_solver::SubSolver, a, params::SolverParam)
    error(
        "You need to implement the function to obtain the capacities of your interdiction constraints for your own SubSolver!",
    )
end

"""
    check(sub_solver::SubSolver, params::SolverParam)

Is called at the before the start of optimization to run some (final) checks. 
"""
function check(sub_solver::SubSolver, params::SolverParam)
    error(
        "You need to implement this function to run the initial checks for your solver for your own SubSolver!",
    )
end

"""
    compute_lower_bound_master_contribution(sub_solver::SubSolver, params::SolverParam, time_limit)

Computes the minimal objective value that can be contributed to the master problem.
In other words, the lower bound to any solution w.r.t. the master objective function.

# Arguments
- 'sub_solver::SubSolver': The solver employed.
- 'params::SolverParam': The Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns
- 'lbm': The computed lower bound.
"""
function compute_lower_bound_master_contribution(sub_solver::SubSolver, params::SolverParam, time_limit)
    error(
        "You need to implement the function to obtain the minimal contribution to the master objective your solver can provide for your own SubSolver!",
    )
end

"""
    name(sub_solver::SubSolver)

Return the name (unique id) of your solver.
"""
function name(sub_solver::SubSolver)
    error(
        "You need to implement this function for your own SubSolver to return its name (unique identifier)!",
    )
end


"""
    separation!(sub_solver::SubSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)

Solves the sub_problem as separation problem for the Benders sub_problem in GBC generation.

# Arguments

- 'svals::SubSolver': The current value of the estimation for the objective term in master.
- 'gvals': The scaling factor for the cost solution.
- 'kvals': Mapping of resource to price.
- 'param::SolverParam': Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns
- 'sub_solver': The SubSolution found
"""
function separation!(sub_solver::SubSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)
    error("You need to implement the separation procedure for your own SubSolver!")
end

"""
    set_nthreads(sub_solver::SubSolver, n)

If multi thread is supported by your solver, please set the number of threads here. 
"""
function set_nthreads(sub_solver::SubSolver, n)
    printstyled(
        "Ignored the threads setting for subsolver $(name(sol)) because not supported";
        color=:orange,
    )
end

"""
    solve_sub_for_x(sub_solver::SubSolver, xvals, params::SolverParam, time_limit)

Solves the sub_problem for given solution of master variables.

# Arguments

- 'sub_solver::SubSolver': The solver employed.
- 'xvals': The values of the master linking variables for the current master solution.
- 'params::SolverParam': The Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns
- 'exists:Bool': True if solution was found. Otherwise, false. Note that time out causes a 'TimeoutException'.
- 'osol': The optimal objective value.
- 'y_sol': Is the solution of the second level y (interdiction) variables (as dict mapping resource to value).

If no solution exists, return false, 0, Dict()
Note that osol = 0 corresponds to removing the g variables from the objective. The subproblem should be solved to optimality if a solution exists.
"""
function solve_sub_for_x(sub_solver::SubSolver, xvals, params::SolverParam, time_limit)
    error(
        "You need to implement the function solving your second level problem for your own SubSolver!",
    )
end
