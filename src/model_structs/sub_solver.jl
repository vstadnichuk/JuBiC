using JuMP

abstract type SubSolver end

"""
Result of a connector pricing/separation call.

`obj_compare` is the objective of the returned feasible pricing solution and is
therefore an upper bound for the minimization pricing problem. `obj_bound` is a
certified lower bound for the same problem. Exact subsolvers return equal
values. Heuristic implementations that cannot provide a certificate must return
`-Inf` as `obj_bound`; JuBiC will then report an infinite approximation error.
"""
struct SubSolution
    vio::Bool  # vio=true iff we found a new violated constraint
    obj_first_level::Number  # The first level obj. value (of the solution)
    obj_second_level::Number  # The second level obj. value appearing in first level (of the solution)
    obj_compare::Number  # feasible pricing objective (upper bound for minimization)
    obj_bound::Number  # certified pricing lower bound; -Inf means unavailable
    A_sub::Any  # The set of resources that are contained in found solution. 
end

# Existing exact subsolvers need not repeat their objective as a bound.
SubSolution(vio, obj_first_level, obj_second_level, obj_compare, A_sub) =
    SubSolution(vio, obj_first_level, obj_second_level, obj_compare, obj_compare, A_sub)

############ Custom error types for subsolvers (and functions related to subsolvers)  ############
struct TimeoutException <: Exception
    message::String
end
Base.showerror(io::IO, err::TimeoutException) = print(io, err.message)

struct NumericalIssueException <: Exception
    message::String
    status::String
end
Base.showerror(io::IO, err::NumericalIssueException) = print(io, err.message)

struct MibSFailureException <: Exception
    message::String
end
Base.showerror(io::IO, err::MibSFailureException) = print(io, err.message)

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
This initialization routine must be solved to proven optimality; returning only
a heuristic incumbent is invalid because the value is used as a global lower
bound in every connector model. Heuristic behavior is permitted only in
`separation!`, where it must be accompanied by a certified objective bound.

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
When implementing the function, you can assume that 'gvals' and 'kvals' are non negative. 

# Arguments

- 'svals::SubSolver': The current value of the estimation for the objective term in master.
- 'gvals': The scaling factor for the cost solution.
- 'kvals': Mapping of resource to price.
- 'param::SolverParam': Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns
- `SubSolution`: `obj_compare` is the feasible pricing objective (an upper
  bound for minimization), and `obj_bound` is a certified lower bound. Exact
  implementations may use the five-argument constructor; inexact
  implementations without a finite certificate must pass `-Inf` explicitly.
"""
function separation!(sub_solver::SubSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)
    error("You need to implement the separation procedure for your own SubSolver!")
end

"""
    separation_exact!(sub_solver, sval, gvals, kvals, params, time_limit)

Run connector pricing to proven optimality. `ConnectorLP` uses this method for
feasibility-cut separation, where an approximate ray is not sufficient.
The default delegates to `separation!`. A subsolver whose ordinary pricing may
stop before optimality can overload this method to provide the exact pricing
needed for feasibility-cut generation. Whether ordinary pricing was exact is
determined from the returned `SubSolution` bounds, not from a subsolver type or
mode flag.
"""
function separation_exact!(sub_solver::SubSolver, sval, gvals, kvals::Dict, params::SolverParam, time_limit)
    return separation!(sub_solver, sval, gvals, kvals, params, time_limit)
end


"""
    separation_BlC!(sub_solver::SubSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)

Solves the separation problem for the Benders-like sub_problem in BlC generation. 
This involves finding the feasible solution that maximizes 'obj_second_level - sum kvals'.
Attention: You need to ensure that the found second-level solution is bilevel feasible, as otherwise the behaviour is undefined.
When implementing the function, you can assume that 'kvals' are non negative. 

# Arguments

- 'svals::SubSolver': The current value of the estimation for the objective term in master.
- 'kvals': Mapping of resource to price.
- 'param::SolverParam': Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns
- 'sub_solver': The SubSolution found
"""
function separation_BlC!(sub_solver::SubSolver, sval, kvals::Dict, param::SolverParam, time_limit)
    error("You need to implement the separation procedure for your own SubSolver!")
end


"""
    set_nthreads(sub_solver::SubSolver, n)

If multi thread is supported by your solver, please set the number of threads here. 
"""
function set_nthreads(sub_solver::SubSolver, n)
    printstyled(
        "Ignored the threads setting for subsolver $(name(sub_solver)) because not supported";
        color=:orange,
    )
end

"""
    set_singlethread(sub_solver::SubSolver)

Force the subsolver to use a single thread. Solvers that do not support
multi-threading can implement this as a no-op.
"""
function set_singlethread(sub_solver::SubSolver)
    set_nthreads(sub_solver, 1)
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
- 'osol': The optimal objective value (w.r.t. second-level function)
- 'osol_L1': The objective value (w.r.t. the first-level function) for found solution
- 'y_sol': Is the solution of the second level y (interdiction) variables (as dict mapping resource to value).

If no solution exists, return false, 0, Dict()
Note that osol = 0 corresponds to removing the g variables from the objective. The subproblem should be solved to optimality if a solution exists.
"""
function solve_sub_for_x(sub_solver::SubSolver, xvals, params::SolverParam, time_limit)
    error(
        "You need to implement the function solving your second level problem for your own SubSolver!",
    )
end

"""
    solve_sub_for_x_optimistic(sub_solver, xvals, params, time_limit)

Optionally solve the follower problem for fixed first-level values using
optimistic bilevel tie-breaking. The method first minimizes the follower
objective and then minimizes the follower contribution to the first-level
objective over the complete set of follower-optimal solutions.

The return tuple is identical to `solve_sub_for_x`: `(exists, follower_value,
first_level_contribution, y_solution)`. This is an optional interface: the
generic function intentionally has no fallback method. GBC checks
applicability only during final incumbent evaluation and falls back to the
ordinary cached fixed-`x` response when a subsolver does not implement it.
"""
function solve_sub_for_x_optimistic end

"""
    supports_bilevel_subproblem_solver(sub_solver::SubSolver)

Return whether the subsolver can solve the bilevel subproblems required by the
BlCLag solver. By default, custom subsolvers are assumed not to support this.
"""
function supports_bilevel_subproblem_solver(sub_solver::SubSolver)
    return false
end
