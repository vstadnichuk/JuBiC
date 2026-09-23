using JuMP

abstract type SubSolver end

struct SubSolution
    vio::Bool  # vio=true iff we found a new violated constraint
    obj_first_level::Number  # The first level obj. value (of the solution)
    obj_second_level::Number  # The second level obj. value appearing in first level (of the solution)
    obj_compare::Number  # The full subproblem objective value used for the violation test
    A_sub::Any  # The set of resources that are contained in found solution. 
end

############ Custom error types for subsolvers (and functions related to subsolvers)  ############
struct TimeoutException <: Exception
    message::String
end
Base.showerror(io::IO, err::TimeoutException) = print(io, err.message)

struct NumericalIssueException <: Exception
    message::String
    status::String
    context::Dict{String,Any}
end

"""
    _round_y_values_and_objectives(sol, y_vals, optL2, optL1)

Convert numerically non-binary follower linking variables to binary values and
re-evaluate both objective expressions with those rounded values.  This keeps
the returned follower solution and objective values consistent when a solver
returns a tiny fractional value for a binary variable.
"""
function _round_y_values_and_objectives(sol, y_vals, optL2, optL1)
    rounded = Dict(a => Float64(y_vals[a]) for a in keys(y_vals))
    steps = String[]
    for a in keys(rounded)
        raw = rounded[a]
        if abs(raw) <= 1e-8
            rounded[a] = 0.0
        elseif abs(raw - 1.0) <= 1e-8
            rounded[a] = 1.0
        else
            rounded_value = min(max(Float64(ceil(raw)), 0.0), 1.0)
            rounded[a] = rounded_value
            push!(steps, "y[$(a)]=$(raw) -> $(rounded_value)")
        end
    end

    isempty(steps) && return y_vals, optL2, optL1

    @warn "Subsolver $(name(sol)) returned non-binary follower y-values. Rounded values: $(join(steps, "; ")). Re-evaluating both follower objective expressions on the rounded y-values."

    y_by_variable = Dict(sol.y_vars[a] => rounded[a] for a in keys(rounded))
    evaluate_with_rounded_y(expr) = JuMP.constant(expr) + sum(
        coefficient * get(y_by_variable, variable, Float64(value(variable)))
        for (coefficient, variable) in JuMP.linear_terms(expr);
        init=0.0,
    )

    rounded_optL2 = evaluate_with_rounded_y(sol.c_objterm)
    rounded_optL1 = evaluate_with_rounded_y(sol.r_objterm)
    @debug "Subsolver $(name(sol)) changed follower objective values after y-rounding: optL2 $(optL2) -> $(rounded_optL2), optL1 $(optL1) -> $(rounded_optL1)."
    return rounded, rounded_optL2, rounded_optL1
end
NumericalIssueException(message::String, status::String) =
    NumericalIssueException(message, status, Dict{String,Any}())
Base.showerror(io::IO, err::NumericalIssueException) = print(io, err.message)

"""
    _round_integer_objective_if_close(value, parameter, sub_name, description)

When `integer_obj=true`, round an objective-derived value to its nearest
integer only if it is already within the numerical tolerance.  Objective
values that are not close to an integer are left unchanged; integer-objective
mode must not turn a small, genuinely fractional value into one unit.
"""
const INTEGER_OBJECTIVE_ROUND_TOL = 1e-4

function _round_integer_objective_if_close(
    value::Real,
    parameter::SolverParam,
    sub_name,
    description::AbstractString,
)
    value_float = Float64(value)
    integer_obj = hasproperty(parameter, :integer_obj) && getproperty(parameter, :integer_obj)
    integer_obj || return value_float

    integer_value = round(value_float)
    if abs(value_float - integer_value) <= INTEGER_OBJECTIVE_ROUND_TOL
        if value_float != integer_value
            @debug "$(sub_name) integerizes $(description): $(value_float) -> $(integer_value)."
        end
        return Float64(integer_value)
    end
    return value_float
end

struct MibSFailureException <: Exception
    message::String
end
Base.showerror(io::IO, err::MibSFailureException) = print(io, err.message)

"""
    _flush_model_updates!(model)

Flush pending model modifications after a batch of temporary constraints has
been deleted. JuMP normally flushes these changes at the next `optimize!`, but
the Gurobi wrapper maintains additional bookkeeping for deleted columns and
constraints. Calling its internal update routine here keeps that bookkeeping
and Gurobi's native model synchronized before the persistent model is reused.

For non-Gurobi models, or models without an attached optimizer, this is a
no-op.
"""
function _flush_model_updates!(model)
    try
        optimizer = JuMP.unsafe_backend(model)
        if optimizer isa Gurobi.Optimizer
            # Use Gurobi.jl's update path rather than calling GRBupdatemodel
            # directly: the wrapper also adjusts its internal indices after
            # deletions and resets its pending-change flags.
            Gurobi._update_if_necessary(optimizer; force=true)
        end
    catch err
        # Some auxiliary model types (notably BilevelJuMP models) do not expose
        # a JuMP optimizer backend directly. They are flushed by their own
        # solver interface, so leave those models untouched here.
        if err isa MethodError || err isa UndefRefError
            return nothing
        end
        rethrow()
    end
    return nothing
end

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
When implementing the function, you can assume that 'gvals' and 'kvals' are non negative. 

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
        "Ignored the threads setting for subsolver $(name(sol)) because not supported";
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
    supports_bilevel_subproblem_solver(sub_solver::SubSolver)

Return whether the subsolver can solve the bilevel subproblems required by the
BlCLag solver. By default, custom subsolvers are assumed not to support this.
"""
function supports_bilevel_subproblem_solver(sub_solver::SubSolver)
    return false
end
