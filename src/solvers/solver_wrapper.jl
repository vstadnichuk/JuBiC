# a wrapper for the solvers we support internally
using JuMP, Gurobi


abstract type SolverWrapper end

struct GurobiSolver <: SolverWrapper
    env::Any  # The enviroment for the currernt solver 
end



"""
    get_next_optimizer(s::SolverWrapper)

Generate the next optimizer object for a model using your solver. 
"""
function get_next_optimizer(s::SolverWrapper)
    print("You need to overload this function for your solver variant!")
end


function get_next_optimizer(s::GurobiSolver)
    # TODO Warning: Gurobi.Env are NOT thread-safe. If two models both use the same environment you must not solve them simultaneously on different threads.
    return Gurobi.Optimizer(s.env)
end
