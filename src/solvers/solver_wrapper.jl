# a wrapper for the solvers we support internally
using JuMP, Gurobi


abstract type SolverWrapper end

struct GurobiSolver <: SolverWrapper
    env::Any  # The enviroment for the currernt solver 
end

"""
    silent_gurobi_env()

Create a Gurobi environment with console output disabled. This suppresses the
license and initialization lines emitted by `Gurobi.Env()` itself; model-level
`set_silent(...)` only applies after the environment already exists.
"""
function silent_gurobi_env()
    return Gurobi.Env(Dict{String,Any}("OutputFlag" => 0))
end

"""
    GurobiSolver(; silent=true)

Convenience constructor for a Gurobi-based solver wrapper. By default we create
the underlying `Gurobi.Env` in silent mode to avoid confusing initialization
messages for end users.
"""
function GurobiSolver(; silent::Bool=true)
    env = silent ? silent_gurobi_env() : Gurobi.Env()
    return GurobiSolver(env)
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

"""
    set_seed!(model, solver, seed)

Apply a solver-specific random seed to the passed JuMP model. If the wrapped
solver does not support explicit seeding, this function is a no-op.
"""
function set_seed!(model::JuMP.Model, s::SolverWrapper, seed)
    return nothing
end

function set_seed!(model::JuMP.Model, s::GurobiSolver, seed)
    set_optimizer_attribute(model, "Seed", seed)
    return nothing
end
