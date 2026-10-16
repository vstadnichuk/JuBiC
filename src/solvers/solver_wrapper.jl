# a wrapper for the solvers we support internally
using JuMP, Gurobi


"""
Abstract wrapper for MIP solvers used by JuBiC.

Concrete wrappers must provide `get_next_optimizer` so JuBiC can create fresh
JuMP optimizer objects for master models, connector LPs, and helper models.
"""
abstract type SolverWrapper end

"""
Gurobi-based implementation of `SolverWrapper`.

The wrapper stores a `Gurobi.Env` and creates `Gurobi.Optimizer` objects through
`get_next_optimizer`.
"""
mutable struct GurobiSolver <: SolverWrapper
    # The environment used by the master/default model factory.
    env::Any
    # Worker environments are created before threaded separation starts.  A
    # model assigned to a worker must keep using that worker's environment.
    worker_envs::Vector{Any}
    silent::Bool
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
    return GurobiSolver(env, Any[], silent)
end

# Preserve the historical one-argument constructor used by client code.
GurobiSolver(env) = GurobiSolver(env, Any[], false)

"""Create `n` worker environments on the calling (initialization) thread."""
function ensure_worker_envs!(s::GurobiSolver, n::Integer)
    n >= 1 || throw(ArgumentError("The number of Gurobi worker environments must be positive."))
    while length(s.worker_envs) < n
        env = s.silent ? silent_gurobi_env() : Gurobi.Env()
        push!(s.worker_envs, env)
    end
    return s.worker_envs
end

"""Return an optimizer bound to the requested pre-created worker environment."""
function get_worker_optimizer(s::GurobiSolver, worker_id::Integer)
    1 <= worker_id <= length(s.worker_envs) ||
        throw(ArgumentError("Worker environment $(worker_id) has not been initialized."))
    return Gurobi.Optimizer(s.worker_envs[worker_id])
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
