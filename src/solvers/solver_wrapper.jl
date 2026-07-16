# a wrapper for the solvers we support internally
using JuMP, Gurobi
using Base.Threads


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
    silent::Bool
    precreated_envs::Vector{Any}
    env_lock::ReentrantLock
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

function _create_gurobi_env(silent::Bool)
    return silent ? silent_gurobi_env() : Gurobi.Env()
end

"""
    GurobiSolver(; silent=true)

Convenience constructor for a Gurobi-based solver wrapper. By default we create
the underlying `Gurobi.Env` in silent mode to avoid confusing initialization
messages for end users.
"""
function GurobiSolver(; silent::Bool=true, precreate_envs::Integer=0)
    solver = GurobiSolver(silent, Any[], ReentrantLock())
    reserve_environments!(solver, precreate_envs)
    return solver
end

"""
    reserve_environments!(solver::SolverWrapper, n)

Prepare `n` solver-side environments in advance. Solver wrappers that do not
manage explicit environment objects can implement this as a no-op.
"""
function reserve_environments!(solver::SolverWrapper, n::Integer)
    return nothing
end

function reserve_environments!(solver::GurobiSolver, n::Integer)
    n <= 0 && return nothing
    threadid() == 1 || throw(
        ArgumentError(
            "Gurobi environments must be created on the main Julia thread. Requested $(n) additional environments from worker thread $(threadid()).",
        ),
    )
    lock(solver.env_lock) do
        for _ in 1:n
            push!(solver.precreated_envs, _create_gurobi_env(solver.silent))
        end
    end
    return nothing
end



"""
    get_next_optimizer(s::SolverWrapper)

Generate the next optimizer object for a model using your solver. 
"""
function get_next_optimizer(s::SolverWrapper)
    print("You need to overload this function for your solver variant!")
end

function _take_gurobi_env!(solver::GurobiSolver)
    env = nothing
    lock(solver.env_lock) do
        if !isempty(solver.precreated_envs)
            env = pop!(solver.precreated_envs)
        end
    end
    if !isnothing(env)
        return env
    end

    threadid() == 1 || throw(
        ArgumentError(
            "JuBiC tried to create a new Gurobi environment on worker thread $(threadid()). " *
            "Gurobi environment creation is not thread-safe. Pre-create environments on the main thread " *
            "and avoid constructing new Gurobi-backed JuMP models inside threaded code.",
        ),
    )
    return _create_gurobi_env(solver.silent)
end

function get_next_optimizer(s::GurobiSolver)
    # Each optimizer receives an exclusive environment. JuBiC intentionally does
    # not share one Gurobi.Env across concurrently solved models.
    return Gurobi.Optimizer(_take_gurobi_env!(s))
end

"""
    _run_post_gurobi_cleanup!()

Force Julia to run finalizers for unreachable Gurobi-backed JuMP objects at a
controlled point between benchmark cases. This avoids deferred model/env
teardown from happening later inside another active Gurobi callback, which has
caused `GRBfree` access violations in long sequential benchmark runs.
"""
function _run_post_gurobi_cleanup!()
    GC.gc()
    GC.gc()
    return nothing
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
