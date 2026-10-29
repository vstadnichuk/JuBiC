# A flexible implementation of a storage for run statistics and some function making working with them easier

using CSV, DataFrames

"""
Container for statistics collected during a JuBiC solve.

Solvers write status information, timings, bounds, objective values, and
experiment-specific metadata into the `data` dictionary. `solve_instance!`
returns this object.
"""
struct RunStats
    data::Dict{String,Any}
    RunStats() = new(Dict{String,Any}())  # Constructor to initialize an empty instance
end

"""
    new_stat!(stats::RunStats, name::String, value::Any)

Adds a new entry to this statistic. If a statistic with this name already exists, use _add_stat!_ instead. 
"""
function new_stat!(stats::RunStats, name::String, value::Any)
    if haskey(stats.data, name)
        throw(
            ArgumentError(
                "You tried to add a new statistic named $(name) that already exists. Please use _add_stat!_ instead!",
            ),
        )
    end
    stats.data[name] = value
end

"""
    add_stat!(stats::RunStats, name::String, value::Any)

Extends the value of an existing statistic by appending the new value with + opperator. 
"""
function add_stat!(stats::RunStats, name::String, value::Any)
    if !haskey(stats.data, name)
        throw(
            ArgumentError(
                "You tried to add a value $(value) to a statistic named $(name) that does not exists. Please use _new_stat!_ instead!",
            ),
        )
    end
    stats.data[name] += value
end

"""Initialize the disjoint ConnectorLP timing statistics."""
function _ensure_connector_time_stats!(stats::RunStats)
    for key in ("ConnectorLPTimeLP", "ConnectorLPTimePricing", "ConnectorLPTimePareto", "ConnectorLPTimeTotal")
        haskey(stats.data, key) || new_stat!(stats, key, 0.0)
    end
    return nothing
end

"""Add one disjoint ConnectorLP timing component and its total."""
function _add_connector_time!(stats::RunStats, key::String, value::Real)
    _ensure_connector_time_stats!(stats)
    add_stat!(stats, key, Float64(value))
    key == "ConnectorLPTimeTotal" || add_stat!(stats, "ConnectorLPTimeTotal", Float64(value))
    return nothing
end
