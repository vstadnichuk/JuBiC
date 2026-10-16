"""Disk-disciplined randomized sensitivity sweep for multimodal HNDP instances.

Only `results/batch_summary.csv` is written below the requested output root.
"""

using JuBiC
using Random

include("hndp_experiment_runner.jl")

const SWEEP_TOPOLOGIES = ["sioux_falls", "ema", "berlin_mitte_center"]
const SWEEP_SCENARIOS = ["multimodal_bike", "multimodal_expansion"]

_pick(rng, values) = values[rand(rng, 1:length(values))]

function _make_sweep_specs(; samples_per_cell::Int=20, seed::Int=20260819)
    rng = MersenneTwister(seed)
    specs = Dict{String,Any}[]
    for topology in SWEEP_TOPOLOGIES, scenario in SWEEP_SCENARIOS
        for _ in 1:samples_per_cell
            station_fraction = topology == "berlin_mitte_center" ?
                _pick(rng, [10, 20, 30, 40]) :
                _pick(rng, [10, 20, 40, 60, 80])
            delta_min = _pick(rng, [0, 1000, 2000, 3000, 4000])
            delta_max = min(10_000, delta_min + _pick(rng, [1000, 2000, 3000, 4000, 5000]))
            spec = Dict{String,Any}(
                "name" => "sweep_$(topology)_$(scenario)",
                "instance_type" => scenario,
                "topologies" => [topology],
                "nusers" => [_pick(rng, [100, 250, 500, 1000, 2000])],
                "parameter_seeds" => [seed + length(specs) + 1],
                "station_fraction_percent" => station_fraction,
                "expansion_station_fraction_percent" => _pick(rng, [5, 10, 20, 30]),
                "station_budget_fraction_percent" => _pick(rng, [10, 25, 50, 75, 100]),
                "arc_budget_fraction_percent" => _pick(rng, [5, 10, 25, 50, 75]),
                "expansion_station_count" => _pick(rng, [1, 2, 3, 4, 6, 8]),
                "car_speed_kmh" => _pick(rng, [25, 35, 45, 55]),
                "bike_speed_kmh" => _pick(rng, [10, 15, 20, 25]),
                "car_degree_penalty_bp" => _pick(rng, [0, 1000, 2000, 3000]),
                "car_capacity_boost_bp" => _pick(rng, [0, 1000, 2000]),
                "bike_low_capacity_boost_bp" => _pick(rng, [0, 1000, 2000]),
                "transit_delta_min_bp" => delta_min,
                "transit_delta_max_bp" => delta_max,
                "bike_station_profit" => _pick(rng, [1, 2, 4]),
                "bike_station_cost" => _pick(rng, [0, 1, 2, 4]),
                "expansion_arc_profit" => _pick(rng, [1, 2, 4]),
                "expansion_arc_cost" => _pick(rng, [0, 1, 2, 4]),
            )
            push!(specs, spec)
        end
    end
    return specs
end

function _sweep_models()
    return [
        Dict{String,Any}("name" => "sd_fixed_network_path", "model_type" => "sd", "big_m_mode" => "fixed_network_path", "indicator_constraints" => false, "bound_duals" => true),
        Dict{String,Any}("name" => "path_all_accelerations", "model_type" => "path", "enumeration_time_limit" => 600.0, "parallelize" => true, "use_decision_arc_dominance" => true),
    ]
end

function _append_sweep_row!(summary_path::AbstractString, row::Dict{String,Any})
    JuBiC._append_batch_summary_csv!(summary_path, row)
    # The general batch helper keeps recovery copies.  This sweep is designed
    # to retain only the aggregate CSV, so remove those transient copies after
    # each successful append.
    result_dir = dirname(summary_path)
    for entry in readdir(result_dir; join=true)
        if endswith(entry, ".bak") || occursin("__snapshot_", basename(entry))
            rm(entry; force=true)
        end
    end
end

function run_multimodal_sweep!(output_root::AbstractString; samples_per_cell::Int=20, seed::Int=20260819)
    result_dir = joinpath(output_root, "results")
    mkpath(result_dir)
    summary_path = joinpath(result_dir, "batch_summary.csv")
    specs = _make_sweep_specs(; samples_per_cell=samples_per_cell, seed=seed)
    models = _sweep_models()
    params = Dict{String,Any}("name" => "gurobi_10min", "mip_solver" => "Gurobi", "runtime" => 600.0, "threads_master" => 8, "threads_sub_con" => 8, "parallel_separation" => true, "seed" => seed)
    total = length(specs) * length(models)
    println("Starting multimodal sweep: $(length(specs)) instances × $(length(models)) models = $(total) jobs")
    completed = 0
    started = time()
    for spec in specs
        generated = only(generate_hndp_networks(Dict{String,Any}("instances" => [spec])))
        for model_spec in models
            experiment_id = _hndp_experiment_id(generated.name, model_spec, params)
            aggregate_row = Dict{String,Any}()
            try
                _, aggregate_row = _run_hndp_experiment(generated, model_spec, params, output_root, false)
            catch err
                aggregate_row = _build_hndp_error_row(generated, model_spec, params, experiment_id, err)
            end
            aggregate_row["experiment_id"] = experiment_id
            aggregate_row["instance_name"] = generated.name
            _append_sweep_row!(summary_path, aggregate_row)
            completed += 1
            println("[$completed/$total] $(generated.name) / $(model_spec["name"]) ($(round((time() - started) / 3600; digits=2)) h)")
        end
    end
    println("Sweep complete. Results: $(summary_path)")
    return summary_path
end

length(ARGS) >= 1 || error("Usage: julia --project=. examples/HNDP/run_multimodal_sweep.jl <output_root> [samples_per_cell]")
output_root = ARGS[1]
samples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 20
run_multimodal_sweep!(output_root; samples_per_cell=samples)
