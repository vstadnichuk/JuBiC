using JuBiC
using CSV
using DataFrames
using Dates
using Statistics

include("hndp_model_generation.jl")

const SWEEP_RUNTIME_LIMIT = 600.0
const SWEEP_PRICING_GAPS = [0.01, 0.05, 0.20]
const SWEEP_MODES = [CONNECTOR_UNDERESTIMATION, CONNECTOR_OVERESTIMATION]

function _sweep_instance_config()
    balanced = Dict("max_cost" => 30, "max_risk" => 30, "max_weight" => 30,
        "construction_cost_min" => 2, "construction_cost_max" => 10)
    risk_heavy = Dict("max_cost" => 12, "max_risk" => 100, "max_weight" => 25,
        "construction_cost_min" => 10, "construction_cost_max" => 30)
    cost_heavy = Dict("max_cost" => 100, "max_risk" => 15, "max_weight" => 40,
        "construction_cost_min" => 1, "construction_cost_max" => 5)

    function csp(name, topology, profile; constrained, user_mode="shared", nusers=2)
        return merge(
            Dict{String,Any}(
                "name" => name,
                "instance_type" => "constrained_shortest_path",
                "topologies" => [topology],
                "nusers" => [nusers],
                "length_constrained" => [constrained],
                "length_slack" => [0.25],
                "two_stage" => [false],
                "user_parameter_mode" => [user_mode],
                "availability_budget_fraction" => [0.35],
            ),
            profile,
        )
    end

    function competition(name, topology, profile, beta; constrained, user_mode="shared", nusers=2)
        return merge(
            Dict{String,Any}(
                "name" => name,
                "instance_type" => "competition",
                "topologies" => [topology],
                "nusers" => [nusers],
                "length_constrained" => [constrained],
                "length_slack" => [0.25],
                "competitor_cost_factor" => [beta],
                "user_parameter_mode" => [user_mode],
                "availability_budget_fraction" => [0.35],
            ),
            profile,
        )
    end

    return Dict{String,Any}(
        "parameter_seeds" => [7],
        "instances" => [
            csp("sioux_balanced", "sioux_falls", balanced; constrained=false),
            csp("sioux_risk", "sioux_falls", risk_heavy; constrained=true, user_mode="per_user"),
            csp("anaheim_balanced", "anaheim", balanced; constrained=false, nusers=1),
            csp("anaheim_cost", "anaheim", cost_heavy; constrained=true, user_mode="per_user", nusers=1),
            csp("ema_risk", "ema", risk_heavy; constrained=false, nusers=1),
            csp("friedrich_cost", "friedrichshain_center", cost_heavy; constrained=true, nusers=1),
            competition("layered_sioux_balanced", "layered_sioux_falls", balanced, 0.75; constrained=false),
            competition("layered_ema_risk", "layered_ema", risk_heavy, 1.25; constrained=true, nusers=1),
        ],
    )
end

function _sweep_param(solver, output_dir, mode)
    return GBCparam(
        solver, false, output_dir, "lp", JuBiC.RunStats(), SWEEP_RUNTIME_LIMIT,
        42, 4, 1, false, PARETO_NONE, true, false, true, 1e9, 0,
        false, 1e-4, 1e-4, true, true, mode,
    )
end

_stat(stats, key, default=missing) = get(stats.data, key, default)
_finite_number(value) = value isa Real && isfinite(value)

function _run_sweep_gbc(generated, output_dir; approximate, mip_gap=nothing,
    mode=CONNECTOR_UNDERESTIMATION)
    solver = GurobiSolver()
    instance = build_hndp_gbc_instance(
        generated.instance,
        solver;
        partial_decomposition=true,
        subproblem_method=HNDP_SUBPROBLEM_MIP,
        heuristic_subsolver=approximate,
        heuristic_mip_gap=mip_gap,
    )
    params = _sweep_param(solver, output_dir, mode)
    params.stats.data["enable_output_logs"] = false
    try
        return solve_instance!(instance, params)
    finally
        instance = nothing
        params = nothing
        solver = nothing
        JuBiC._run_post_gurobi_cleanup!()
    end
end

function _base_row(generated)
    metadata = generated.metadata
    return Dict{String,Any}(
        "instance" => generated.name,
        "instance_type" => metadata["instance_type"],
        "topology" => metadata["topology_family"],
        "nusers" => metadata["nusers"],
        "nnodes" => metadata["nnodes"],
        "narcs" => metadata["narcs"],
        "decision_arcs" => metadata["decision_arcs"],
        "length_constrained" => metadata["length_constrained"],
        "max_cost" => metadata["max_cost"],
        "max_risk" => metadata["max_risk"],
        "max_weight" => metadata["max_weight"],
        "construction_cost_min" => metadata["construction_cost_min"],
        "construction_cost_max" => metadata["construction_cost_max"],
    )
end

function _stats_row(generated, stats; run_kind, mode, mip_gap, exact_opt=missing)
    row = _base_row(generated)
    merge!(row, Dict{String,Any}(
        "run_kind" => run_kind,
        "mode" => mode,
        "requested_mip_gap" => mip_gap,
        "exact_opt" => exact_opt,
        "opt" => _stat(stats, "Opt"),
        "opt_status" => _stat(stats, "Opt_status"),
        "gbc_status" => _stat(stats, "GBCStatus"),
        "gbc_result_status" => _stat(stats, "GBCResultStatus"),
        "solution_type" => _stat(stats, "GBCSolutionType"),
        "runtime" => _stat(stats, "runtime"),
        "lower" => _stat(stats, "ObjectiveIntervalLower"),
        "upper" => _stat(stats, "ObjectiveIntervalUpper"),
        "interval_width" => _stat(stats, "ObjectiveIntervalWidth"),
        "interval_certified" => _stat(stats, "ObjectiveIntervalCertified", false),
        "interval_width_bound" => _stat(stats, "ObjectiveIntervalWidthBound"),
        "width_bounded_by_pricing" => _stat(
            stats,
            "ObjectiveIntervalWidthBoundedByPricingError",
            false,
        ),
        "optimistic_evaluation_complete" => _stat(
            stats,
            "FinalOptimisticEvaluationComplete",
            false,
        ),
        "optimistic_evaluation_time" => _stat(stats, "FinalOptimisticEvaluationTime"),
        "master_bound" => _stat(stats, "MasterObjectiveBound"),
        "master_gap" => _stat(stats, "MasterAbsoluteGap"),
        "connector_error" => _stat(stats, "ConnectorApproximationErrorBound"),
        "final_connector_error" => _stat(stats, "FinalConnectorErrorSum"),
        "pricing_max_error" => _stat(stats, "ConnectorPricingMaxError"),
        "inexact_pricing_calls" => _stat(stats, "NInexactPricingCalls", 0),
        "n_opt_cuts" => _stat(stats, "NOptCuts", 0),
        "n_feas_cuts" => _stat(stats, "NFeasCuts", 0),
        "interval_contains_opt" => missing,
        "error" => missing,
    ))
    return row
end

function _error_row(generated, err; run_kind, mode, mip_gap, exact_opt=missing)
    row = _base_row(generated)
    merge!(row, Dict{String,Any}(
        "run_kind" => run_kind, "mode" => mode,
        "requested_mip_gap" => mip_gap, "exact_opt" => exact_opt,
        "opt" => missing, "opt_status" => "Error", "gbc_status" => "Error",
        "gbc_result_status" => "Error", "solution_type" => missing,
        "runtime" => missing, "lower" => missing, "upper" => missing,
        "interval_width" => missing, "interval_certified" => false,
        "interval_width_bound" => missing, "width_bounded_by_pricing" => false,
        "optimistic_evaluation_complete" => false,
        "optimistic_evaluation_time" => missing,
        "master_bound" => missing, "master_gap" => missing,
        "connector_error" => missing, "final_connector_error" => missing,
        "pricing_max_error" => missing,
        "inexact_pricing_calls" => missing, "n_opt_cuts" => missing,
        "n_feas_cuts" => missing, "interval_contains_opt" => missing,
        "error" => sprint(showerror, err),
    ))
    return row
end

function _write_sweep_rows(path, rows)
    CSV.write(path, DataFrame(rows))
end

function _interval_contains(row, optimum)
    lower, upper = row["lower"], row["upper"]
    lower isa Real && upper isa Real || return false
    (isnan(lower) || isnan(upper)) && return false
    tolerance = 1e-5 * max(1.0, abs(Float64(optimum)))
    return lower <= optimum + tolerance && optimum <= upper + tolerance
end

function _write_sweep_report(path, rows)
    exact_rows = [r for r in rows if r["run_kind"] == "exact"]
    accepted = [r for r in exact_rows if r["exact_opt"] !== missing]
    approximate = [r for r in rows if r["run_kind"] == "approximation"]
    interval_failures = count(r -> r["interval_contains_opt"] === false, approximate)
    optimistic_failures = count(r -> r["optimistic_evaluation_complete"] !== true, approximate)
    width_bound_failures = count(
        r -> r["width_bounded_by_pricing"] !== true,
        approximate,
    )

    open(path, "w") do io
        println(io, "# HNDP certified connector approximation sweep")
        println(io)
        println(io, "Generated: $(Dates.now())")
        println(io)
        println(io, "- Generated instances: $(length(exact_rows))")
        println(io, "- Exact optima proven within 600 seconds: $(length(accepted))")
        println(io, "- Instances excluded without a proven exact optimum: $(length(exact_rows) - length(accepted))")
        println(io, "- Approximation runs completed: $(length(approximate))")
        println(io, "- Certified interval failures: $(interval_failures)")
        println(io, "- Incomplete optimistic final evaluations: $(optimistic_failures)")
        println(io, "- Missing pricing-error width bounds: $(width_bound_failures)")
        println(io)
        println(io, "| mode | requested gap | runs | heuristic results | mean runtime (s) | median runtime (s) | mean width | mean lower-bound shortfall |")
        println(io, "|---|---:|---:|---:|---:|---:|---:|---:|")
        for mode in ("underestimation", "overestimation")
            for gap in SWEEP_PRICING_GAPS
                group = [r for r in approximate if r["mode"] == mode && r["requested_mip_gap"] == gap && r["error"] === missing]
                runtimes = Float64[r["runtime"] for r in group if _finite_number(r["runtime"])]
                widths = Float64[r["interval_width"] for r in group if _finite_number(r["interval_width"])]
                shortfalls = Float64[r["exact_opt"] - r["lower"] for r in group if _finite_number(r["exact_opt"]) && _finite_number(r["lower"])]
                heuristic_count = count(r -> r["solution_type"] == "Heuristic", group)
                mean_runtime = isempty(runtimes) ? NaN : mean(runtimes)
                median_runtime = isempty(runtimes) ? NaN : median(runtimes)
                mean_width = isempty(widths) ? NaN : mean(widths)
                mean_shortfall = isempty(shortfalls) ? NaN : mean(shortfalls)
                println(io, "| $(mode) | $(gap) | $(length(group)) | $(heuristic_count) | $(round(mean_runtime; digits=3)) | $(round(median_runtime; digits=3)) | $(round(mean_width; digits=4)) | $(round(mean_shortfall; digits=4)) |")
            end
        end
    end
end

function run_hndp_certified_approximation_sweep(output_root::AbstractString)
    mkpath(output_root)
    result_path = joinpath(output_root, "sweep_results.csv")
    report_path = joinpath(output_root, "sweep_report.md")
    rows = Dict{String,Any}[]
    generated_instances = generate_hndp_networks(_sweep_instance_config())

    for generated in generated_instances
        println("[exact] ", generated.name)
        exact_stats = try
            _run_sweep_gbc(
                generated,
                joinpath(output_root, "runs", generated.name, "exact");
                approximate=false,
            )
        catch err
            push!(rows, _error_row(generated, err; run_kind="exact", mode="exact", mip_gap=0.0))
            _write_sweep_rows(result_path, rows)
            continue
        end

        exact_row = _stats_row(generated, exact_stats; run_kind="exact", mode="exact", mip_gap=0.0)
        exact_is_proven = exact_row["opt_status"] == "Optimal" && _finite_number(exact_row["opt"])
        exact_opt = exact_is_proven ? Float64(exact_row["opt"]) : missing
        exact_row["exact_opt"] = exact_opt
        exact_row["interval_contains_opt"] = exact_is_proven ? _interval_contains(exact_row, exact_opt) : missing
        push!(rows, exact_row)
        _write_sweep_rows(result_path, rows)
        exact_is_proven || continue

        for gap in SWEEP_PRICING_GAPS
            for mode in SWEEP_MODES
                mode_name = string(mode)
                println("[approx] ", generated.name, " mode=", mode_name, " gap=", gap)
                row = try
                    stats = _run_sweep_gbc(
                        generated,
                        joinpath(output_root, "runs", generated.name, "$(mode_name)_$(gap)");
                        approximate=true,
                        mip_gap=gap,
                        mode=mode,
                    )
                    _stats_row(
                        generated,
                        stats;
                        run_kind="approximation",
                        mode=mode_name,
                        mip_gap=gap,
                        exact_opt=exact_opt,
                    )
                catch err
                    _error_row(
                        generated,
                        err;
                        run_kind="approximation",
                        mode=mode_name,
                        mip_gap=gap,
                        exact_opt=exact_opt,
                    )
                end

                if row["error"] === missing
                    row["interval_contains_opt"] = _interval_contains(row, exact_opt)
                end
                push!(rows, row)
                _write_sweep_rows(result_path, rows)

                if row["interval_certified"] !== true || row["interval_contains_opt"] === false
                    _write_sweep_report(report_path, rows)
                    failed_lower = row["lower"]
                    failed_upper = row["upper"]
                    error(
                        "Missing or invalid certified interval for $(generated.name), mode=$(mode_name), " *
                        "gap=$(gap): exact optimum $(exact_opt) is not in " *
                        "[$(failed_lower), $(failed_upper)]. Results saved to $(result_path).",
                    )
                end
                if row["width_bounded_by_pricing"] !== true ||
                   row["optimistic_evaluation_complete"] !== true
                    _write_sweep_report(report_path, rows)
                    error(
                        "Missing optimistic pricing-error width certificate for " *
                        "$(generated.name), mode=$(mode_name), gap=$(gap).",
                    )
                end
                observed_width = row["interval_width"]
                certified_width_bound = row["interval_width_bound"]
                width_tolerance = 1e-5 * max(1.0, abs(Float64(certified_width_bound)))
                if observed_width > certified_width_bound + width_tolerance
                    _write_sweep_report(report_path, rows)
                    error(
                        "Observed interval width $(observed_width) exceeds its " *
                        "pricing-error bound $(certified_width_bound) for " *
                        "$(generated.name), mode=$(mode_name), gap=$(gap).",
                    )
                end
            end
        end
    end

    _write_sweep_report(report_path, rows)
    println("Sweep results: ", result_path)
    println("Sweep report: ", report_path)
    return rows
end

if abspath(PROGRAM_FILE) == @__FILE__
    default_output = joinpath(
        JuBiC.repo_local_tempdir("benchmarks", "hndp_certified_approximation"),
        Dates.format(Dates.now(), "yyyymmdd_HHMMSS"),
    )
    output_root = isempty(ARGS) ? default_output : abspath(ARGS[1])
    run_hndp_certified_approximation_sweep(output_root)
end
