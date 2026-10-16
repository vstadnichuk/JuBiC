"""Reproduce the published EMA BlC/BlCLag benchmark.

Run from the JuBiC repository root with:
    julia --project=. benchmark_run/reproduce_ema_big_runs.jl

The raw run output is written below benchmark_run/runs/ unless
JUBIC_REPRO_OUTPUT is set.  BlCLag is reproduced with warm starts only.
"""

ENV["JUBIC_TOPOLOGY"] = get(ENV, "JUBIC_TOPOLOGY", "ema")
ENV["JUBIC_USERS"] = get(ENV, "JUBIC_USERS", "50,75,100,125,150,175,200")
ENV["JUBIC_DECISION_ARCS"] = get(ENV, "JUBIC_DECISION_ARCS", "43")
ENV["JUBIC_SINGLE_LAYER_MODES"] = get(ENV, "JUBIC_SINGLE_LAYER_MODES", "decision_only_negative")
ENV["JUBIC_THREADS"] = get(ENV, "JUBIC_THREADS", "8")
ENV["JUBIC_BRANCHING_RULE"] = get(ENV, "JUBIC_BRANCHING_RULE", "fixed_linking_order")
ENV["JUBIC_GUROBI_MINIMAL"] = get(ENV, "JUBIC_GUROBI_MINIMAL", "0")
ENV["JUBIC_RUN_BLC"] = get(ENV, "JUBIC_RUN_BLC", "1")
ENV["JUBIC_RUN_BLCLAG"] = get(ENV, "JUBIC_RUN_BLCLAG", "1")
ENV["JUBIC_BLCLAG_WARM_ONLY"] = get(ENV, "JUBIC_BLCLAG_WARM_ONLY", "1")

include(joinpath(@__DIR__, "..", "examples", "HNDP", "reproduce_fixed_order_big_runs.jl"))
