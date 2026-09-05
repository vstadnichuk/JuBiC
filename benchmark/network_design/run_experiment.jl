using JuBiC

include(joinpath(@__DIR__, "..", "..", "examples", "HNDP", "hndp_experiment_runner.jl"))

const ROOT = @__DIR__
const OUTPUT = get(ENV, "JUBIC_HNDP_OUTPUT", joinpath(ROOT, "results"))

run_hndp_experiments!(
    joinpath(ROOT, "instances.json"),
    joinpath(ROOT, "models.json"),
    joinpath(ROOT, "params.json");
    output_root=OUTPUT,
    resume=true,
)
