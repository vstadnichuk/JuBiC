using JuBiC

include(joinpath(@__DIR__, "..", "..", "examples", "HNDP", "hndp_experiment_runner.jl"))

const ROOT = @__DIR__
const OUTPUT = get(ENV, "JUBIC_HNDP_OUTPUT", joinpath(ROOT, "results"))

function _network_design_model_filter(generated_network, model_spec)
    model_name = String(get(model_spec, "name", ""))
    if model_name != "ND_WS_ASTAR"
        return true
    end

    # The A* implementation requires nonnegative costs in every active
    # objective state. For the HNDP competition instances, nonnegative risk
    # and routing costs imply nonnegative master, follower, and connector
    # arc costs (the generated linking penalties are nonnegative).
    hndp = generated_network.instance
    return all(
        isfinite(Float64(user.mrisk[src(edge), dst(edge)])) &&
        isfinite(Float64(user.mcost[src(edge), dst(edge)])) &&
        user.mrisk[src(edge), dst(edge)] >= 0 &&
        user.mcost[src(edge), dst(edge)] >= 0
        for user in hndp.users for edge in edges(hndp.mygraph)
    )
end

run_hndp_experiments!(
    joinpath(ROOT, "instances.json"),
    joinpath(ROOT, "models.json"),
    joinpath(ROOT, "params.json");
    output_root=OUTPUT,
    resume=true,
    model_filter=_network_design_model_filter,
)
