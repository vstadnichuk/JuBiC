using CSV
using DataFrames
using Statistics

include("hndp_experiment_runner.jl")

"""
    run_hndp_calibration!(instance_cfg, model_cfg, param_cfg; output_root, resume=true)

Run the HNDP experiment pipeline and return the resulting summary table.
"""
function run_hndp_calibration!(
    instance_cfg::AbstractString,
    model_cfg::AbstractString,
    param_cfg::AbstractString;
    output_root::AbstractString,
    resume::Bool=true,
)
    run_hndp_experiments!(
        instance_cfg,
        model_cfg,
        param_cfg;
        output_root=output_root,
        resume=resume,
    )
    return CSV.read(joinpath(output_root, "results", "batch_summary.csv"), DataFrame)
end

"""
    summarize_hndp_calibration(df)

Print a compact runtime summary grouped by instance and model.
"""
function summarize_hndp_calibration(df::DataFrame)
    if !("instance_name" in names(df) && "model_name" in names(df) && "runtime" in names(df))
        error("The calibration table must contain at least the columns instance_name, model_name, and runtime.")
    end

    selected = select(
        df,
        :instance_name,
        :model_name,
        :runtime,
        :Opt,
        :Opt_status,
        :gap,
        :instance_nusers,
        :instance_length_slack,
    )
    sort!(selected, [:instance_name, :runtime], rev=[false, true])
    return selected
end
