using JSON

include("hndp_experiment_runner.jl")

function _usage()
    return """
    Usage:
      julia --project=. examples/HNDP/run_hndp_manifest.jl <manifest_dir> <output_root> [--dry-run] [--no-resume]

    The manifest directory must contain:
      instances.json
      models.json
      params.json
    """
end

function _count_manifest_runs(manifest_dir::AbstractString)
    instance_path = joinpath(manifest_dir, "instances.json")
    model_path = joinpath(manifest_dir, "models.json")
    param_path = joinpath(manifest_dir, "params.json")

    ninstances = 0
    visit_hndp_networks(instance_path, _ -> begin
        ninstances += 1
    end)
    model_cfg = JSON.parsefile(model_path)
    param_cfg = JSON.parsefile(param_path)

    nmodels = length(get(model_cfg, "models", Any[]))
    nparams = length(get(param_cfg, "params", Any[]))
    return ninstances, nmodels, nparams, ninstances * nmodels * nparams
end

function main(args)
    length(args) >= 2 || error(_usage())

    manifest_dir = args[1]
    output_root = args[2]
    dry_run = "--dry-run" in args
    resume = !("--no-resume" in args)

    instance_path = joinpath(manifest_dir, "instances.json")
    model_path = joinpath(manifest_dir, "models.json")
    param_path = joinpath(manifest_dir, "params.json")

    for path in (instance_path, model_path, param_path)
        isfile(path) || error("Required manifest file not found: $(path)")
    end

    if dry_run
        ninstances, nmodels, nparams, total = _count_manifest_runs(manifest_dir)
        println("Manifest loads successfully.")
        println("Generated instances: $(ninstances)")
        println("Model specs: $(nmodels)")
        println("Parameter specs: $(nparams)")
        println("Total solve jobs: $(total)")
        return nothing
    end

    run_hndp_experiments!(
        instance_path,
        model_path,
        param_path;
        output_root=output_root,
        resume=resume,
    )
    return nothing
end

main(ARGS)
