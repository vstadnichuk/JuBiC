using JSON

include("hndp_experiment_runner.jl")

function _usage()
    return """
    Usage:
      julia --project=. examples/HNDP/run_hndp_manifest.jl <config_file_or_manifest_dir> <output_root> [--dry-run] [--no-resume]

    Supply one JSON file containing instances, models, and params, or a
    directory containing:
      instances.json
      models.json
      params.json
    """
end

function _manifest_paths(config_path::AbstractString)
    isfile(config_path) && return (config_path, config_path, config_path)
    return (joinpath(config_path, "instances.json"),
            joinpath(config_path, "models.json"),
            joinpath(config_path, "params.json"))
end

function _count_manifest_runs(config_path::AbstractString)
    instance_path, model_path, param_path = _manifest_paths(config_path)

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

    config_path = args[1]
    output_root = args[2]
    dry_run = "--dry-run" in args
    resume = !("--no-resume" in args)

    instance_path, model_path, param_path = _manifest_paths(config_path)

    for path in (instance_path, model_path, param_path)
        isfile(path) || error("Required manifest file not found: $(path)")
    end

    if dry_run
        ninstances, nmodels, nparams, total = _count_manifest_runs(config_path)
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

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
