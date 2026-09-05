using JuBiC
using Gurobi
using CSV
using CodecZlib

const ROOT = @__DIR__
const INSTANCE_ROOT = normpath(joinpath(ROOT, "..", "..", "examples", "BOBILib", "instances", "benchmark-set", "benchmark-sets"))
const OUTPUT = joinpath(ROOT, "results.csv")
const WORK_ROOT = joinpath(ROOT, "working_instances")
const TIME_LIMIT = 3600.0
const MAX_INSTANCES = parse(Int, get(ENV, "JUBIC_BOBILIB_MAX_INSTANCES", "0"))

function materialize_mps(path::String, work_dir::String)
    endswith(path, ".mps") && return path
    endswith(path, ".mps.gz") || error("Unsupported MPS file: $path")
    target = joinpath(work_dir, splitext(basename(path))[1])
    open(path, "r") do input
        stream = GzipDecompressorStream(input)
        open(target, "w") do output
            write(output, read(stream))
        end
        close(stream)
    end
    return target
end

function add_metadata!(stats, instance_name, solver_name, mps_path, aux_path)
    new_stat!(stats, "benchmark_instance", instance_name)
    new_stat!(stats, "benchmark_solver", solver_name)
    new_stat!(stats, "benchmark_mps", relpath(mps_path, INSTANCE_ROOT))
    new_stat!(stats, "benchmark_aux", relpath(aux_path, INSTANCE_ROOT))
    new_stat!(stats, "benchmark_time_limit", TIME_LIMIT)
    return stats
end

function solve_gbc(mps_path, aux_path, instance_name, output_dir)
    stats = JuBiC.RunStats()
    instance = get_GBC_instance(mps_path, aux_path, Gurobi.Optimizer;
        partial_decomposition=true, preprocessing=true, stats=stats)
    params = GBCparam(
        GurobiSolver(), false, output_dir, "lp", stats, TIME_LIMIT,
        42, 8, 1, true, PARETO_OPTIMALITY_ONLY, true, false, true,
        1e9, 0, false,
    )
    add_metadata!(stats, instance_name, "GBC_WS_PD", mps_path, aux_path)
    solve_instance!(instance, params)
    return stats
end

function solve_mibs(mps_path, aux_path, instance_name, output_dir)
    stats = JuBiC.RunStats()
    instance = get_MibS_instance(mps_path, aux_path; stats=stats)
    params = MibSparam(false, output_dir, TIME_LIMIT, stats)
    add_metadata!(stats, instance_name, "MiBS", mps_path, aux_path)
    solve_instance!(instance, params)
    return stats
end

function main()
    mkpath(WORK_ROOT)
    stats_list = JuBiC.RunStats[]
    pairs = Tuple{String,String,String}[]
    for (root, _, files) in walkdir(INSTANCE_ROOT)
        for file in files
            startswith(file, ".") && continue
            endswith(file, ".aux") || continue
            aux_path = joinpath(root, file)
            stem = replace(file, ".aux" => "")
            plain = joinpath(root, stem * ".mps")
            compressed = joinpath(root, stem * ".mps.gz")
            mps_path = isfile(plain) ? plain : compressed
            isfile(mps_path) || continue
            push!(pairs, (stem, mps_path, aux_path))
        end
    end
    sort!(pairs; by=first)
    if MAX_INSTANCES > 0
        pairs = first(pairs, min(MAX_INSTANCES, length(pairs)))
    end

    for (index, (name, source_mps, aux_path)) in enumerate(pairs)
        println("[$index/$(length(pairs))] $name")
        mps_path = materialize_mps(source_mps, WORK_ROOT)
        instance_dir = joinpath(WORK_ROOT, name)
        mkpath(instance_dir)
        for solver in ("GBC", "MiBS")
            try
                stats = solver == "GBC" ?
                    solve_gbc(mps_path, aux_path, name, instance_dir) :
                    solve_mibs(mps_path, aux_path, name, instance_dir)
                push!(stats_list, stats)
            catch err
                stats = JuBiC.RunStats()
                add_metadata!(stats, name, solver, source_mps, aux_path)
                new_stat!(stats, "RunStatus", "Error")
                new_stat!(stats, "Error", sprint(showerror, err))
                push!(stats_list, stats)
                @error "Failed on $name with $solver" exception=(err, catch_backtrace())
            end
            CSV.write(OUTPUT, JuBiC.stats_to_dataframe(stats_list); delim=';', decimal=',')
        end
    end
    println("Wrote $(length(stats_list)) rows to $OUTPUT")
end

main()
