"""Run the BOBILib benchmark with MiBS and ND_WS_PD.

MiBS is run with its standard JuBiC wrapper. ND_WS_PD uses one master thread,
one parallel worker, and 15 Gurobi threads in that worker. Runs are executed
sequentially and appended to the same semicolon-separated CSV. Existing rows
for a solver/instance pair are reused.
"""

using JuBiC
using Gurobi
using CSV
using CodecZlib
using DataFrames

const REPO_ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const DEFAULT_INSTANCE_ROOT = joinpath(
    REPO_ROOT, "examples", "BOBILib", "instances", "benchmark-set",
    "benchmark-sets", "benchmark-interdiction",
)
const DEFAULT_OUTPUT = joinpath(@__DIR__, "results_bobilib_nd_ws_pd_sub15.csv")
const DEFAULT_SOURCE = ""
const TIME_LIMIT = 3600.0
const MIBS_LABEL = "MiBS"
const GBC_LABEL = "GBC_ND_WS_PD_SUB15"

function option(args, name, default)
    i = findfirst(==(name), args)
    isnothing(i) && return default
    i < length(args) || error("Missing value for $name")
    return args[i + 1]
end

function materialize_mps(path, work_dir)
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

function append_stats!(output, stats)
    row = JuBiC.stats_to_dataframe([stats])
    if isfile(output)
        old = CSV.read(output, DataFrame; delim=';', decimal=',')
        CSV.write(output, vcat(old, row; cols=:union); delim=';', decimal=',')
    else
        CSV.write(output, row; delim=';', decimal=',')
    end
end

function add_metadata!(stats, name, solver, mps, aux, root)
    new_stat!(stats, "benchmark_instance", name)
    new_stat!(stats, "benchmark_solver", solver)
    new_stat!(stats, "benchmark_mps", relpath(mps, root))
    new_stat!(stats, "benchmark_aux", relpath(aux, root))
    new_stat!(stats, "benchmark_time_limit", TIME_LIMIT)
end

function source_instance_names(source)
    isempty(source) && return nothing
    isfile(source) || error("Source CSV does not exist: $source")
    old = CSV.read(source, DataFrame; delim=';', decimal=',')
    :benchmark_instance in propertynames(old) || error("Source CSV lacks benchmark_instance")
    return Set(String.(old.benchmark_instance))
end

function discover_pairs(instance_root)
    pairs = Tuple{String,String,String}[]
    for (root, _, files) in walkdir(instance_root)
        basename(root) == "PaxHeader" && continue
        for file in files
            endswith(file, ".aux") || continue
            stem = replace(file, ".aux" => "")
            plain = joinpath(root, stem * ".mps")
            compressed = joinpath(root, stem * ".mps.gz")
            mps = isfile(plain) ? plain : compressed
            isfile(mps) && push!(pairs, (stem, mps, joinpath(root, file)))
        end
    end
    sort!(pairs; by=first)
    return pairs
end

function solve_mibs(name, mps, aux, work_dir, root)
    stats = JuBiC.RunStats()
    input_dir = joinpath(work_dir, "inputs", "MiBS")
    log_dir = joinpath(work_dir, "logs", "MiBS", name)
    mkpath(input_dir); mkpath(log_dir)
    local_mps = materialize_mps(mps, input_dir)
    instance = get_MibS_instance(local_mps, aux; stats=stats)
    params = MibSparam(false, log_dir, TIME_LIMIT, stats)
    add_metadata!(stats, name, MIBS_LABEL, mps, aux, root)
    solve_instance!(instance, params)
    return stats
end

function solve_nd_ws_pd(name, mps, aux, work_dir, root)
    stats = JuBiC.RunStats()
    input_dir = joinpath(work_dir, "inputs", "ND_WS_PD")
    log_dir = joinpath(work_dir, "logs", "ND_WS_PD", name)
    mkpath(input_dir); mkpath(log_dir)
    local_mps = materialize_mps(mps, input_dir)
    instance = get_GBC_instance(local_mps, aux, Gurobi.Optimizer;
        partial_decomposition=true, preprocessing=true, stats=stats)
    # One master thread, one worker, and 15 Gurobi threads in that worker.
    # ConnectorLP uses one thread in this configuration.
    params = GBCparam(GurobiSolver(), false, log_dir, "lp", stats,
        TIME_LIMIT, 42, 1, 15, true, PARETO_OPTIMALITY_ONLY,
        true, false, true, 1e9, 0, true; parallel_subsolvers=1)
    add_metadata!(stats, name, GBC_LABEL, mps, aux, root)
    solve_instance!(instance, params)
    return stats
end

function failed_stats(name, solver, mps, aux, root, err)
    stats = JuBiC.RunStats()
    add_metadata!(stats, name, solver, mps, aux, root)
    status = err isa JuBiC.UnsupportedInstanceFormatError ? "UnsupportedFormat" : "Error"
    new_stat!(stats, "RunStatus", status)
    new_stat!(stats, "GBCStatus", status)
    new_stat!(stats, "Error", sprint(showerror, err))
    return stats
end

function main()
    instance_root = normpath(option(ARGS, "--interdiction-root", DEFAULT_INSTANCE_ROOT))
    output = normpath(option(ARGS, "--output", DEFAULT_OUTPUT))
    source_arg = option(ARGS, "--source-csv", DEFAULT_SOURCE)
    isdir(instance_root) || error("Interdiction instance directory does not exist: $instance_root")
    mkpath(dirname(output))

    existing = isfile(output) ? CSV.read(output, DataFrame; delim=';', decimal=',') : DataFrame()
    completed = Set{Tuple{String,String}}()
    if :benchmark_solver in propertynames(existing) && :benchmark_instance in propertynames(existing)
        for row in eachrow(existing)
            push!(completed, (String(row.benchmark_solver), String(row.benchmark_instance)))
        end
    end

    pairs = discover_pairs(instance_root)
    requested = source_instance_names(source_arg)
    isnothing(requested) || filter!(p -> p[1] in requested, pairs)
    isempty(pairs) && error("No benchmark instances were found below $instance_root")

    work_root = mktempdir()
    for (index, (name, source_mps, aux_path)) in enumerate(pairs)
        startswith(lowercase(name), "general") && continue
        println("BOBILib [$index/$(length(pairs))] $name")
        for solver in (MIBS_LABEL, GBC_LABEL)
            (solver, name) in completed && continue
            stats = try
                solver == MIBS_LABEL ?
                    solve_mibs(name, source_mps, aux_path, work_root, instance_root) :
                    solve_nd_ws_pd(name, source_mps, aux_path, work_root, instance_root)
            catch err
                @error "Failed on $name with $solver" exception=(err, catch_backtrace())
                failed_stats(name, solver, source_mps, aux_path, instance_root, err)
            end
            append_stats!(output, stats)
            push!(completed, (solver, name))
        end
    end
end

main()
