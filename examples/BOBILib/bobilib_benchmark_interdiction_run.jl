using Dates

const BENCHMARK_INTERDICTION_ROOT = joinpath(
    @__DIR__,
    "instances",
    "benchmark-set",
    "benchmark-sets",
    "benchmark-interdiction",
)

if !haskey(ENV, "BOBILIB_OUTPUT_ROOT")
    ENV["BOBILIB_OUTPUT_ROOT"] = joinpath(
        "tmp_compare",
        "runs",
        "bobilib_benchmark_interdiction_clean_$(Dates.format(now(), "yyyymmdd_HHMMSS"))",
    )
end

include(joinpath(@__DIR__, "bobilib_knapsack_benchmark.jl"))

function _parse_resume_csv_line(line::AbstractString)
    values = String[]
    buffer = IOBuffer()
    in_quotes = false
    i = firstindex(line)
    while i <= lastindex(line)
        c = line[i]
        if c == '"'
            if in_quotes && i < lastindex(line) && line[nextind(line, i)] == '"'
                print(buffer, '"')
                i = nextind(line, i)
            else
                in_quotes = !in_quotes
            end
        elseif c == ',' && !in_quotes
            push!(values, String(take!(buffer)))
        else
            print(buffer, c)
        end
        i = nextind(line, i)
    end
    push!(values, String(take!(buffer)))
    return values
end

function _read_resume_csv_records(path::AbstractString)
    content = read(path, String)
    records = String[]
    buffer = IOBuffer()
    in_quotes = false
    i = firstindex(content)
    while i <= lastindex(content)
        c = content[i]
        if c == '"'
            print(buffer, c)
            if in_quotes && i < lastindex(content) && content[nextind(content, i)] == '"'
                i = nextind(content, i)
                print(buffer, content[i])
            else
                in_quotes = !in_quotes
            end
        elseif (c == '\n' || c == '\r') && !in_quotes
            record = String(take!(buffer))
            isempty(strip(record)) || push!(records, record)
            if c == '\r' && i < lastindex(content) && content[nextind(content, i)] == '\n'
                i = nextind(content, i)
            end
        else
            print(buffer, c)
        end
        i = nextind(content, i)
    end
    record = String(take!(buffer))
    isempty(strip(record)) || push!(records, record)
    return records
end

function _load_existing_results!()
    isempty(RESULT_ROWS) || return length(RESULT_ROWS)
    isfile(RESULT_CSV) || return 0

    records = _read_resume_csv_records(RESULT_CSV)
    isempty(records) && return 0
    header = _parse_resume_csv_line(first(records))
    for record in Iterators.drop(records, 1)
        values = _parse_resume_csv_line(record)
        row = Dict{String,Any}()
        for (idx, key) in enumerate(header)
            row[key] = idx <= length(values) ? values[idx] : ""
        end
        if isempty(get(row, "experiment", "")) || isempty(get(row, "instance", "")) || isempty(get(row, "solver_config", ""))
            @warn "Skipping malformed resume row without complete experiment/instance/solver_config metadata."
            continue
        end
        push!(RESULT_ROWS, row)
    end
    return length(RESULT_ROWS)
end

function _has_existing_result(instance::BOBIInstance, variant::String, solver_config::String)
    return any(
        row -> get(row, "experiment", "") == variant &&
               get(row, "instance", "") == instance.name &&
               get(row, "solver_config", "") == solver_config,
        RESULT_ROWS,
    )
end

function _discover_benchmark_interdiction_instances(root::AbstractString=BENCHMARK_INTERDICTION_ROOT)
    instances = BOBIInstance[]
    for file in sort!(readdir(root))
        endswith(lowercase(file), ".aux") || continue
        startswith(file, "._") && continue
        aux_path = joinpath(root, file)
        mps_name = _mps_name_from_aux(aux_path)
        mps_name === nothing && continue
        startswith(mps_name, "._") && continue
        mps_path = joinpath(root, mps_name)
        if !isfile(mps_path) && isfile(mps_path * ".gz")
            mps_path *= ".gz"
        end
        isfile(mps_path) || error("Missing MPS file for AUX $(aux_path): expected $(mps_name) or $(mps_name).gz")
        name = splitext(file)[1]
        push!(instances, BOBIInstance("benchmark-interdiction", -1, name, mps_path, aux_path))
    end
    return instances
end

function _has_upper_level_coupling_constraints(mps_path::String, aux_path::String)
    mps_data = JuBiC._read_mps(_plain_mps_path(mps_path))
    aux_data = JuBiC._read_aux(aux_path)
    return JuBiC._has_upper_level_coupling_constraints(mps_data, aux_data)
end

function _solve_gbc_with_connector_options(
    mps_path,
    aux_path;
    partial_decomposition::Bool,
    out_dir::String,
    connector_add_current_solution_cut::Bool,
    subsolver_numerical_preprocessing::Bool,
)
    mps_path = _plain_mps_path(mps_path)
    stats = JuBiC.RunStats()
    inst = get_GBC_instance(mps_path, aux_path, Gurobi.Optimizer; partial_decomposition=partial_decomposition, preprocessing=true, stats=stats)
    params = GBCparam(
        GurobiSolver(),
        false,
        out_dir,
        "lp",
        stats,
        TIME_LIMIT,
        42,
        8,
        8,
        true,
        PARETO_OPTIMALITY_ONLY,
        true,
        false,
        true,
        1e9,
        0,
        false,
        1e-4,
        1e-4,
        connector_add_current_solution_cut,
        subsolver_numerical_preprocessing,
    )
    new_stat!(params.stats, "enable_output_logs", false)
    new_stat!(params.stats, "connector_add_current_solution_cut", connector_add_current_solution_cut)
    new_stat!(params.stats, "subsolver_numerical_preprocessing", subsolver_numerical_preprocessing)
    new_stat!(params.stats, "gbc_partial_decomposition", partial_decomposition)
    new_stat!(params.stats, "gbc_warmstart", true)
    new_stat!(params.stats, "gbc_pareto_mode", "PARETO_OPTIMALITY_ONLY")
    new_stat!(params.stats, "gbc_bigMwithLC", false)
    new_stat!(params.stats, "gbc_trim_coeff", true)
    return solve_instance!(inst, params)
end

function _base_row_meta(
    instance::BOBIInstance,
    variant::String,
    solver_config::String,
    mps_path::String,
    aux_path::String;
    connector_add_current_solution_cut=nothing,
    subsolver_numerical_preprocessing=nothing,
    partial_decomposition=nothing,
)
    return Dict{String,Any}(
        "experiment" => variant,
        "solver_config" => solver_config,
        "family" => instance.family,
        "item_count" => instance.item_count,
        "instance" => instance.name,
        "source_mps" => mps_path,
        "source_aux" => aux_path,
        "runtime_limit" => TIME_LIMIT,
        "has_upper_level_coupling_constraints" => false,
        "connector_add_current_solution_cut" => connector_add_current_solution_cut,
        "subsolver_numerical_preprocessing" => subsolver_numerical_preprocessing,
        "gbc_partial_decomposition" => partial_decomposition,
        "gbc_warmstart" => solver_config == "mibs" ? nothing : true,
        "gbc_pareto_mode" => solver_config == "mibs" ? nothing : "PARETO_OPTIMALITY_ONLY",
        "gbc_bigMwithLC" => solver_config == "mibs" ? nothing : false,
        "gbc_trim_coeff" => solver_config == "mibs" ? nothing : true,
    )
end

function _append_exception_row!(row_meta::Dict{String,Any}, err)
    row_meta["RunStatus"] = "Exception"
    row_meta["Error"] = sprint(showerror, err)
    _append_csv!(RESULT_CSV, row_meta)
end

function _run_single_clean_benchmark(instance::BOBIInstance, variant::String, mps_path::String, aux_path::String)
    run_dir = joinpath(OUTPUT_ROOT, "runs", variant, instance.family, instance.name, "mibs")
    row_meta = _base_row_meta(instance, variant, "mibs", mps_path, aux_path)
    if _has_existing_result(instance, variant, "mibs")
        println("  skipping existing mibs")
    else
        try
            stats = _solve_mibs(mps_path, aux_path; out_dir=run_dir)
            _append_csv!(RESULT_CSV, _stats_row(stats, row_meta))
        catch err
            _append_exception_row!(row_meta, err)
        end
    end

    for (solver_config, add_current_row, numerical_preprocessing) in (
        ("gbc_partial", false, false),
        ("gbc_partial_new_options", true, true),
    )
        run_dir = joinpath(OUTPUT_ROOT, "runs", variant, instance.family, instance.name, solver_config)
        row_meta = _base_row_meta(
            instance,
            variant,
            solver_config,
            mps_path,
            aux_path;
            connector_add_current_solution_cut=add_current_row,
            subsolver_numerical_preprocessing=numerical_preprocessing,
            partial_decomposition=true,
        )
        if _has_existing_result(instance, variant, solver_config)
            println("  skipping existing $(solver_config)")
            continue
        end
        try
            stats = _solve_gbc_with_connector_options(
                mps_path,
                aux_path;
                partial_decomposition=true,
                out_dir=run_dir,
                connector_add_current_solution_cut=add_current_row,
                subsolver_numerical_preprocessing=numerical_preprocessing,
            )
            _append_csv!(RESULT_CSV, _stats_row(stats, row_meta))
        catch err
            _append_exception_row!(row_meta, err)
        end
    end
end

function main_benchmark_interdiction()
    mkpath(OUTPUT_ROOT)
    loaded = _load_existing_results!()
    loaded > 0 && println("Loaded $(loaded) existing rows from $(RESULT_CSV); missing configurations will be resumed.")
    instances = _discover_benchmark_interdiction_instances()
    eligible = BOBIInstance[]
    skipped = 0
    for inst in instances
        if _has_upper_level_coupling_constraints(inst.mps_path, inst.aux_path)
            skipped += 1
        else
            push!(eligible, inst)
        end
    end
    println("Selected $(length(eligible)) benchmark-interdiction instances without upper-level coupling constraints; skipped $(skipped). Output: $(RESULT_CSV)")
    println("Solver configurations: mibs, gbc_partial, gbc_partial_new_options")
    for (idx, inst) in enumerate(eligible)
        println("[$idx/$(length(eligible))] benchmark-interdiction $(inst.name)")
        _run_single_clean_benchmark(inst, "benchmark_interdiction_clean", inst.mps_path, inst.aux_path)
    end
    println("Finished. Results: $(RESULT_CSV)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_benchmark_interdiction()
end
