using JuBiC
using JuMP
using Gurobi
using Random
using Dates
using SHA

const BOBI_ROOT = joinpath(@__DIR__, "instances")
const OUTPUT_ROOT = get(ENV, "BOBILIB_OUTPUT_ROOT", joinpath("tmp_compare", "runs", "bobilib_knapsack_$(Dates.format(now(), "yyyymmdd_HHMMSS"))"))
const RESULT_CSV = joinpath(OUTPUT_ROOT, "bobilib_knapsack_results.csv")
const GENERATED_ROOT = joinpath(OUTPUT_ROOT, "generated_instances")
const UNPACKED_ROOT = joinpath(GENERATED_ROOT, "unpacked_mps")
const TIME_LIMIT = 600.0
const MAX_PER_TYPE_SIZE = 2
const MULTI_FOLLOWERS_PER_GROUP = 3
const RESULT_ROWS = Vector{Dict{String,Any}}()

struct BOBIInstance
    family::String
    item_count::Int
    name::String
    mps_path::String
    aux_path::String
end

function _csv_value(v)
    v === nothing && return ""
    s = string(v)
    return occursin(r"[,\n\"]", s) ? "\"" * replace(s, "\"" => "\"\"") * "\"" : s
end

function _append_csv!(path::AbstractString, row::Dict{String,Any})
    mkpath(dirname(path))
    push!(RESULT_ROWS, copy(row))
    keys_sorted = sort!(collect(reduce(union, (keys(r) for r in RESULT_ROWS); init=Set{String}())))
    tmp_path = path * ".tmp"
    open(tmp_path, "w") do io
        println(io, join(keys_sorted, ","))
        for result_row in RESULT_ROWS
            println(io, join((_csv_value(get(result_row, k, "")) for k in keys_sorted), ","))
        end
    end
    mv(tmp_path, path; force=true)
end

function _stats_row(stats, metadata::Dict{String,Any})
    row = copy(metadata)
    for (k, v) in stats.data
        row[string(k)] = v
    end
    return row
end

function _discover_instances(root::AbstractString=BOBI_ROOT)
    instances = BOBIInstance[]
    for (dir, _, files) in walkdir(root)
        family = _infer_family(dir)
        family === nothing && continue
        for file in files
            endswith(lowercase(file), ".aux") || continue
            aux_path = joinpath(dir, file)
            mps_name = _mps_name_from_aux(aux_path)
            mps_name === nothing && continue
            mps_path = joinpath(dir, mps_name)
            if !isfile(mps_path) && isfile(mps_path * ".gz")
                mps_path *= ".gz"
            end
            isfile(mps_path) || continue
            base = splitext(file)[1]
            push!(instances, BOBIInstance(family, _infer_item_count(family, base), base, mps_path, aux_path))
        end
    end
    sort!(instances, by=x -> (x.family, x.item_count, x.name))
    return instances
end

function _infer_family(dir::AbstractString)
    parts = Set(lowercase.(splitpath(normpath(dir))))
    "cclw" in parts && return "cclw"
    "inter-kp" in parts && return "inter-kp"
    "kp" in parts && return "kp"
    return nothing
end

function _mps_name_from_aux(aux_path::AbstractString)
    lines = readlines(aux_path)
    idx = findfirst(line -> lowercase(strip(line)) == "@mps", lines)
    if idx !== nothing && idx < length(lines)
        return strip(lines[idx + 1])
    end
    idx = findfirst(line -> endswith(lowercase(strip(line)), ".mps"), lines)
    return idx === nothing ? nothing : strip(lines[idx])
end

function _plain_mps_path(path::String)
    endswith(lowercase(path), ".gz") || return path
    rel = relpath(path, BOBI_ROOT)
    out = joinpath(UNPACKED_ROOT, rel[1:end-3])
    if !isfile(out) || filesize(out) == 0 || mtime(out) < mtime(path)
        mkpath(dirname(out))
        input_path = replace(abspath(path), "'" => "''")
        output_path = replace(abspath(out), "'" => "''")
        script = """
        \$inputPath = '$input_path'
        \$outputPath = '$output_path'
        \$inputStream = [System.IO.File]::OpenRead(\$inputPath)
        try {
            \$gzipStream = [System.IO.Compression.GzipStream]::new(\$inputStream, [System.IO.Compression.CompressionMode]::Decompress)
            try {
                \$outputStream = [System.IO.File]::Create(\$outputPath)
                try {
                    \$gzipStream.CopyTo(\$outputStream)
                } finally {
                    \$outputStream.Dispose()
                }
            } finally {
                \$gzipStream.Dispose()
            }
        } finally {
            \$inputStream.Dispose()
        }
        """
        run(`powershell -NoProfile -Command $script`)
    end
    return out
end

function _infer_item_count(family::String, name::String)
    if family == "cclw"
        m = match(r"interdiction(\d+)-", name)
        return isnothing(m) ? -1 : parse(Int, m.captures[1])
    elseif family == "inter-kp"
        m = match(r"K(\d{2,3})(\d{2})W", name)
        return isnothing(m) ? -1 : parse(Int, m.captures[1])
    elseif family == "kp"
        m = match(r"interKP-(\d+)-", name)
        return isnothing(m) ? -1 : parse(Int, m.captures[1])
    end
    return -1
end

function _selected_instances(instances)
    selected = BOBIInstance[]
    for family in sort(unique(x.family for x in instances))
        fam = filter(x -> x.family == family, instances)
        for item_count in sort(unique(x.item_count for x in fam))
            item_count < 0 && continue
            append!(selected, first(filter(x -> x.item_count == item_count, fam), min(MAX_PER_TYPE_SIZE, count(x -> x.item_count == item_count, fam))))
        end
    end
    return selected
end

function _randomize_first_level_objective(mps_path::String, aux_path::String, out_dir::String; seed::Int=42)
    mps_path = _plain_mps_path(mps_path)
    mps = JuBiC._read_mps(mps_path)
    aux = JuBiC._read_aux(aux_path)
    rng = MersenneTwister(seed)
    for obj_row in mps.rows_natural
        entries = get(mps.columns, obj_row, Tuple{String,Number}[])
        isempty(entries) && continue
        max_abs = maximum(abs(float(c)) for (_, c) in entries; init=1.0)
        mps.columns[obj_row] = [(name, round(rand(rng) * max_abs; digits=6)) for (name, _) in entries]
    end
    mkpath(out_dir)
    name = splitext(basename(mps_path))[1] * "_random_leader"
    new_mps = joinpath(out_dir, name * ".mps")
    new_aux = joinpath(out_dir, name * ".aux")
    _write_mps_data(new_mps, mps; name=name)
    cp(aux_path, new_aux; force=true)
    return new_mps, new_aux
end

function _write_mps_data(path::String, mps; name::String=mps.name)
    open(path, "w") do io
        println(io, "NAME          $(name)")
        println(io, "ROWS")
        for r in mps.rows_natural; println(io, " N  $(r)"); end
        for r in mps.rows_less_than; println(io, " L  $(r)"); end
        for r in mps.rows_greater_than; println(io, " G  $(r)"); end
        for r in mps.rows_equal; println(io, " E  $(r)"); end
        println(io, "COLUMNS")
        rows = vcat(mps.rows_natural, mps.rows_less_than, mps.rows_greater_than, mps.rows_equal)
        vars = sort!(collect(keys(mps.bounds)))
        for var in vars
            entries = [(row, coeff) for row in rows for (name, coeff) in get(mps.columns, row, Tuple{String,Number}[]) if name == var]
            for (row, coeff) in entries
                println(io, "    " * rpad(var, 20) * rpad(row, 20) * " " * string(coeff))
            end
        end
        println(io, "RHS")
        for row in rows
            rhs = get(mps.rhs, row, 0.0)
            rhs == 0.0 && continue
            println(io, "    " * rpad("rhs", 20) * rpad(row, 20) * " " * string(rhs))
        end
        println(io, "BOUNDS")
        for var in vars
            b = mps.bounds[var]
            if b.is_binary
                println(io, " BV bnd                 $(var)")
            elseif b.is_integer
                if b.lower_bound !== nothing
                    println(io, " LI bnd                 $(var) $(b.lower_bound)")
                end
                if b.upper_bound !== nothing
                    println(io, " UI bnd                 $(var) $(b.upper_bound)")
                end
                if b.lower_bound === nothing && b.upper_bound === nothing
                    println(io, " LI bnd                 $(var) 0.0")
                end
            else
                if b.lower_bound !== nothing
                    println(io, " LO bnd                 $(var) $(b.lower_bound)")
                end
                if b.upper_bound !== nothing
                    println(io, " UP bnd                 $(var) $(b.upper_bound)")
                end
            end
        end
        println(io, "ENDATA")
    end
end

function _solve_gbc(mps_path, aux_path; partial_decomposition::Bool, out_dir::String)
    mps_path = _plain_mps_path(mps_path)
    stats = JuBiC.RunStats()
    inst = get_GBC_instance(mps_path, aux_path, Gurobi.Optimizer; partial_decomposition=partial_decomposition, preprocessing=true, stats=stats)
    params = GBCparam(GurobiSolver(), false, out_dir, "lp", stats, TIME_LIMIT, 42, 8, 8, true, PARETO_OPTIMALITY_ONLY, true, false, true, 1e9, 0, false)
    new_stat!(params.stats, "enable_output_logs", false)
    return solve_instance!(inst, params)
end

function _solve_mibs(mps_path, aux_path; out_dir::String)
    mps_path = _plain_mps_path(mps_path)
    stats = JuBiC.RunStats()
    inst = get_MibS_instance(mps_path, aux_path; stats=stats)
    params = MibSparam(false, out_dir, TIME_LIMIT, stats)
    new_stat!(params.stats, "enable_output_logs", false)
    return solve_instance!(inst, params)
end

function _run_single(instance::BOBIInstance, variant::String, mps_path::String, aux_path::String)
    for (solver_name, partial) in (("gbc_partial", true), ("gbc_no_partial", false))
        run_dir = joinpath(OUTPUT_ROOT, "runs", variant, instance.family, instance.name, solver_name)
        row_meta = Dict{String,Any}(
            "experiment" => variant,
            "solver_config" => solver_name,
            "family" => instance.family,
            "item_count" => instance.item_count,
            "instance" => instance.name,
            "source_mps" => mps_path,
            "source_aux" => aux_path,
            "runtime_limit" => TIME_LIMIT,
        )
        try
            stats = _solve_gbc(mps_path, aux_path; partial_decomposition=partial, out_dir=run_dir)
            _append_csv!(RESULT_CSV, _stats_row(stats, row_meta))
        catch err
            row_meta["RunStatus"] = "Exception"
            row_meta["Error"] = sprint(showerror, err)
            _append_csv!(RESULT_CSV, row_meta)
        end
    end
    run_dir = joinpath(OUTPUT_ROOT, "runs", variant, instance.family, instance.name, "mibs")
    row_meta = Dict{String,Any}(
        "experiment" => variant,
        "solver_config" => "mibs",
        "family" => instance.family,
        "item_count" => instance.item_count,
        "instance" => instance.name,
        "source_mps" => mps_path,
        "source_aux" => aux_path,
        "runtime_limit" => TIME_LIMIT,
    )
    try
        stats = _solve_mibs(mps_path, aux_path; out_dir=run_dir)
        _append_csv!(RESULT_CSV, _stats_row(stats, row_meta))
    catch err
        row_meta["RunStatus"] = "Exception"
        row_meta["Error"] = sprint(showerror, err)
        _append_csv!(RESULT_CSV, row_meta)
    end
end

function _run_multi(group::Vector{BOBIInstance})
    isempty(group) && return
    base = first(group)
    group_digest = bytes2hex(sha1(join((x.name for x in group), "_")))[1:8]
    group_name = "$(base.family)_n$(base.item_count)_$(group_digest)"
    generated_dir = joinpath(GENERATED_ROOT, "multi", group_name)
    mkpath(generated_dir)
    open(joinpath(generated_dir, "manifest.txt"), "w") do io
        println(io, "family=$(base.family)")
        println(io, "item_count=$(base.item_count)")
        for inst in group
            println(io, "follower=$(inst.name),mps=$(inst.mps_path),aux=$(inst.aux_path)")
        end
    end

    for (solver_name, partial) in (("gbc_partial", true), ("gbc_no_partial", false))
        run_dir = joinpath(OUTPUT_ROOT, "runs", "multi_random_leader", base.family, group_name, solver_name)
        row_meta = Dict{String,Any}(
            "experiment" => "multi_random_leader",
            "solver_config" => solver_name,
            "family" => base.family,
            "item_count" => base.item_count,
            "instance" => group_name,
            "nfollowers" => length(group),
            "runtime_limit" => TIME_LIMIT,
        )
        try
            stats = JuBiC.RunStats()
            instance = _build_multi_gbc_instance(group; partial_decomposition=partial, stats=stats)
            if !partial
                export_dir = joinpath(generated_dir, "jubic_gbc")
                try
                    output_GBC_solver_instance(instance, export_dir)
                    row_meta["generated_jubic_instance"] = export_dir
                catch export_err
                    row_meta["generated_jubic_export_error"] = sprint(showerror, export_err)
                end
            end
            params = GBCparam(GurobiSolver(), false, run_dir, "lp", stats, TIME_LIMIT, 42, 8, 8, true, PARETO_OPTIMALITY_ONLY, true, false, true, 1e9, 0, false)
            new_stat!(params.stats, "enable_output_logs", false)
            solved = solve_instance!(instance, params)
            _append_csv!(RESULT_CSV, _stats_row(solved, row_meta))
        catch err
            row_meta["RunStatus"] = "Exception"
            row_meta["Error"] = sprint(showerror, err)
            _append_csv!(RESULT_CSV, row_meta)
        end
    end

    run_dir = joinpath(OUTPUT_ROOT, "runs", "multi_random_leader", base.family, group_name, "mibs_merged")
    row_meta = Dict{String,Any}(
        "experiment" => "multi_random_leader",
        "solver_config" => "mibs_merged",
        "family" => base.family,
        "item_count" => base.item_count,
        "instance" => group_name,
        "nfollowers" => length(group),
        "runtime_limit" => TIME_LIMIT,
    )
    try
        stats = JuBiC.RunStats()
        gbc_instance = _build_multi_gbc_instance(group; partial_decomposition=false, stats=stats)
        mibs_instance = transform_GBC_to_MibS(gbc_instance, GurobiSolver())
        export_dir = joinpath(generated_dir, "mibs_merged")
        try
            output_MibS_instance(mibs_instance, group_name * "_merged", export_dir)
            row_meta["generated_mibs_instance"] = export_dir
        catch export_err
            row_meta["generated_mibs_export_error"] = sprint(showerror, export_err)
        end
        params = MibSparam(false, run_dir, TIME_LIMIT, stats)
        new_stat!(params.stats, "enable_output_logs", false)
        solved = solve_instance!(mibs_instance, params)
        _append_csv!(RESULT_CSV, _stats_row(solved, row_meta))
    catch err
        row_meta["RunStatus"] = "Exception"
        row_meta["Error"] = sprint(showerror, err)
        _append_csv!(RESULT_CSV, row_meta)
    end
end

function _build_multi_gbc_instance(group::Vector{BOBIInstance}; partial_decomposition::Bool, stats::JuBiC.RunStats)
    base_mps = JuBiC._read_mps(_plain_mps_path(first(group).mps_path))
    base_aux = JuBiC._read_aux(first(group).aux_path)
    _randomize_mps_objective!(base_mps; seed=first(group).item_count + length(group))
    solver = Gurobi.Optimizer
    model_upper = Model(solver)
    upper_vars = Dict{String,Any}()
    upper_vars_partial = Dict{String,Any}()
    for (name, bound) in base_mps.bounds
        name in base_aux.variables && continue
        lb = isnothing(bound.lower_bound) ? -Inf : bound.lower_bound
        ub = isnothing(bound.upper_bound) ? Inf : bound.upper_bound
        var_ref = bound.is_binary ? @variable(model_upper, base_name=name, binary=true) :
            (bound.is_integer ? @variable(model_upper, base_name=name, integer=true, lower_bound=lb, upper_bound=ub) :
             @variable(model_upper, base_name=name, lower_bound=lb, upper_bound=ub))
        upper_vars[name] = var_ref
        upper_vars_partial[name] = var_ref
    end
    @assert length(base_mps.rows_natural) <= 1
    if length(base_mps.rows_natural) == 1
        obj_row = base_mps.rows_natural[1]
        @objective(model_upper, Min, sum(upper_vars[name] * c for (name, c) in get(base_mps.columns, obj_row, []) if haskey(upper_vars, name); init=0.0))
    else
        @objective(model_upper, Min, 0.0)
    end

    all_A = String[]
    xdict = Dict{String,VariableRef}()
    subs = SubSolver[]
    sub_names = String[]

    for (idx, inst) in enumerate(group)
        mps = JuBiC._read_mps(_plain_mps_path(inst.mps_path))
        aux = JuBiC._read_aux(inst.aux_path)
        _randomize_mps_objective!(mps; seed=idx + inst.item_count)
        regular_linking, regular_linking_data, link_vars = JuBiC._classify_variables(mps, aux, true)
        link_vars_aux = String[var_name * "_F$(idx)" * JuBiC.VARIABLE_AUX_SUFFIX for var_name in link_vars]
        link_map = Dict(var_name => var_name * "_F$(idx)" for var_name in link_vars)
        A_local = vcat(String["F$(idx)::" * a for a in regular_linking],
                       String["F$(idx)::" * a for a in link_vars],
                       String["F$(idx)::" * a for a in link_vars_aux])
        append!(all_A, A_local)

        model_lower = Model(solver)
        lower_vars = Dict{String,Any}()
        for (name, bound) in mps.bounds
            if name in aux.variables
                lb = isnothing(bound.lower_bound) ? -Inf : bound.lower_bound
                ub = isnothing(bound.upper_bound) ? Inf : bound.upper_bound
                lower_vars[name] = bound.is_binary ? @variable(model_lower, base_name="F$(idx)_$name", binary=true) :
                    (bound.is_integer ? @variable(model_lower, base_name="F$(idx)_$name", integer=true, lower_bound=lb, upper_bound=ub) :
                     @variable(model_lower, base_name="F$(idx)_$name", lower_bound=lb, upper_bound=ub))
            end
        end
        @variable(model_lower, link_lower[A_local], Bin)
        @variable(model_lower, upper_copy[A_local], Bin)
        aux_A = String["F$(idx)::" * a for a in link_vars_aux]
        @variable(model_upper, link_upper_aux[aux_A], Bin)

        y_vars = JuMP.Containers.DenseAxisArray{VariableRef}(undef, A_local)
        local_link_varsC = JuMP.Containers.DenseAxisArray{VariableRef}(undef, A_local)
        for a in A_local
            raw = replace(a, "F$(idx)::" => "")
            if raw in regular_linking
                y_vars[a] = lower_vars[regular_linking_data[raw][1]]
                xdict[a] = upper_vars[raw]
            else
                y_vars[a] = link_lower[a]
                if endswith(raw, JuBiC.VARIABLE_AUX_SUFFIX)
                    original = replace(raw, "_F$(idx)" * JuBiC.VARIABLE_AUX_SUFFIX => "")
                    xdict[a] = link_upper_aux[a]
                    @constraint(model_upper, link_upper_aux[a] + upper_vars[original] == 1)
                else
                    xdict[a] = upper_vars[raw]
                end
            end
            local_link_varsC[a] = upper_copy[a]
            @constraint(model_lower, link_lower[a] <= upper_copy[a])
        end

        lower_obj = @expression(model_lower, sum(aux.objective[name] * var for (name, var) in lower_vars; init=0.0))
        @objective(model_lower, Min, lower_obj)
        master_sub_obj = length(mps.rows_natural) == 1 ?
            sum(lower_vars[name] * c for (name, c) in get(mps.columns, mps.rows_natural[1], []) if haskey(lower_vars, name); init=AffExpr(0)) :
            AffExpr(0)
        _add_multi_constraints!(model_upper, model_lower, mps, aux, upper_vars, lower_vars, link_lower, upper_copy, A_local, partial_decomposition, upper_vars_partial)
        subname = "F$(idx)_$(inst.name)"
        push!(sub_names, subname)
        push!(subs, SubSolverJuMP(subname, model_lower, A_local, local_link_varsC, y_vars, master_sub_obj, lower_obj, timelimit -> (false, 0)))
        set_silent(model_lower)
    end
    set_silent(model_upper)
    master = Master(model_upper, all_A, xdict, sub_names)
    new_stat!(stats, "instance", "multi_$(first(group).family)_$(first(group).item_count)")
    new_stat!(stats, "instanceGen_multifollower", true)
    new_stat!(stats, "instanceGen_nfollowers", length(group))
    new_stat!(stats, "instanceGen_use_partial_decomposition", partial_decomposition)
    return Instance(master, subs)
end

function _randomize_mps_objective!(mps; seed::Int)
    rng = MersenneTwister(seed)
    for obj_row in mps.rows_natural
        entries = get(mps.columns, obj_row, Tuple{String,Number}[])
        max_abs = maximum(abs(float(c)) for (_, c) in entries; init=1.0)
        mps.columns[obj_row] = [(name, round(rand(rng) * max_abs; digits=6)) for (name, _) in entries]
    end
end

function _add_multi_constraints!(model_upper, model_lower, mps, aux, upper_vars, lower_vars, link_lower, upper_copy, A_local, partial_decomposition, upper_vars_partial)
    function add_rows(rows, sense)
        for row in rows
            sub_model_is_upper = !(row in aux.constraints)
            if sub_model_is_upper
                lhs = sum(upper_vars[name] * c for (name, c) in get(mps.columns, row, []) if haskey(upper_vars, name); init=0.0)
                sense == 'L' ? @constraint(model_upper, lhs <= mps.rhs[row]) :
                sense == 'G' ? @constraint(model_upper, lhs >= mps.rhs[row]) :
                               @constraint(model_upper, lhs == mps.rhs[row])
            else
                lhs = sum(lower_vars[name] * c for (name, c) in get(mps.columns, row, []) if haskey(lower_vars, name); init=0.0)
                for (name, c) in get(mps.columns, row, [])
                    haskey(upper_vars, name) || continue
                    for a in A_local
                        endswith(a, "::" * name) && add_to_expression!(lhs, c, link_lower[a])
                    end
                end
                sense == 'L' ? @constraint(model_lower, lhs <= mps.rhs[row]) :
                sense == 'G' ? @constraint(model_lower, lhs >= mps.rhs[row]) :
                               @constraint(model_lower, lhs == mps.rhs[row])
                if partial_decomposition
                    lhs_partial = sum(upper_vars_partial[name] * c for (name, c) in get(mps.columns, row, []) if haskey(upper_vars_partial, name); init=0.0)
                    sense == 'L' ? @constraint(model_upper, lhs_partial <= mps.rhs[row]) :
                    sense == 'G' ? @constraint(model_upper, lhs_partial >= mps.rhs[row]) :
                                   @constraint(model_upper, lhs_partial == mps.rhs[row])
                end
            end
        end
    end
    add_rows(mps.rows_less_than, 'L')
    add_rows(mps.rows_greater_than, 'G')
    add_rows(mps.rows_equal, 'E')
end

function main()
    mkpath(OUTPUT_ROOT)
    instances = _selected_instances(_discover_instances())
    println("Selected $(length(instances)) BOBILib instances. Output: $(RESULT_CSV)")
    for inst in instances
        println("Original $(inst.family) $(inst.name)")
        _run_single(inst, "original_interdiction", inst.mps_path, inst.aux_path)
        randomized_dir = joinpath(GENERATED_ROOT, "random_leader", inst.family, inst.name)
        random_mps, random_aux = _randomize_first_level_objective(inst.mps_path, inst.aux_path, randomized_dir; seed=inst.item_count)
        _run_single(inst, "random_leader", random_mps, random_aux)
    end
    for family in sort(unique(x.family for x in instances))
        fam = filter(x -> x.family == family, instances)
        for item_count in sort(unique(x.item_count for x in fam))
            group = first(filter(x -> x.item_count == item_count, fam), min(MULTI_FOLLOWERS_PER_GROUP, count(x -> x.item_count == item_count, fam)))
            length(group) >= 2 || continue
            println("Multi $(family) n=$(item_count) followers=$(length(group))")
            _run_multi(group)
        end
    end
    println("Finished. Results: $(RESULT_CSV)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
