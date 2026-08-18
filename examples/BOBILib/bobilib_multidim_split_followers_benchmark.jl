using Dates
using SHA

const MULTIDIM_ROOT = joinpath(
    @__DIR__,
    "instances",
    "multidimensional-knapsack-interd",
    "multidimensional-knapsack-interdiction",
)

if !haskey(ENV, "BOBILIB_OUTPUT_ROOT")
    ENV["BOBILIB_OUTPUT_ROOT"] = joinpath(
        "tmp_compare",
        "runs",
        "bobilib_multidim_split_followers_$(Dates.format(now(), "yyyymmdd_HHMMSS"))",
    )
end

include(joinpath(@__DIR__, "bobilib_knapsack_benchmark.jl"))

const MULTIDIM_TIME_LIMIT = parse(Float64, get(ENV, "BOBILIB_TIME_LIMIT", string(TIME_LIMIT)))

function _discover_multidim_instances(root::AbstractString=MULTIDIM_ROOT)
    instances = BOBIInstance[]
    for aux_path in sort!(collect(Base.Filesystem.walkdir(root)) |> dirs -> reduce(vcat, (joinpath(dir, f) for (dir, _, files) in dirs for f in files if endswith(lowercase(f), ".aux") && !startswith(f, "._")); init=String[]))
        mps_name = _mps_name_from_aux(aux_path)
        mps_name === nothing && continue
        startswith(mps_name, "._") && continue
        mps_path = joinpath(dirname(aux_path), mps_name)
        if !isfile(mps_path) && isfile(mps_path * ".gz")
            mps_path *= ".gz"
        end
        isfile(mps_path) || error("Missing MPS file for AUX $(aux_path): expected $(mps_name) or $(mps_name).gz")
        family = splitpath(dirname(aux_path))[end]
        name = splitext(basename(aux_path))[1]
        push!(instances, BOBIInstance(family, _infer_item_count(family, name), name, mps_path, aux_path))
    end
    return instances
end

function _filtered_multidim_instances()
    instances = _discover_multidim_instances()
    pattern = get(ENV, "BOBILIB_MULTIDIM_PATTERN", "")
    if !isempty(pattern)
        instances = filter(inst -> occursin(pattern, inst.name), instances)
    end
    family = get(ENV, "BOBILIB_MULTIDIM_FAMILY", "")
    if !isempty(family)
        instances = filter(inst -> inst.family == family, instances)
    end
    max_instances = tryparse(Int, get(ENV, "BOBILIB_MULTIDIM_MAX", "0"))
    if max_instances !== nothing && max_instances > 0
        instances = first(instances, min(max_instances, length(instances)))
    end
    return instances
end

_is_interdiction_row(row::AbstractString) = startswith(row, "I")

function _multidim_knapsack_rows(aux_data)
    rows = [row for row in aux_data.constraints if !_is_interdiction_row(row)]
    isempty(rows) && error("No follower knapsack rows found in AUX data. Expected non-I* follower rows such as KF0.")
    return rows
end

function _multidim_interdiction_rows(aux_data)
    return [row for row in aux_data.constraints if _is_interdiction_row(row)]
end

function _create_bobi_var!(model::JuMP.Model, name::String, bound)
    lb = isnothing(bound.lower_bound) ? -Inf : bound.lower_bound
    ub = isnothing(bound.upper_bound) ? Inf : bound.upper_bound
    if bound.is_binary
        return @variable(model, base_name = name, binary = true)
    elseif bound.is_integer
        return @variable(model, base_name = name, integer = true, lower_bound = lb, upper_bound = ub)
    else
        return @variable(model, base_name = name, lower_bound = lb, upper_bound = ub)
    end
end

function _add_rows_from_mps!(model, rows, sense, mps_data, var_lookup::Dict)
    for row in rows
        lhs = sum(var_lookup[name] * coef for (name, coef) in get(mps_data.columns, row, Tuple{String,Float64}[]) if haskey(var_lookup, name); init=0.0)
        rhs = mps_data.rhs[row]
        if sense == 'L'
            @constraint(model, lhs <= rhs)
        elseif sense == 'G'
            @constraint(model, lhs >= rhs)
        elseif sense == 'E'
            @constraint(model, lhs == rhs)
        else
            error("Unsupported row sense $(sense)")
        end
    end
end

function _add_one_multidim_follower_rows!(model_lower, mps_data, row_set::Set{String}, var_lookup::Dict)
    _add_rows_from_mps!(model_lower, [row for row in mps_data.rows_less_than if row in row_set], 'L', mps_data, var_lookup)
    _add_rows_from_mps!(model_lower, [row for row in mps_data.rows_greater_than if row in row_set], 'G', mps_data, var_lookup)
    _add_rows_from_mps!(model_lower, [row for row in mps_data.rows_equal if row in row_set], 'E', mps_data, var_lookup)
end

function _build_multidim_split_followers_instance(mps_path::String, aux_path::String; partial_decomposition::Bool=true, stats::JuBiC.RunStats=JuBiC.RunStats())
    mps_data = JuBiC._read_mps(_plain_mps_path(mps_path))
    aux_data = JuBiC._read_aux(aux_path)
    JuBiC._assert_no_upper_level_coupling_constraints(mps_data, aux_data; importer="multidimensional split-follower GBC builder")
    obj_bias = JuBiC._preprocess_model(mps_data, aux_data)

    knapsack_rows = _multidim_knapsack_rows(aux_data)
    interdiction_rows = _multidim_interdiction_rows(aux_data)

    model_upper = Model(Gurobi.Optimizer)
    upper_vars = Dict{String,Any}()
    for (name, bound) in mps_data.bounds
        name in aux_data.variables && continue
        upper_vars[name] = _create_bobi_var!(model_upper, name, bound)
    end

    regular_linking, regular_linking_data, link_vars = JuBiC._classify_variables(mps_data, aux_data, true)
    link_vars_aux = String[var_name * JuBiC.VARIABLE_AUX_SUFFIX for var_name in link_vars]
    all_link_vars = [link_vars; link_vars_aux]
    A = [regular_linking; link_vars; link_vars_aux]

    @variable(model_upper, link_upper_aux[link_vars_aux], Bin)
    @constraint(model_upper, [a in link_vars], link_upper_aux[a * JuBiC.VARIABLE_AUX_SUFFIX] + upper_vars[a] == 1)

    xdict = merge(
        Dict{String,VariableRef}(a => upper_vars[a] for a in regular_linking),
        Dict{String,VariableRef}(a => upper_vars[a] for a in link_vars),
        Dict{String,VariableRef}(a => link_upper_aux[a] for a in link_vars_aux),
    )

    if length(mps_data.rows_natural) == 1 && haskey(mps_data.columns, mps_data.rows_natural[1])
        obj_row = mps_data.rows_natural[1]
        @objective(model_upper, Min, sum(upper_vars[name] * coef for (name, coef) in get(mps_data.columns, obj_row, []) if haskey(upper_vars, name); init=obj_bias))
    else
        @objective(model_upper, Min, obj_bias)
    end

    upper_row_lookup = Dict{String,Any}(upper_vars)
    _add_rows_from_mps!(model_upper, [row for row in mps_data.rows_less_than if !(row in aux_data.constraints)], 'L', mps_data, upper_row_lookup)
    _add_rows_from_mps!(model_upper, [row for row in mps_data.rows_greater_than if !(row in aux_data.constraints)], 'G', mps_data, upper_row_lookup)
    _add_rows_from_mps!(model_upper, [row for row in mps_data.rows_equal if !(row in aux_data.constraints)], 'E', mps_data, upper_row_lookup)

    sub_names = String[]
    subs = []
    for (idx, knapsack_row) in enumerate(knapsack_rows)
        sub_name = "KF$(idx)_$(knapsack_row)"
        model_lower = Model(Gurobi.Optimizer)
        lower_vars = Dict{String,Any}()
        for (name, bound) in mps_data.bounds
            name in aux_data.variables || continue
            lower_vars[name] = _create_bobi_var!(model_lower, "$(sub_name)_$(name)", bound)
        end

        @variable(model_lower, link_lower[all_link_vars], Bin)
        @variable(model_lower, upper_copy[A], Bin)
        @constraint(model_lower, [a in link_vars], link_lower[a * JuBiC.VARIABLE_AUX_SUFFIX] + link_lower[a] == 1)
        @constraint(model_lower, [a in all_link_vars], link_lower[a] <= upper_copy[a])

        y_vars = JuMP.Containers.DenseAxisArray{VariableRef}(undef, A)
        for a in A
            y_vars[a] = if a in regular_linking
                lower_vars[regular_linking_data[a][1]]
            else
                link_lower[a]
            end
        end

        lower_obj = @expression(model_lower, sum(aux_data.objective[name] * var for (name, var) in lower_vars; init=0.0))
        @objective(model_lower, Min, lower_obj)

        master_sub_obj = if length(mps_data.rows_natural) == 1 && haskey(mps_data.columns, mps_data.rows_natural[1])
            sum(lower_vars[name] * coef for (name, coef) in get(mps_data.columns, mps_data.rows_natural[1], []) if haskey(lower_vars, name); init=AffExpr(0))
        else
            AffExpr(0)
        end

        row_set = Set([knapsack_row; interdiction_rows])
        var_lookup = merge(
            Dict{String,Any}(lower_vars),
            Dict{String,Any}(name => link_lower[name] for name in link_vars),
            Dict{String,Any}(name => upper_copy[name] for name in regular_linking),
        )
        _add_one_multidim_follower_rows!(model_lower, mps_data, row_set, var_lookup)

        push!(sub_names, sub_name)
        push!(subs, SubSolverJuMP(sub_name, model_lower, A, y_vars, master_sub_obj, lower_obj, timelimit -> (false, 0)))
        set_silent(model_lower)
    end

    new_stat!(stats, "instanceGen_multidim_split_followers", true)
    new_stat!(stats, "instanceGen_nfollowers", length(knapsack_rows))
    new_stat!(stats, "instanceGen_knapsack_rows", join(knapsack_rows, ";"))
    new_stat!(stats, "instanceGen_interdiction_rows", length(interdiction_rows))
    new_stat!(stats, "instanceGen_use_partial_decomposition", partial_decomposition)
    new_stat!(stats, "instanceGen_obj_bias", obj_bias)

    master = Master(model_upper, A, xdict, sub_names)
    return Instance(master, subs)
end

function _write_split_manifest(path::String, instance::BOBIInstance, knapsack_rows, generated_dir::String)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, "source_mps=$(instance.mps_path)")
        println(io, "source_aux=$(instance.aux_path)")
        println(io, "family=$(instance.family)")
        println(io, "instance=$(instance.name)")
        println(io, "split_rule=one follower per non-I* AUX follower constraint")
        println(io, "nfollowers=$(length(knapsack_rows))")
        println(io, "followers=$(join(knapsack_rows, ","))")
        println(io, "generated_dir=$(generated_dir)")
    end
end

function _solve_multidim_split_gbc(instance::BOBIInstance, generated_dir::String, run_dir::String)
    stats = JuBiC.RunStats()
    inst = _build_multidim_split_followers_instance(instance.mps_path, instance.aux_path; partial_decomposition=true, stats=stats)
    export_dir = joinpath(generated_dir, "gbc_multifollower_export")
    try
        output_GBC_instance(inst, instance.name * "_split_followers", export_dir, GurobiSolver(); anonymous_names=false)
    catch err
        @warn "Could not export split-follower GBC instance for $(instance.name)" exception=(err, catch_backtrace())
    end
    params = GBCparam(GurobiSolver(), false, run_dir, "lp", stats, MULTIDIM_TIME_LIMIT, 42, 8, 8, true, PARETO_OPTIMALITY_ONLY, true, false, true, 1e9, 0, false)
    new_stat!(params.stats, "enable_output_logs", false)
    return solve_instance!(inst, params)
end

function _solve_multidim_split_mibs(instance::BOBIInstance, generated_dir::String, run_dir::String)
    stats = JuBiC.RunStats()
    gbc_instance = _build_multidim_split_followers_instance(instance.mps_path, instance.aux_path; partial_decomposition=false, stats=stats)
    mibs_instance = transform_GBC_to_MibS(gbc_instance, GurobiSolver())
    export_dir = joinpath(generated_dir, "mibs_aggregated_export")
    try
        output_MibS_instance(mibs_instance, instance.name * "_split_followers_aggregated", export_dir; anonymous_names=false)
    catch err
        @warn "Could not export aggregated MiBS instance for $(instance.name)" exception=(err, catch_backtrace())
    end
    params = MibSparam(false, run_dir, MULTIDIM_TIME_LIMIT, stats)
    new_stat!(params.stats, "enable_output_logs", false)
    return solve_instance!(mibs_instance, params)
end

function _run_multidim_split_instance(instance::BOBIInstance)
    mps_data = JuBiC._read_mps(_plain_mps_path(instance.mps_path))
    aux_data = JuBiC._read_aux(instance.aux_path)
    if JuBiC._has_upper_level_coupling_constraints(mps_data, aux_data)
        println("  skipping $(instance.name): upper-level constraints contain follower variables")
        return
    end
    knapsack_rows = _multidim_knapsack_rows(aux_data)
    generated_dir = joinpath(GENERATED_ROOT, "multidim_split_followers", instance.family, instance.name)
    _write_split_manifest(joinpath(generated_dir, "manifest.txt"), instance, knapsack_rows, generated_dir)

    for (solver_name, solve_fn) in (
        ("gbc_partial", _solve_multidim_split_gbc),
        ("mibs_aggregated", _solve_multidim_split_mibs),
    )
        run_dir = joinpath(OUTPUT_ROOT, "runs", "multidim_split_followers", instance.family, instance.name, solver_name)
        row_meta = Dict{String,Any}(
            "experiment" => "multidim_split_followers",
            "solver_config" => solver_name,
            "family" => instance.family,
            "item_count" => instance.item_count,
            "instance" => instance.name,
            "source_mps" => instance.mps_path,
            "source_aux" => instance.aux_path,
            "runtime_limit" => MULTIDIM_TIME_LIMIT,
            "nfollowers" => length(knapsack_rows),
        )
        try
            stats = solve_fn(instance, generated_dir, run_dir)
            _append_csv!(RESULT_CSV, _stats_row(stats, row_meta))
        catch err
            row_meta["RunStatus"] = "Exception"
            row_meta["Error"] = sprint(showerror, err)
            _append_csv!(RESULT_CSV, row_meta)
        end
    end
end

function main_multidim_split_followers()
    mkpath(OUTPUT_ROOT)
    instances = _filtered_multidim_instances()
    println("Selected $(length(instances)) multidimensional split-follower instances. Output: $(RESULT_CSV)")
    println("Solver configurations: gbc_partial, mibs_aggregated")
    if Threads.nthreads() == 1
        @warn "Julia is running with one thread. Start this script with multiple Julia threads, e.g. julia --threads=8 --project=. examples/BOBILib/bobilib_multidim_split_followers_benchmark.jl"
    end
    for (idx, inst) in enumerate(instances)
        println("[$idx/$(length(instances))] multidim-split $(inst.family) $(inst.name)")
        _run_multidim_split_instance(inst)
    end
    println("Finished. Results: $(RESULT_CSV)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    main_multidim_split_followers()
end
