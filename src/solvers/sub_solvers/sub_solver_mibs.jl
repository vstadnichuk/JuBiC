using JuMP, BilevelJuMP, MibS_jll
import MathOptInterface as MOI

"""
An implementation of a `SubSolver` that uses a `BilevelModel` and MibS.
It mirrors the `SubSolverJuMP` logic, but preserves bilevel feasibility in the
BlC/BlCLag separation routines by keeping the follower objective in the lower
level and optimizing the separation objective in the upper level.
"""
struct SubSolverMiBS{T} <: SubSolver
    name::String
    bi_model::BilevelModel
    A::Vector{T}
    link_varsC::Any
    y_vars::Any
    r_objterm::GenericAffExpr
    c_objterm::GenericAffExpr
    lower_name_map::Dict{String,String}
    upper_name_map::Dict{String,String}
end

function _create_link_varsC(bi_model::BilevelModel, sub_name, A, y_vars)
    link_varsC = @variable(Upper(bi_model), [a in A], Bin, base_name = "x_copy_$(sub_name)")
    @constraint(Lower(bi_model), [a in A], y_vars[a] <= link_varsC[a])
    return link_varsC
end

function SubSolverMiBS(
    name,
    bi_model,
    A::Vector{T},
    y_vars,
    r_objterm,
    c_objterm,
) where T
    link_varsC = _create_link_varsC(bi_model, name, A, y_vars)
    return SubSolverMiBS{T}(
        name,
        bi_model,
        A,
        link_varsC,
        y_vars,
        r_objterm,
        c_objterm,
        Dict{String,String}(),
        Dict{String,String}(),
    )
end

function SubSolverMiBS(
    name,
    bi_model,
    A::Vector{T},
    link_varsC,
    y_vars,
    r_objterm,
    c_objterm,
) where T
    return SubSolverMiBS{T}(
        name,
        bi_model,
        A,
        link_varsC,
        y_vars,
        r_objterm,
        c_objterm,
        Dict{String,String}(),
        Dict{String,String}(),
    )
end

function SubSolverMiBS(
    name,
    bi_model,
    A::Vector{T},
    link_varsC,
    y_vars,
    link_constraints_capacities,
    r_objterm,
    c_objterm,
) where T
    @warn "The deprecated SubSolverMiBS constructor with explicit link_varsC and link_constraints_capacities was used for subsolver $(name). The passed linking copies and capacities are ignored; JuBiC now generates internal linking copies and assumes unit capacities."
    return SubSolverMiBS(name, bi_model, A, y_vars, r_objterm, c_objterm)
end

function capacity_linking(sol::SubSolverMiBS, a, params::SolverParam)
    return 1
end

_mibs_safe_name(name::AbstractString) = replace(String(name), ' ' => '_')

function _initialize_name_maps!(sol::SubSolverMiBS)
    if !isempty(sol.lower_name_map) && !isempty(sol.upper_name_map)
        return nothing
    end

    empty!(sol.lower_name_map)
    empty!(sol.upper_name_map)
    for v in all_variables(Lower(sol.bi_model))
        original_name = JuMP.name(v)
        sol.lower_name_map[_mibs_safe_name(original_name)] = original_name
    end
    for v in all_variables(Upper(sol.bi_model))
        original_name = JuMP.name(v)
        sol.upper_name_map[_mibs_safe_name(original_name)] = original_name
    end
    return nothing
end

function check(sol::SubSolverMiBS, params::SolverParam)
    if !(length(sol.link_varsC) == length(sol.y_vars))
        error(
            "In subsolver $(name(sol)), the number of linking variables does not match the number of second level variables in linking constraints! |link_vars| = $(length(sol.link_varsC)) and |y_vars| = $(length(sol.y_vars)).",
        )
    elseif !(length(sol.y_vars) == length(sol.A))
        error(
            "In subsolver $(name(sol)), the number of second level variables in linking constraints does not match the number of resources. |A| = $(length(sol.A)) and |y_vars| = $(length(sol.y_vars)).",
        )
    end

    for a in sol.A
        if !is_valid(Upper(sol.bi_model), sol.link_varsC[a])
            error("Linking variable for resource $(a) is not contained in upper level of subsolver $(name(sol))!")
        elseif !is_valid(Lower(sol.bi_model), sol.y_vars[a])
            error("Second level variable (y) for resource $(a) is not contained in lower level of subsolver $(name(sol))!")
        end
    end
end

function compute_lower_bound_master_contribution(sol::SubSolverMiBS, params::SolverParam, time_limit)
    @objective(Upper(sol.bi_model), Min, sol.r_objterm)

    try
        if should_debbug_print(params)
            write_to_file(
                sol.bi_model,
                "$(params.output_folder_path)/ubbound_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    solution = _solve_mip(sol, params.debbug_out, time_limit)
    if !solution.status
        if solution.timeout
            throw(TimeoutException("We hit the time limit while solving MiBS subsolver $(sol.name)."))
        end
        throw(MibSFailureException("MiBS failed to solve submodel $(sol.name) to optimality when computing the lower bound."))
    end

    @debug "The found lower bound solution for $(sol.name) is $(solution.objective)."
    return solution.objective
end

function name(ss::SubSolverMiBS)
    return ss.name
end

function separation!(sol::SubSolverMiBS, sval, gvals, kvals::Dict, params::SolverParam, time_limit)
    var_non_zero_tolerance = 1e-3

    k_term = sum(kvals[a] * sol.link_varsC[a] for a in sol.A)
    new_obj = sol.c_objterm * gvals + sol.r_objterm + k_term
    @objective(Upper(sol.bi_model), Min, new_obj)

    try
        if should_debbug_print(params)
            write_to_file(
                sol.bi_model,
                "$(params.output_folder_path)/sub_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end
    solution = _solve_mip(sol, params.debbug_out, time_limit)

    if !solution.status
        if solution.timeout
            throw(TimeoutException("We hit the time limit while solving MiBS subsolver $(sol.name) in GBC separation."))
        end
        throw(MibSFailureException("MiBS failed while solving subproblem $(sol.name) for separating GBC cuts."))
    end

    opt_obj = solution.objective
    @debug "Optimal solution of sub $(sol.name) is $(opt_obj)"

    is_violate = Bool(!(sval < opt_obj + var_non_zero_tolerance))
    @debug "The violated status of sub $(sol.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    r = eval_affexpr_by_name(sol.r_objterm, sol.bi_model, solution.all_upper, solution.all_lower)
    c = eval_affexpr_by_name(sol.c_objterm, sol.bi_model, solution.all_upper, solution.all_lower)
    as = indices_where_coeff_pos(solution.all_upper, sol.link_varsC)
    return SubSolution(is_violate, r, c, opt_obj, as)
end

function separation_BlC!(sub_solver::SubSolverMiBS, sval, kvals::Dict, params::SolverParam, time_limit)
    var_non_zero_tolerance = 1e-3

    k_term = sum(kvals[a] * (1 - sub_solver.link_varsC[a]) for a in sub_solver.A)
    k_penalty = sum(var_non_zero_tolerance * (1 - sub_solver.link_varsC[a]) for a in sub_solver.A)
    new_obj = sub_solver.c_objterm - k_term - k_penalty
    @objective(Upper(sub_solver.bi_model), Max, new_obj)

    try
        if should_debbug_print(params)
            write_to_file(
                sub_solver.bi_model,
                "$(params.output_folder_path)/sub_$(sub_solver.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sub_solver.name) to file. error message $err"
    end
    solution = _solve_mip(sub_solver, params.debbug_out, time_limit)

    if !solution.status
        if solution.timeout
            throw(TimeoutException("We hit the time limit while solving MiBS subsolver $(sub_solver.name) for separating BlC cuts."))
        end
        throw(MibSFailureException("MiBS failed while solving subproblem $(sub_solver.name) for separating BlC cuts."))
    end

    opt_obj = eval_affexpr_by_name(
        sub_solver.c_objterm - k_term,
        sub_solver.bi_model,
        solution.all_upper,
        solution.all_lower,
    )
    @debug "Optimal solution of sub $(sub_solver.name) is $(opt_obj)"

    is_violate = Bool(!(sval > opt_obj - var_non_zero_tolerance))
    @debug "The violated status of sub $(sub_solver.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    r = eval_affexpr_by_name(sub_solver.r_objterm, sub_solver.bi_model, solution.all_upper, solution.all_lower)
    c = eval_affexpr_by_name(sub_solver.c_objterm, sub_solver.bi_model, solution.all_upper, solution.all_lower)
    as = indices_where_coeff_pos(solution.all_upper, sub_solver.link_varsC; non_used=true)
    return SubSolution(is_violate, r, c, opt_obj, as)
end

function set_nthreads(sol::SubSolverMiBS, n)
    @warn "We currently do not directly support setting the number of threads for the MiBS subsolver."
end

function supports_bilevel_subproblem_solver(sol::SubSolverMiBS)
    return true
end

function solve_sub_for_x(sol::SubSolverMiBS, xvals, params::SolverParam, time_limit)
    @objective(Upper(sol.bi_model), Min, sol.c_objterm)

    @constraint(Upper(sol.bi_model), fixc[a=sol.A], sol.link_varsC[a] == round(xvals[a]))

    try
        if should_debbug_print(params)
            write_to_file(
                sol.bi_model,
                "$(params.output_folder_path)/subfix_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    @debug "Subproblem $(sol.name) MIP was adjusted. Start MibS solver to solve it."
    solution = _solve_mip(sol, params.debbug_out, time_limit)

    try
        if !solution.status
            if solution.timeout
                throw(TimeoutException("We hit the time limit while solving MiBS subsolver $(sol.name) for fixed first-level values."))
            end
            throw(MibSFailureException("MiBS failed while solving subproblem $(sol.name) for fixed first-level values."))
        end

        y_vals = index_value_map(sol.y_vars, solution.all_lower)
        osol = eval_affexpr_by_name(sol.c_objterm, sol.bi_model, solution.all_upper, solution.all_lower)
        @debug "We found a solution for subproblem $(sol.name) with value $(osol)."

        try
            if should_debbug_print(params)
                write_to_file(
                    sol.bi_model,
                    "$(params.output_folder_path)/subfixCLEAN_$(sol.name).$(params.file_format_output)",
                )
            end
        catch err
            @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
        end

        osol_L1 = eval_affexpr_by_name(sol.r_objterm, sol.bi_model, solution.all_upper, solution.all_lower)
        return true, osol, osol_L1, y_vals
    finally
        for a in sol.A
            delete(sol.bi_model, fixc[a])
        end
        unregister(sol.bi_model, :fixc)
    end
end

function eval_affexpr_by_name(
    A::GenericAffExpr,
    model::BilevelModel,
    upper_values::AbstractDict{String,<:Real},
    lower_values::AbstractDict{String,<:Real};
    require_all::Bool=true,
)
    s = JuMP.constant(A)
    seen = Set{String}()
    merged_values = Dict{String,Float64}()
    for (nm, val) in upper_values
        merged_values[nm] = float(val)
    end
    for (nm, val) in lower_values
        merged_values[nm] = float(val)
    end

    for (nm, val) in merged_values
        vref = JuMP.variable_by_name(model, nm)
        vref === nothing && continue
        coeff = JuMP.coefficient(A, vref)
        if coeff != 0.0
            s += coeff * val
            push!(seen, nm)
        end
    end

    if require_all
        for v in all_variables(model)
            coeff = JuMP.coefficient(A, v)
            if coeff != 0.0 && !(JuMP.name(v) in seen)
                error("Missing value for variable $(JuMP.name(v)) used in the expression.")
            end
        end
    end
    return s
end

_indexed_refs(B::AbstractDict) = collect(B)
_indexed_refs(B) = [(i, B[i]) for i in axes(B, 1)]

function indices_where_coeff_pos(values::Dict{String,<:Real}, B; non_used=false)
    varnonzero_tolerance = 4e-4
    pairs_B = _indexed_refs(B)
    I = typeof(first(first(pairs_B)))
    result = I[]
    for (i, v) in pairs_B
        nm = JuMP.name(v)
        if !haskey(values, nm)
            error("No value provided for variable '$nm'.")
        end
        val = values[nm]
        if val < -varnonzero_tolerance
            @error "The value of variable $nm was negative in MiBS solution, but variable should be binary!"
        end
        if non_used
            if val < 1 - varnonzero_tolerance
                push!(result, i)
            end
        else
            if val > varnonzero_tolerance
                push!(result, i)
            end
        end
    end
    return result
end

function index_value_map(B, values::AbstractDict{String,<:Real})
    pairs_B = _indexed_refs(B)
    I = typeof(first(first(pairs_B)))
    result = Dict{I,Float64}()
    for (i, v) in pairs_B
        nm = JuMP.name(v)
        haskey(values, nm) || error("No value provided for variable '$nm'.")
        result[i] = float(values[nm])
    end
    return result
end

function _call_mibs_local(
    mps_filename::AbstractString,
    aux_filename::AbstractString,
    runtime_limit::Real,
)
    output_path = joinpath(dirname(mps_filename), "mibs_output.txt")
    error_path = joinpath(dirname(mps_filename), "mibs_errors.txt")
    par_path = joinpath(dirname(mps_filename), "mibs.par")
    _write_mibs_parameter_file!(
        par_path,
        mps_filename,
        aux_filename,
        MibSparam(false, dirname(mps_filename), runtime_limit),
    )
    return _call_mibs_with_param_file(
        par_path,
        output_path,
        error_path,
        runtime_limit,
    )
end

function _parse_mibs_output_local(
    output::AbstractString,
    new_model::MOI.FileFormats.MPS.Model,
    lower_variables::Vector{MOI.VariableIndex},
    lower_name_map::AbstractDict{String,String}=Dict{String,String}(),
    upper_name_map::AbstractDict{String,String}=Dict{String,String}(),
)
    lines = split(output, '\n')
    found_status = false
    found_timeout = false
    objective_value = NaN

    upper = Dict{Int,Float64}()
    lower = Dict{Int,Float64}()

    all_var = MOI.get(new_model, MOI.ListOfVariableIndices())

    cntu = 1
    cntd = 1

    dict_lower_name = Dict{Int,String}()
    dict_lower_value = Dict{String,Float64}()
    dict_upper_name = Dict{Int,String}()
    dict_upper_value = Dict{String,Float64}()
    dict_upper_index_to_model = Dict{Int,MOI.VariableIndex}()
    dict_lower_index_to_model = Dict{Int,MOI.VariableIndex}()
    dict_all = Dict{MOI.VariableIndex,Float64}()

    for y in all_var
        nameofvar = MOI.get(new_model, MOI.VariableName(), y)
        if y in lower_variables
            dict_lower_name[cntd] = nameofvar
            dict_lower_value[get(lower_name_map, nameofvar, nameofvar)] = 0.0
            dict_lower_index_to_model[cntd] = y
            cntd += 1
        else
            dict_upper_name[cntu] = nameofvar
            dict_upper_value[get(upper_name_map, nameofvar, nameofvar)] = 0.0
            dict_upper_index_to_model[cntu] = y
            cntu += 1
        end
        dict_all[y] = 0.0
    end

    for line in lines
        if !found_status
            if occursin("Optimal solution", line)
                found_status = true
            elseif occursin("Reached time limit", line)
                found_timeout = true
            end
            if !found_status
                continue
            end
        end

        m = match(r"(.+)\[([0-9]+)\] \= (.+)", line)
        if m === nothing
            m = match(r"Cost \= (.+)", line)
            if m !== nothing
                objective_value = parse(Float64, m[1])
            end
            continue
        end

        var_prefix = strip(m[1])
        column = parse(Int, m[2])
        value = parse(Float64, m[3])

        if haskey(dict_upper_name, column) && dict_upper_name[column] == "$(var_prefix)[$column]"
            upper[column] = value
            nameofvar = get(upper_name_map, dict_upper_name[column], dict_upper_name[column])
            dict_upper_value[nameofvar] = value
            indexofvar = dict_upper_index_to_model[column]
            dict_all[indexofvar] = value
        elseif haskey(dict_lower_name, column) && dict_lower_name[column] == "$(var_prefix)[$column]"
            lower[column] = value
            nameofvar = get(lower_name_map, dict_lower_name[column], dict_lower_name[column])
            dict_lower_value[nameofvar] = value
            indexofvar = dict_lower_index_to_model[column]
            dict_all[indexofvar] = value
        elseif haskey(dict_upper_name, column)
            upper[column] = value
            nameofvar = get(upper_name_map, dict_upper_name[column], dict_upper_name[column])
            dict_upper_value[nameofvar] = value
            indexofvar = dict_upper_index_to_model[column]
            dict_all[indexofvar] = value
        elseif haskey(dict_lower_name, column)
            lower[column] = value
            nameofvar = get(lower_name_map, dict_lower_name[column], dict_lower_name[column])
            dict_lower_value[nameofvar] = value
            indexofvar = dict_lower_index_to_model[column]
            dict_all[indexofvar] = value
        else
            error("Could not map MiBS output line '$(line)' to an upper or lower variable.")
        end
    end

    return (
        status=found_status,
        timeout=found_timeout,
        objective=objective_value,
        nonzero_upper=upper,
        nonzero_lower=lower,
        all_upper=dict_upper_value,
        all_lower=dict_lower_value,
        all_var=dict_all,
        log_output=output,
    )
end

function _solve_mip(sol::SubSolverMiBS, printlog::Bool, time_limit::Real)
    @debug "Starting MiBS subsolver $(sol.name)"
    path = repo_local_tempdir("mibs_subsolver", _mibs_safe_name(sol.name); prefix="mibs_sub")
    try
        mps_filename = joinpath(path, "model.mps")
        aux_filename = joinpath(path, "model.aux")
        new_model, variables, objective, constraints, sense =
            BilevelJuMP._build_single_model(sol.bi_model, true)
        _initialize_name_maps!(sol)
        MOI.write_to_file(new_model, mps_filename)
        BilevelJuMP._write_auxillary_file(
            new_model,
            variables,
            objective,
            constraints,
            sense,
            aux_filename,
        )

        output, err = _call_mibs_local(mps_filename, aux_filename, time_limit)
        if !isempty(err)
            if occursin("MiBS exceeded the JuBiC wall-clock limit", err)
                throw(TimeoutException("We hit the time limit while solving MiBS subsolver $(sol.name)."))
            end
            throw(MibSFailureException("MibS returned:\n\n$(err)"))
        end
        if isempty(output)
            error("MibS failed to return output.")
        end
        if printlog
            print(output)
            debug_dir = repo_local_tempdir("mibs_subsolver_logs", _mibs_safe_name(sol.name); prefix="mibs_log")
            cp(mps_filename, joinpath(debug_dir, "model.mps"); force=true)
            cp(aux_filename, joinpath(debug_dir, "model.aux"); force=true)
            write(joinpath(debug_dir, "mibs_output.txt"), output)
            write(joinpath(debug_dir, "mibs_errors.txt"), err)
            @info "Stored MiBS subsolver debug files for $(sol.name) in $(debug_dir)."
        end

        return _parse_mibs_output_local(
            output,
            new_model,
            variables,
            sol.lower_name_map,
            sol.upper_name_map,
        )
    finally
        rm(path; force=true, recursive=true)
    end
end
