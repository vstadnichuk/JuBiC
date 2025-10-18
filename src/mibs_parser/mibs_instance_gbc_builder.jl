const VARIABLE_AUX_SUFFIX = "_JuBiCaux"

"""
    get_GBC_instance(mps_file_path::String, aux_file_path::String, optimizer) -> Instance

Builds an instance from the given MPS and AUX files.
"""
function get_GBC_instance(mps_file_path::String, aux_file_path::String, optimizer)
    mps_data = _read_mps(mps_file_path)
    aux_data = _read_aux(aux_file_path)

    nsub = "Sub"
    upper_vars = Dict{String,Any}()
    lower_vars = Dict{String,Any}()
    model_upper = Model(optimizer)
    model_lower = Model(optimizer)

    # create the variables
    for (name, bound) in mps_data.bounds
        sub_model = (name in aux_data.variables) ? model_lower : model_upper

        lb = isnothing(bound.lower_bound) ? -Inf : bound.lower_bound
        ub = isnothing(bound.upper_bound) ? Inf : bound.upper_bound

        var_ref = if bound.is_binary
            @variable(sub_model, base_name = name, binary = true)
        elseif bound.is_integer
            @variable(sub_model, base_name = name, integer = true, lower_bound = lb, upper_bound = ub)
        else
            @variable(sub_model, base_name = name, lower_bound = lb, upper_bound = ub)
        end

        if name in aux_data.variables
            lower_vars[name] = var_ref
        else
            upper_vars[name] = var_ref
        end
    end

    function _is_var_in_lower_model(var_name::String)
        for row in aux_data.constraints
            if any(var_name == name for (name, _) in mps_data.columns[row])
                return true
            end
        end
        return false
    end

    link_vars = [v for v in keys(upper_vars) if _is_var_in_lower_model(v)]
    link_vars_aux = string.(link_vars, VARIABLE_AUX_SUFFIX)
    A = [link_vars; link_vars_aux]

    # check whether all linking variables are binary
    for var in link_vars
        @assert is_binary(upper_vars[var]) "Linking variable $var is not binary"
    end

    @variable(model_upper, link_upper_aux[link_vars_aux], Bin)
    @variable(model_lower, link_lower[A], Bin)
    @variable(model_lower, upper_copy[A], Bin)

    # linking constraints
    @constraint(model_upper, [a in link_vars],
        link_upper_aux[string(a, VARIABLE_AUX_SUFFIX)] + upper_vars[a] == 1)
    @constraint(model_lower, [a in link_vars],
        link_lower[string(a, VARIABLE_AUX_SUFFIX)] + link_lower[a] == 1)
    @constraint(model_lower, [a in A], link_lower[a] <= upper_copy[a])

    xdict = merge(
        Dict(a => upper_vars[a] for a in link_vars),
        Dict(a => link_upper_aux[a] for a in link_vars_aux)
    )

    # objectives
    lower_obj = @expression(model_lower, sum(aux_data.objective[name] * var for (name, var) in lower_vars))
    @objective(model_lower, Min, 1 * lower_obj)
    @objective(model_upper, Min,
        sum(upper_vars[name] * n for (name, n) in mps_data.columns[mps_data.rows_natural[1]] if name in keys(upper_vars)))
    master_sub_obj = sum(lower_vars[name] * n for (name, n) in mps_data.columns[mps_data.rows_natural[1]] if name in keys(lower_vars))

    # constraints
    function _add_constraints!(row_names::Vector{String}, sense::Char)
        for row in row_names
            sub_model = (row in aux_data.constraints ? model_lower : model_upper)
            lhs = if sub_model == model_upper
                sum(upper_vars[name] * n for (name, n) in mps_data.columns[row]; init=0)
            else
                sum(lower_vars[name] * n for (name, n) in mps_data.columns[row] if name in keys(lower_vars); init=0) +
                sum(link_lower[name] * n for (name, n) in mps_data.columns[row] if name in keys(upper_vars); init=0)
            end
            rhs = mps_data.rhs[row]
            if sense == 'L'
                @constraint(sub_model, lhs <= rhs)
            elseif sense == 'G'
                @constraint(sub_model, lhs >= rhs)
            elseif sense == 'E'
                @constraint(sub_model, lhs == rhs)
            end
        end
    end

    _add_constraints!(mps_data.rows_less_than, 'L')
    _add_constraints!(mps_data.rows_greater_than, 'G')
    _add_constraints!(mps_data.rows_equal, 'E')

    master = Master(model_upper, A, xdict, [nsub])
    link_constraints_capacities = Dict(a => 1.0 for a in A)

    subS = SubSolverJuMP(
        nsub,
        model_lower,
        A,
        upper_copy,
        link_lower,
        link_constraints_capacities,
        master_sub_obj,
        lower_obj,
        timelimit -> (false, 0)
    )

    set_silent(model_lower)
    return Instance(master, [subS])
end
