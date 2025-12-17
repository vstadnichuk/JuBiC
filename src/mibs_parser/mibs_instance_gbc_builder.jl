const VARIABLE_AUX_SUFFIX = "_JuBiC-aux"
const VARIABLE_PARTIAL_SUFFIX = "_JuBiC-partial"


"""
    get_GBC_instance(mps_file_path::String, aux_file_path::String, optimizer) -> Instance

Builds an instance from the given MPS and AUX files.
"""
function get_GBC_instance(mps_file_path::String, aux_file_path::String, optimizer; partial_decomposition::Bool=true, preprocessing::Bool=true)
    mps_data = _read_mps(mps_file_path)
    aux_data = _read_aux(aux_file_path)

    obj_bias = preprocessing ? _preprocess_model(mps_data, aux_data) : 0.0

    # initialize models and variable containers
    nsub = "Sub"
    model_upper = Model(optimizer)
    model_lower = Model(optimizer)
    upper_vars = Dict{String,Any}()
    lower_vars = Dict{String,Any}()
    upper_vars_partial = Dict{String,Any}()

    # create variables
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
            upper_vars_partial[name] = var_ref
        end

        if partial_decomposition && sub_model === model_lower
            if bound.is_binary
                lb = 0.0
                ub = 1.0
            end
            var_partial = @variable(model_upper, base_name = name * VARIABLE_PARTIAL_SUFFIX, lower_bound = lb, upper_bound = ub)
            upper_vars_partial[name] = var_partial
        end
    end

    # classify and create linking variables
    regular_linking, regular_linking_data, link_vars = _classify_variables(mps_data, aux_data, preprocessing)
    link_vars_aux = String[var_name * VARIABLE_AUX_SUFFIX for var_name in link_vars]
    all_link_vars = [link_vars; link_vars_aux]
    A = [regular_linking; link_vars; link_vars_aux]

    # check whether all linking variables are binary
    for var in link_vars
        @assert is_binary(upper_vars[var]) "Linking variable $var is not binary"
    end

    # create auxiliary linking variables
    @variable(model_upper, link_upper_aux[link_vars_aux], Bin)
    @variable(model_lower, link_lower[all_link_vars], Bin)
    @variable(model_lower, upper_copy[A], Bin)

    # build the y variables (lower level and linked variables)
    y_vars = JuMP.Containers.DenseAxisArray{VariableRef}(undef, A)
    for var_name in A
        y_vars[var_name] = if var_name in regular_linking
            lower_vars[regular_linking_data[var_name][1]]
        else
            link_lower[var_name]
        end
    end

    # add linking constraints
    @constraint(model_upper, [a in link_vars], link_upper_aux[a*VARIABLE_AUX_SUFFIX] + upper_vars[a] == 1)
    @constraint(model_lower, [a in link_vars], link_lower[a*VARIABLE_AUX_SUFFIX] + link_lower[a] == 1)
    @constraint(model_lower, [a in all_link_vars], link_lower[a] <= upper_copy[a])

    # create dictionary for mapping each resource to its upper level variable
    xdict = merge(
        Dict{String,VariableRef}(a => upper_vars[a] for a in regular_linking),
        Dict{String,VariableRef}(a => upper_vars[a] for a in link_vars),
        Dict{String,VariableRef}(a => link_upper_aux[a] for a in link_vars_aux)
    )

    # set lower level objective function
    lower_obj = @expression(model_lower, sum(aux_data.objective[name] * var for (name, var) in lower_vars; init=0))
    @objective(model_lower, Min, 1 * lower_obj)

    # set upper level objective function
    @assert length(mps_data.rows_natural) <= 1 "Multiple natural rows are not supported"
    if length(mps_data.rows_natural) == 1 && haskey(mps_data.columns, mps_data.rows_natural[1])
        @objective(model_upper, Min, sum(upper_vars[name] * n for (name, n) in mps_data.columns[mps_data.rows_natural[1]] if name in keys(upper_vars); init=obj_bias))
        master_sub_obj = sum(lower_vars[name] * n for (name, n) in mps_data.columns[mps_data.rows_natural[1]] if name in keys(lower_vars); init=AffExpr(0))
    else
        @objective(model_upper, Min, obj_bias)
        master_sub_obj = AffExpr(0)
    end

    # add constraints
    function _add_constraints!(row_names::Vector{String}, sense::Char)
        for row in row_names
            sub_model = (row in aux_data.constraints) ? model_lower : model_upper

            lhs = if sub_model === model_upper
                sum(upper_vars[name] * n for (name, n) in mps_data.columns[row])
            else
                sum(lower_vars[name] * n for (name, n) in mps_data.columns[row] if name in keys(lower_vars); init=0) +
                sum(link_lower[name] * n for (name, n) in mps_data.columns[row] if name in keys(upper_vars) && !(name in regular_linking); init=0) +
                sum(upper_copy[name] * n for (name, n) in mps_data.columns[row] if name in keys(upper_vars) && name in regular_linking; init=0)
            end

            rhs = mps_data.rhs[row]

            if sense == 'L'
                @constraint(sub_model, lhs <= rhs)
            elseif sense == 'G'
                @constraint(sub_model, lhs >= rhs)
            elseif sense == 'E'
                @constraint(sub_model, lhs == rhs)
            end

            if partial_decomposition && sub_model === model_lower
                lhs_partial = sum(upper_vars_partial[name] * n for (name, n) in mps_data.columns[row])
                if sense == 'L'
                    @constraint(model_upper, lhs_partial <= rhs)
                elseif sense == 'G'
                    @constraint(model_upper, lhs_partial >= rhs)
                elseif sense == 'E'
                    @constraint(model_upper, lhs_partial == rhs)
                end
            end
        end
    end

    _add_constraints!(mps_data.rows_less_than, 'L')
    _add_constraints!(mps_data.rows_greater_than, 'G')
    _add_constraints!(mps_data.rows_equal, 'E')

    # build and return the instance
    master = Master(model_upper, A, xdict, [nsub])

    link_constraints_capacities = merge(
        Dict{String,Float64}(var => regular_linking_data[var][2] for var in regular_linking),
        Dict{String,Float64}(a => 1.0 for a in all_link_vars)
    )

    subS = SubSolverJuMP(
        nsub,
        model_lower,
        A,
        upper_copy,
        y_vars,
        link_constraints_capacities,
        master_sub_obj,
        lower_obj,
        timelimit -> (false, 0)
    )

    set_silent(model_lower)
    return Instance(master, [subS])
end
