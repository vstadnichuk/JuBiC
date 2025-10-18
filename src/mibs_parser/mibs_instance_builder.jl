using BilevelJuMP

"""
    get_MibS_instance(mps_file_path::String, aux_file_path::String) -> Instance

Builds a bilevel optimization instance from the given MPS and AUX files.
"""
function get_MibS_instance(mps_file_path::String, aux_file_path::String)
    mps_data = _read_mps(mps_file_path)
    aux_data = _read_aux(aux_file_path)

    model = BilevelModel()

    lower_vars = Dict{String,Any}()
    upper_vars = Dict{String,Any}()
    all_vars = Dict{String,Any}()

    # create variables with bounds
    for (name, bound) in mps_data.bounds
        sub_model = name in aux_data.variables ? Lower(model) : Upper(model)
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
        all_vars[name] = var_ref
    end

    # lower-level objective
    @objective(
        Lower(model),
        Min,
        sum(aux_data.objective[name] * var for (name, var) in lower_vars)
    )

    # upper-level objective
    @objective(
        Upper(model),
        Min,
        sum(all_vars[name] * n for (name, n) in mps_data.columns[mps_data.rows_natural[1]])
    )

    # add constraints based on sense
    function _add_constraints!(row_names::Vector{String}, sense::Char)
        for row in row_names
            sub_model = (row in aux_data.constraints ? Lower(model) : Upper(model))
            lhs = sum(all_vars[name] * n for (name, n) in mps_data.columns[row])
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

    return Instance(MibSMaster(model), nothing)
end
