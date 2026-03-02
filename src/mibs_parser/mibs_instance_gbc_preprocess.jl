const LINKING_CONSTRAINT_SUFFIX = "_JuBiC-negation"


"""
    _is_var_in_lower_model(var_name::String, mps_data::MPSData, aux_data::AUXData) -> Bool

Check if a variable is present in the lower model based on the auxiliary data.
"""
function _is_var_in_lower_model(var_name::String, mps_data::MPSData, aux_data::AUXData)
    return any(row -> any(var_name == name for (name, _) in mps_data.columns[row]), aux_data.constraints)
end

"""
    _is_var_in_lower_model_without_linking(var_name::String, mps_data::MPSData, aux_data::AUXData) -> Bool

Check if a variable is present in the lower model except a single linking constraint.
"""
function _is_var_in_lower_model_without_linking(var_name::String, mps_data::MPSData, aux_data::AUXData)
    other_var_name = nothing
    C = 0.0

    for row in aux_data.constraints
        if any(var_name == name for (name, _) in mps_data.columns[row])
            if other_var_name !== nothing
                return true, nothing
            end

            if mps_data.rhs[row] != 0
                return true, nothing
            end

            if !(row in mps_data.rows_less_than)
                return true, nothing
            end

            if length(mps_data.columns[row]) != 2
                return true, nothing
            end

            column = mps_data.columns[row]
            factor = column[1][1] == var_name ? column[1][2] : column[2][2]
            other_factor = column[1][1] == var_name ? column[2][2] : column[1][2]
            other_var_name = column[1][1] == var_name ? column[2][1] : column[1][1]
            C = -factor / other_factor

            if other_factor <= 0 || C < 0
                return true, nothing
            end
        end
    end

    return false, (other_var_name, C)
end

"""
    _classify_variables(mps_data::MPSData, aux_data::AUXData, preprocessing::Bool) -> Tuple{Vector{String}, Dict{String,Tuple{String,Number}}, Vector{String}}

Classifies variables into regular linking variables and link variables based on their presence in the lower model and preprocessing flag.
"""
function _classify_variables(mps_data::MPSData, aux_data::AUXData, preprocessing::Bool)
    regular_linking = String[]
    regular_linking_data = Dict{String,Tuple{String,Number}}()
    link_vars = String[]

    for var_name in keys(mps_data.bounds)
        if !(var_name in aux_data.variables) && _is_var_in_lower_model(var_name, mps_data, aux_data)
            if preprocessing
                not_binary_linking, data = _is_var_in_lower_model_without_linking(var_name, mps_data, aux_data)
                if not_binary_linking
                    push!(link_vars, var_name)
                else
                    other_var_name, coeff = data
                    regular_linking_data[var_name] = (other_var_name, coeff)
                    push!(regular_linking, var_name)
                end
            else
                push!(link_vars, var_name)
            end
        end
    end

    return regular_linking, regular_linking_data, link_vars
end

"""
    _preprocess_model(mps_data::MPSData, aux_data::AUXData) -> Float64

Preprocesses the MPS model by transforming certain constraints and variables. It returns the constant added to the objective function due to preprocessing.
"""
function _preprocess_model(mps_data::MPSData, aux_data::AUXData)
    obj_bias = 0.0

    # transform a*x + y <= a into -a*\bar{x} + y <= 0
    for (var_name, bound) in mps_data.bounds
        if bound.is_binary && !(var_name in aux_data.variables)
            is_binary_linking = false

            for row in aux_data.constraints
                if any(var_name == name for (name, _) in mps_data.columns[row])
                    if is_binary_linking
                        is_binary_linking = false
                        break
                    end

                    if !(row in mps_data.rows_less_than)
                        break
                    end

                    if length(mps_data.columns[row]) != 2
                        break
                    end

                    column = mps_data.columns[row]
                    factor = (column[1][1] == var_name) ? column[1][2] : column[2][2]
                    if mps_data.rhs[row] <= 0 || factor <= 0 || factor != mps_data.rhs[row]
                        break
                    end

                    is_binary_linking = true
                end
            end

            if is_binary_linking
                delete!(mps_data.bounds, var_name)

                new_var_name = var_name * LINKING_CONSTRAINT_SUFFIX
                mps_data.bounds[new_var_name] = VariableBound(true, false, 0.0, 1.0)

                for (column, entries) in mps_data.columns
                    idx = findfirst(c -> c[1] == var_name, entries)
                    if idx !== nothing
                        _, coefficient = entries[idx]

                        filter!(e -> e[1] != var_name, mps_data.columns[column])
                        push!(mps_data.columns[column], (new_var_name, -coefficient))

                        if column in mps_data.rows_natural
                            obj_bias += coefficient
                        else
                            mps_data.rhs[column] -= coefficient
                        end
                    end
                end
            end
        end
    end

    # TODO: add other preprocessing steps

    return obj_bias
end
