mutable struct VariableBound
    is_binary::Bool  # True if variable is binary (BV)
    is_integer::Bool  # True if variable is integer (LI, UI)
    lower_bound::Union{Nothing,Float64}  # Lower bound value or nothing if not set
    upper_bound::Union{Nothing,Float64}  # Upper bound value or nothing if not set
end

struct MPSData
    name::String  # Name of the model
    rows_less_than::Vector{String}  # Row names with 'L' (<=) type
    rows_greater_than::Vector{String}  # Row names with 'G' (>=) type
    rows_equal::Vector{String}  # Row names with 'E' (=) type
    rows_natural::Vector{String}  # Row names with 'N' (objective) type
    columns::Dict{String,Vector{Tuple{String,Number}}}  # Mapping row names to list of (variable name, coefficient)
    rhs::Dict{String,Number}  # Mapping row names to their RHS value
    bounds::Dict{String,VariableBound}  # Mapping variable names to their bounds
end

"""
    _read_mps(file_path::String) -> MPSData

Reads an mps file and returns an MPSData struct containing the parsed data.
"""
function _read_mps(file_path::String)
    SECTION_NAMES = Set(["NAME", "ROWS", "COLUMNS", "RHS", "BOUNDS", "ENDATA"])

    name = ""
    rows_less_than = String[]
    rows_greater_than = String[]
    rows_equal = String[]
    rows_natural = String[]
    columns = Dict{String,Vector{Tuple{String,Number}}}()
    rhs = Dict{String,Number}()
    bounds = Dict{String,VariableBound}()

    current_section = ""
    open(file_path, "r") do file
        for line in eachline(file)
            line = strip(line)
            if isempty(line) || startswith(line, "*")
                continue
            end

            parts = split(line)
            if isempty(parts)
                continue
            end

            if parts[1] in SECTION_NAMES
                current_section = parts[1]
                if current_section == "NAME" && length(parts) >= 2
                    name = parts[2]
                elseif current_section == "ENDATA"
                    break
                end
                continue
            end

            try
                if current_section == "ROWS"
                    if length(parts) >= 2
                        row_type = parts[1][1]
                        row_name = parts[2]
                        if row_type == 'N'
                            push!(rows_natural, row_name)
                        elseif row_type == 'L'
                            push!(rows_less_than, row_name)
                        elseif row_type == 'G'
                            push!(rows_greater_than, row_name)
                        elseif row_type == 'E'
                            push!(rows_equal, row_name)
                        end
                    end
                elseif current_section == "COLUMNS"
                    if parts[1] in ("MARK0000", "MARK0001")
                        continue
                    end
                    var_name = parts[1]
                    row_name = parts[2]
                    coeff = parse(Float64, parts[3])
                    if !haskey(columns, row_name)
                        columns[row_name] = Vector{Tuple{String,Number}}()
                    end
                    push!(columns[row_name], (var_name, coeff))

                    if length(parts) > 3
                        row_name2 = parts[4]
                        coeff2 = parse(Float64, parts[5])
                        if !haskey(columns, row_name2)
                            columns[row_name2] = Vector{Tuple{String,Number}}()
                        end
                        push!(columns[row_name2], (var_name, coeff2))
                    end

                    if !haskey(rhs, row_name)
                        rhs[row_name] = 0.0
                    end
                elseif current_section == "RHS"
                    row_name = parts[2]
                    rhs[row_name] = parse(Float64, parts[3])
                elseif current_section == "BOUNDS"
                    is_binary = false
                    is_integer = false
                    lower_bound = nothing
                    upper_bound = nothing

                    bound_type = parts[1]
                    var_name = parts[3]
                    if bound_type == "LO"
                        lower_bound = parse(Float64, parts[4])
                    elseif bound_type == "UP"
                        upper_bound = parse(Float64, parts[4])
                    elseif bound_type == "FX"
                        lower_bound = parse(Float64, parts[4])
                        upper_bound = parse(Float64, parts[4])
                    elseif bound_type == "BV"
                        is_binary = true
                    elseif bound_type == "LI"
                        is_integer = true
                        lower_bound = parse(Float64, parts[4])
                    elseif bound_type == "UI"
                        is_integer = true
                        upper_bound = parse(Float64, parts[4])
                    elseif bound_type == "SC"
                        error("Semi-continuous variables (SC) are not supported.")
                    end

                    if haskey(bounds, var_name)
                        bounds[var_name].is_binary |= is_binary
                        bounds[var_name].is_integer |= is_integer
                        if lower_bound !== nothing
                            bounds[var_name].lower_bound = lower_bound
                        end
                        if upper_bound !== nothing
                            bounds[var_name].upper_bound = upper_bound
                        end
                    else
                        bounds[var_name] = VariableBound(is_binary, is_integer, lower_bound, upper_bound)
                    end
                end
            catch e
                error("Error parsing line in section $current_section: $line\n$e")
            end
        end
    end

    return MPSData(name, rows_less_than, rows_greater_than, rows_equal, rows_natural, columns, rhs, bounds)
end
