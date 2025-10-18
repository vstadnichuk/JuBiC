struct AUXData
    variables::Set{String}  # Set of variable names
    constraints::Set{String}  # Set of constraint names
    objective::Dict{String,Number}  # Dict mapping variable names to objective coefficients
end

"""
    _read_aux(file_path::String) -> AUXData

Reads an aux file and returns an AUXData struct containing the variables, constraints, and objective coefficients.
"""
function _read_aux(file_path::String)
    # TODO: Add support for Legacy Index-based aux files 
    variables = Set{String}()
    constraints = Set{String}()
    objective = Dict{String,Number}()

    current_section = ""
    open(file_path, "r") do file
        for line in eachline(file)
            line = strip(line)
            if isempty(line) || startswith(line, "*")
                continue
            end
            if startswith(line, "@")
                current_section = line
                continue
            end

            parts = split(line)
            if current_section == "@VARSBEGIN"
                if length(parts) < 2
                    error("Invalid variable line in aux file: $line")
                end
                push!(variables, parts[1])
                objective[parts[1]] = parse(Float64, parts[2])
            elseif current_section == "@CONSTRSBEGIN"
                if length(parts) < 1
                    error("Invalid constraint line in aux file: $line")
                end
                push!(constraints, parts[1])
            end
        end
    end

    return AUXData(variables, constraints, objective)
end
