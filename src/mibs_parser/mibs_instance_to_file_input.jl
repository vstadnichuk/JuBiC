using JuMP, BilevelJuMP

struct ConstraintData
    index::Int
    is_upper::Bool
    sense::Char  # 'L', 'G', 'E'
    vars::Dict{String,Float64}
    rhs::Float64
end

"""
    _get_variable_names(model::BilevelJuMP.BilevelModel) -> (Dict{String,String}, Vector{JuMP.VariableRef}, Vector{JuMP.VariableRef})

Returns a mapping from original variable names to new standardized names, along with lists of upper and lower level variable references.
"""
function _get_variable_names(model::BilevelJuMP.BilevelModel, anonymous_names::Bool)
    var_names = Dict{String,String}()
    upper_vars = Vector{JuMP.VariableRef}()
    lower_vars = Vector{JuMP.VariableRef}()
    lower_var_count = 0
    upper_var_count = 0

    for (i, info) in model.var_info
        if info.level in [BilevelJuMP.LOWER_BOTH, BilevelJuMP.LOWER_ONLY]
            lower_var_count += 1

            var_name = JuMP.name(model.var_lower[i])
            new_var_name = anonymous_names ? "y" * lpad(string(lower_var_count), 5, '0') : var_name

            push!(lower_vars, model.var_lower[i])
            var_names[var_name] = new_var_name
        else
            upper_var_count += 1

            var_name = JuMP.name(model.var_upper[i])
            new_var_name = anonymous_names ? "x" * lpad(string(upper_var_count), 5, '0') : var_name

            push!(upper_vars, model.var_upper[i])
            var_names[var_name] = new_var_name
        end
    end

    return var_names, upper_vars, lower_vars
end

"""
    _get_constraint_data(model::BilevelJuMP.BilevelModel, var_names::Dict{String,String}) -> Vector{ConstraintData}

Returns a list of data for all constraints in the given BilevelModel.
"""
function _get_constraint_data(model::BilevelJuMP.BilevelModel, var_names::Dict{String,String})
    upper_constraints = JuMP.all_constraints(model.upper; include_variable_in_set_constraints=false)
    lower_constraints = JuMP.all_constraints(model.lower; include_variable_in_set_constraints=false)
    all_constraints = vcat(upper_constraints, lower_constraints)

    constraint_data = Vector{ConstraintData}()

    for (index, creff) in enumerate(all_constraints)
        is_upper = creff in upper_constraints
        obj = JuMP.constraint_object(creff)
        set = obj.set

        sense = if set isa MOI.LessThan
            'L'
        elseif set isa MOI.GreaterThan
            'G'
        elseif set isa MOI.EqualTo
            'E'
        else
            'E'
        end

        rhs = if set isa MOI.LessThan
            set.upper
        elseif set isa MOI.GreaterThan
            set.lower
        elseif set isa MOI.EqualTo
            set.value
        else
            0.0
        end

        vars = Dict{String,Float64}()
        for (var, coeff) in obj.func.terms
            var_name = var_names[JuMP.name(var)]
            vars[var_name] = coeff
        end

        push!(constraint_data, ConstraintData(index, is_upper, sense, vars, rhs))
    end

    return constraint_data
end

_format_mps_field(name1::String, name2::String, value) = "    " * rpad(string(name1), 20) * rpad(string(name2), 20) * lpad(string(value), 23) * "\n"
_format_row_name(index::Int) = "R" * lpad(string(index), 5, '0')

"""
    _write_mps_file(io::IO, name::String, constraint_data::Vector{ConstraintData}, upper_obj_terms, model, var_names::Dict{String,String}, upper_vars::Vector{JuMP.VariableRef}, lower_vars::Vector{JuMP.VariableRef})

Writes the MPS file representation of the given model to the provided IO stream.
"""
function _write_mps_file(io::IO, name::String, constraint_data::Vector{ConstraintData}, upper_obj_terms, model, var_names::Dict{String,String}, upper_vars::Vector{JuMP.VariableRef}, lower_vars::Vector{JuMP.VariableRef})
    write(io, "NAME          $(name)\n")
    write(io, "ROWS\n")
    write(io, " N  OBJ\n")
    for constr in constraint_data
        write(io, " $(constr.sense)  $(_format_row_name(constr.index))\n")
    end

    write(io, "COLUMNS\n")
    upper_obj = Dict(JuMP.name(k) => v for (k, v) in upper_obj_terms)
    for var in [lower_vars; upper_vars]
        var_name = var_names[JuMP.name(var)]

        if haskey(upper_obj, JuMP.name(var))
            coeff = upper_obj[JuMP.name(var)]
            write(io, _format_mps_field(var_name, "OBJ", coeff))
        end

        for constr in constraint_data
            if haskey(constr.vars, var_name)
                coeff = constr.vars[var_name]
                write(io, _format_mps_field(var_name, _format_row_name(constr.index), coeff))
            end
        end
    end

    write(io, "RHS\n")
    for constr in constraint_data
        if constr.rhs != 0.0
            write(io, _format_mps_field("rhs", _format_row_name(constr.index), constr.rhs))
        end
    end

    write(io, "BOUNDS\n")
    for var in JuMP.all_variables(model)
        var_name = var_names[JuMP.name(var)]
        is_binary = JuMP.is_binary(var)
        is_integer = JuMP.is_integer(var)
        lb = has_lower_bound(var) ? lower_bound(var) : -Inf
        ub = has_upper_bound(var) ? upper_bound(var) : Inf

        if is_binary
            write(io, " BV bnd $(repeat(" ", 15)) $(var_name)\n")
        elseif is_integer
            if isfinite(lb)
                write(io, " LI bnd $(repeat(" ", 15)) $(var_name) $(lb)\n")
            end
            if isfinite(ub)
                write(io, " UI bnd $(repeat(" ", 15)) $(var_name) $(ub)\n")
            end
            if !isfinite(lb) && !isfinite(ub)
                inf = 1e9
                write(io, " LI bnd $(repeat(" ", 15)) $(var_name) $(-inf)\n")
                write(io, " UI bnd $(repeat(" ", 15)) $(var_name) $(inf)\n")
            end
        else
            if isfinite(lb)
                write(io, " LO bnd $(repeat(" ", 15)) $(var_name) $(lb)\n")
            end
            if isfinite(ub)
                write(io, " UP bnd $(repeat(" ", 15)) $(var_name) $(ub)\n")
            end
            if !isfinite(lb) && !isfinite(ub)
                write(io, " FR bnd $(repeat(" ", 15)) $(var_name)\n")
            end
        end
    end

    write(io, "ENDATA")
end

"""
    _write_aux_file(io::IO, name::String, number_of_lower_vars::Int, number_of_lower_constrs::Int, constraint_data::Vector{ConstraintData}, lower_obj_terms, var_names::Dict{String,String}, lower_vars::Vector{JuMP.VariableRef})

Writes the AUX file representation of the given model to the provided IO stream.
"""
function _write_aux_file(io::IO, name::String, number_of_lower_vars::Int, number_of_lower_constrs::Int, constraint_data::Vector{ConstraintData}, lower_obj_terms, var_names::Dict{String,String}, lower_vars::Vector{JuMP.VariableRef})
    write(io, "@NAME\n$(name)\n")
    write(io, "@MPS\n$(name).mps\n")
    write(io, "@NUMVARS\n$(number_of_lower_vars)\n")
    write(io, "@NUMCONSTRS\n$(number_of_lower_constrs)\n")

    write(io, "@VARSBEGIN\n")
    for var in lower_vars
        var_name = var_names[JuMP.name(var)]
        coeff = haskey(lower_obj_terms, var) ? lower_obj_terms[var] : 0.0
        write(io, "$(var_name) $(coeff)\n")
    end
    write(io, "@VARSEND\n")

    write(io, "@CONSTRSBEGIN\n")
    for constr in constraint_data
        if !constr.is_upper
            write(io, _format_row_name(constr.index) * "\n")
        end
    end
    write(io, "@CONSTRSEND")
end

"""
    output_MibS_instance(instance::Instance, name::String, output_directory::String)

Outputs the given MibS instance to the MibS input file format (MPS and AUX files).
"""
function output_MibS_instance(instance::Instance, name::String, output_directory::String; anonymous_names::Bool=true)
    @assert !isnothing(instance.master) "The instance must have a master problem defined."
    @assert instance.master isa MibSMaster "The master problem must be of type MibSMaster."
    @assert !isnothing(instance.master.model) "The MibSMaster must have a BilevelModel defined."

    if !isdir(output_directory)
        mkpath(output_directory)
    end

    mps_file_path = joinpath(output_directory, name * ".mps")
    aux_file_path = joinpath(output_directory, name * ".aux")

    model = instance.master.model
    var_names, upper_vars, lower_vars = _get_variable_names(model, anonymous_names)
    constraint_data = _get_constraint_data(model, var_names)

    number_of_lower_constrs = length(JuMP.all_constraints(model.lower; include_variable_in_set_constraints=false))
    number_of_lower_vars = length(model.lower_only) + length(model.lower_to_upper_link)

    upper_obj_terms = JuMP.objective_function(model.upper).terms
    lower_obj_terms = JuMP.objective_function(model.lower).terms

    open(mps_file_path, "w") do mps_io
        _write_mps_file(mps_io, name, constraint_data, upper_obj_terms, model, var_names, upper_vars, lower_vars)
    end

    open(aux_file_path, "w") do aux_io
        _write_aux_file(aux_io, name, number_of_lower_vars, number_of_lower_constrs, constraint_data, lower_obj_terms, var_names, lower_vars)
    end
end

"""
    output_GBC_instance(instance::Instance, name::String, output_directory::String)

Outputs the given GBC instance to the MibS input file format (MPS and AUX files).
"""
function output_GBC_instance(instance::Instance, name::String, output_directory::String; anonymous_names::Bool=true)
    @assert !isnothing(instance.master) "The instance must have a master problem defined."
    @assert instance.master isa Master "The master problem must be of type Master."

    mibs_instance = transform_GBC_to_MibS(instance)
    output_MibS_instance(mibs_instance, name, output_directory; anonymous_names=anonymous_names)
end
