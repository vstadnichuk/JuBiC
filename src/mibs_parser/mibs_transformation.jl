using JuMP, BilevelJuMP

"""
    transform_GBC_to_MibS(instance::Instance) -> Instance

Transforms a GBC instance into a MiBS instance.
"""
function transform_GBC_to_MibS(instance::Instance)
    @assert length(instance.subproblems) == 1 "Model is not a bilevel problem (there are $(length(instance.subproblems)) != 1 subproblems)."
    @assert instance.master isa Master "Model is not a bilevel problem (master is not of type Master)."
    @assert instance.subproblems[1] isa SubSolverJuMP "The transformation is currently only implemented for SubSolverJuMP subproblems."  # TODO: Add support for different subproblem types

    bilevel_model = BilevelJuMP.BilevelModel()
    upper_vars_ref = Dict{Any,Any}()
    lower_vars_ref = Dict{Any,Any}()

    function _create_variable!(model, var)
        lb = has_lower_bound(var) ? lower_bound(var) : -Inf
        ub = has_upper_bound(var) ? upper_bound(var) : Inf
        var_name = JuMP.name(var)

        if is_binary(var)
            return @variable(model, base_name = var_name, binary = true)
        elseif is_integer(var)
            return @variable(model, base_name = var_name, integer = true, lower_bound = lb, upper_bound = ub)
        else
            return @variable(model, base_name = var_name, lower_bound = lb, upper_bound = ub)
        end
    end

    function _add_constraints!(model, src_model, vars_ref)
        constraint_types = [
            (MOI.LessThan{Float64}, (con) -> constraint_object(con).set.upper, 'L'),
            (MOI.GreaterThan{Float64}, (con) -> constraint_object(con).set.lower, 'G'),
            (MOI.EqualTo{Float64}, (con) -> constraint_object(con).set.value, 'E')
        ]

        for (con_type, get_rhs, op) in constraint_types
            for con in all_constraints(src_model, AffExpr, con_type)
                expr = sum(coef * vars_ref[var] for (var, coef) in constraint_object(con).func.terms)
                rhs = get_rhs(con)

                if op == 'E'  # handle equality constraints by adding two inequalities
                    @constraint(model, expr <= rhs)
                    @constraint(model, expr >= rhs)
                elseif op == 'L'
                    @constraint(model, expr <= rhs)
                elseif op == 'G'
                    @constraint(model, expr >= rhs)
                end
            end
        end
    end

    # create upper level variables
    for var in all_variables(instance.master.model)
        upper_vars_ref[var] = _create_variable!(Upper(bilevel_model), var)
    end

    # create lower level variables
    for a in instance.master.A
        lower_vars_ref[instance.subproblems[1].link_varsC[a]] = upper_vars_ref[instance.master.link_vars[a]]
    end

    for var in all_variables(instance.subproblems[1].mip_model)
        if !haskey(lower_vars_ref, var)
            lower_vars_ref[var] = _create_variable!(Lower(bilevel_model), var)
        end
    end

    # define upper level objective
    upper_obj = objective_function(instance.master.model)
    @objective(
        Upper(bilevel_model),
        Min,
        upper_obj.constant +  # TODO: The constant part of the objective is not supported by the current BilevelJuMP version. Thus, the solution values may differ by this constant.
        sum(coef * upper_vars_ref[var] for (var, coef) in upper_obj.terms) +
        sum(
            sum(coef * lower_vars_ref[var] for (var, coef) in sub_problem.r_objterm.terms)
            for sub_problem in instance.subproblems
        )  # TODO: Currently only supports the SubSolverJuMP subproblem type
    )

    # define lower level objective
    lower_obj = objective_function(instance.subproblems[1].mip_model)
    @objective(
        Lower(bilevel_model),
        Min,
        sum(coef * lower_vars_ref[var] for (var, coef) in lower_obj.terms)
    )

    # add upper level constraints
    _add_constraints!(Upper(bilevel_model), instance.master.model, upper_vars_ref)

    # add lower level constraints
    _add_constraints!(Lower(bilevel_model), instance.subproblems[1].mip_model, lower_vars_ref)

    return Instance(MibSMaster(bilevel_model), nothing)
end
