using JuMP, BilevelJuMP

"""
    transform_GBC_to_MibS(instance::Instance, solver::SolverWrapper) -> Instance

Transforms a GBC instance into a MiBS instance.
"""
function transform_GBC_to_MibS(instance::Instance, solver::SolverWrapper)
    #@assert length(instance.subproblems) == 1 "Model is not a bilevel problem (there are $(length(instance.subproblems)) != 1 subproblems)."
    @assert instance.master isa Master "Model is not a bilevel problem (master is not of type Master)."
    @assert instance.subproblems[1] isa SubSolverJuMP "The transformation is currently only implemented for SubSolverJuMP subproblems."  # TODO: Add support for different subproblem types

    if length(instance.subproblems) > 1
        @info "As MiBS can only deal with MIB-MIB instances with only one subproblem. But, the provided instance has $(length(instance.subproblems)) subproblems.
        Therefore, we now merge the subproblems into an instance with a single subproblem."
        instance = merge_subproblems(instance, solver)
    end

    bilevel_model = BilevelJuMP.BilevelModel()
    upper_vars_ref = Dict{Any,Any}()
    lower_vars_ref = Dict{Any,Any}()

    

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
    l1term = if upper_obj isa GenericVariableRef upper_obj else sum(coef * upper_vars_ref[var] for (var, coef) in upper_obj.terms) + constant(upper_obj) end # small hack if L1 Obj is only 1 variable
    @objective(
        Upper(bilevel_model),
        Min,
        sum(l1term) +
        sum(
            constant(sub_problem.r_objterm) +
            sum(coef * lower_vars_ref[var] for (var, coef) in sub_problem.r_objterm.terms)
            for sub_problem in instance.subproblems
        )  # TODO: Currently only supports the SubSolverJuMP subproblem type
    )
    # TODO: "fun fact" MiBS ignores constants in objective function O_O

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



"""
    function merge_subproblems(instance::Instance, solver::SolverWrapper) -> Instance

Create a new Instance where there is only one subproblem that is the merged version of all current subproblems. 

"""
function merge_subproblems(instance::Instance, solver::SolverWrapper)
    @assert instance.master isa Master "Model is not a bilevel problem (master is not of type Master)."
    for instance in instance.subproblems
        @assert instance isa SubSolverJuMP "The transformation is currently only implemented for SubSolverJuMP subproblems."
    end

    mainsub = JuMP.Model(() -> get_next_optimizer(solver))
    
    # create link_copy and y vars in mainsub
    @variable(mainsub, xlink[instance.master.A], Bin, base_name="x_mcopy")

    # merge subproblems
    newvars_ref = Dict()
    r_obj_mainsub_ref = Dict()
    c_obj_mainsub_ref = Dict()
    for sub in instance.subproblems
        # add vars to mainsub
        for var in all_variables(sub.mip_model)
            # check if linking vars copy
            key_link = _key_for_value(sub.link_varsC, var, instance.master.A)
            if !isnothing(key_link)
                newvars_ref[var] = xlink[key_link]
                continue
            end

            # check if y var
            """
            key_link = _key_for_value(sub.y_vars,var, instance.master.A)
            if !isnothing(key_link)
                newvars_ref[var] = y_vars_main[key_link]
                continue
            end
            """

            # just model variable
            newvars_ref[var] = _create_variable!(mainsub, var)
        end

        # add constraints to mainsub
        _add_constraints!(mainsub, sub.mip_model, newvars_ref)

        # generate r and c objective terms with new variables
        r_obj_mainsub_ref[name(sub)] = constant(sub.r_objterm) + sum(coef * newvars_ref[var] for (var, coef) in sub.r_objterm.terms)
        c_obj_mainsub_ref[name(sub)] = constant(sub.c_objterm) + sum(coef * newvars_ref[var] for (var, coef) in sub.c_objterm.terms)
    end

    # add objective
    @objective(
        mainsub,
        Min,
        sum(values(c_obj_mainsub_ref))
    )

    # build new structs
    mainsubname = "mainsub"
    nmaster = Master(instance.master.model, instance.master.A, instance.master.link_vars, [mainsubname], instance.master.partial_decomposition, instance.master.objL2) 
    nsub = SubSolverJuMP(mainsubname, mainsub, instance.master.A, xlink, xlink, # TODO: setting y=x here is quite a hack but should work for MiBS. Not tested for other solver
        Dict(a => 1 for a in instance.master.A), sum(values(r_obj_mainsub_ref)), sum(values(c_obj_mainsub_ref)), timelimit -> (false, 0))
    return Instance(nmaster, [nsub])
end

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


function _key_for_value(vars, val, key_list)
    for k in key_list
        if vars[k] == val
            return k    # found it, return the key
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
