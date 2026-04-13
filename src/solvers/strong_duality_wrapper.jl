using JuMP
using MathOptInterface

const MOI = MathOptInterface

"""
    build_strong_duality_mip_instance(master_model, x_link, subproblems, big_m)

Experimental wrapper that builds a monolithic MIP reformulation from a
first-level JuMP model and a list of `SubSolverJuMP` LP followers.

The current trial implementation assumes:
- the master objective sense is minimization
- all followers are LPs once the internally generated `link_varsC` variables are
  ignored
- the only leader-follower coupling in each subproblem is the standard JuBiC
  linking pattern `y[a] <= x[a]`
- `x_link` is a dictionary-like object keyed compatibly with `sub.A`

The returned instance should be solved with JuBiC's existing `MIPparam`.
"""
function build_strong_duality_mip_instance(
    master_model::JuMP.Model,
    x_link,
    subproblems,
    big_m::Function,
)
    @info "The generic strong-duality reformulation is currently an experimental feature."
    @info "Model generation may take some time because JuBiC rebuilds fresh primal and dual blocks for all follower models."

    x_link isa AbstractDict ||
        throw(ArgumentError("The current experimental strong-duality wrapper expects the linking variables `x_link` as a dictionary-like object keyed by the follower resources."))
    all(sub -> sub isa SubSolverJuMP, subproblems) ||
        throw(ArgumentError("The current experimental strong-duality wrapper only supports SubSolverJuMP followers."))
    objective_sense(master_model) == MOI.MIN_SENSE ||
        throw(ArgumentError("The current experimental strong-duality wrapper only supports minimization master models."))

    mip, ref_map = copy_model(master_model)
    x_master = Dict(key => ref_map[var] for (key, var) in pairs(x_link))

    leader_terms = Any[]
    for sub in subproblems
        push!(leader_terms, _add_sd_block!(mip, x_master, sub, big_m))
    end

    current_obj = objective_function(mip)
    @objective(mip, Min, current_obj + sum(leader_terms))
    return Instance(MIPMaster(mip), nothing)
end

function _add_sd_block!(
    mip::JuMP.Model,
    x_master::AbstractDict,
    sub::SubSolverJuMP,
    big_m::Function,
)
    follower = sub.mip_model
    objective_sense(follower) == MOI.MIN_SENSE ||
        throw(ArgumentError("Subproblem $(sub.name) must be a minimization LP to use the experimental strong-duality wrapper."))

    link_keys = Set(keys(sub.link_varsC))
    link_vars = Set(values(sub.link_varsC))
    regular_vars = VariableRef[]
    for var in all_variables(follower)
        if !(var in link_vars)
            if is_binary(var) || is_integer(var)
                throw(ArgumentError("Subproblem $(sub.name) contains integer follower variables. The experimental strong-duality wrapper currently supports LP followers only."))
            end
            push!(regular_vars, var)
        end
    end

    primal_var_map = Dict{VariableRef,VariableRef}()
    for (idx, var) in enumerate(regular_vars)
        basename = JuMP.name(var)
        if isempty(basename)
            basename = "var$(idx)"
        end
        primal_var_map[var] = @variable(mip, base_name = "sd_$(sub.name)_$(basename)")
    end

    rows = Any[]
    for (F, S) in list_of_constraint_types(follower)
        if F == VariableRef
            for cref in all_constraints(follower, F, S)
                obj = constraint_object(cref)
                var = obj.func
                if var in link_vars
                    continue
                end
                push!(rows, _copy_variable_constraint!(mip, primal_var_map, sub, cref, obj))
            end
        elseif F == AffExpr
            for cref in all_constraints(follower, F, S)
                obj = constraint_object(cref)
                if _affexpr_involves_link_vars(obj.func, link_vars)
                    _validate_standard_link_constraint(obj, sub)
                    continue
                end
                push!(rows, _copy_affine_constraint!(mip, primal_var_map, obj))
            end
        else
            throw(ArgumentError("Unsupported follower constraint type $(F) in subproblem $(sub.name)."))
        end
    end

    for a in sub.A
        haskey(x_master, a) || throw(ArgumentError("The master linking variable dictionary does not contain key $(a), which is required by subproblem $(sub.name)."))
        y_var = primal_var_map[sub.y_vars[a]]
        @constraint(mip, y_var <= x_master[a], base_name = "sd_link_$(sub.name)_$(a)")
        push!(
            rows,
            (
                sense = :le,
                coeffs = Dict(y_var => 1.0),
                rhs_constant = 0.0,
                rhs_x_terms = [(a, 1.0)],
            ),
        )
    end

    copied_c_obj = _copy_affexpr_to_model(sub.c_objterm, primal_var_map)
    copied_r_obj = _copy_affexpr_to_model(sub.r_objterm, primal_var_map)

    dual_vars = Dict{Int,VariableRef}()
    for (row_idx, row) in enumerate(rows)
        dual_vars[row_idx] = _new_dual_variable(mip, sub.name, row_idx, row.sense, row.rhs_x_terms)
    end

    c_coeffs = _affexpr_coefficients(sub.c_objterm)
    c_constant = JuMP.constant(sub.c_objterm)
    !iszero(c_constant) &&
        @warn "The follower objective of subproblem $(sub.name) contains a constant term. The experimental SD wrapper keeps it in strong duality, but this case has not been stress-tested yet."

    for primal_var in values(primal_var_map)
        coeff_expr = sum(get(row.coeffs, primal_var, 0.0) * dual_vars[row_idx] for (row_idx, row) in enumerate(rows))
        c_coeff = get(c_coeffs, _invert_lookup(primal_var_map, primal_var), 0.0)
        @constraint(mip, coeff_expr == c_coeff, base_name = "sd_dual_feas_$(sub.name)_$(JuMP.name(primal_var))")
    end

    dual_obj = AffExpr(0.0)
    for (row_idx, row) in enumerate(rows)
        π = dual_vars[row_idx]
        add_to_expression!(dual_obj, row.rhs_constant, π)
        for (akey, coeff) in row.rhs_x_terms
            link_bound = float(big_m(akey, sub.name))
            product_var = _linearize_binary_product!(mip, x_master[akey], π, row.sense, link_bound, sub.name, row_idx, akey)
            add_to_expression!(dual_obj, coeff, product_var)
        end
    end

    @constraint(mip, copied_c_obj == dual_obj, base_name = "sd_strong_duality_$(sub.name)")
    return copied_r_obj
end

function _copy_variable_constraint!(mip, primal_var_map, sub, cref, obj)
    set = obj.set
    new_var = primal_var_map[obj.func]
    if set isa MOI.GreaterThan{Float64}
        @constraint(mip, new_var >= set.lower, base_name = "sd_bound_ge_$(sub.name)_$(JuMP.name(new_var))")
        return (sense = :ge, coeffs = Dict(new_var => 1.0), rhs_constant = set.lower, rhs_x_terms = Tuple{Any,Float64}[])
    elseif set isa MOI.LessThan{Float64}
        @constraint(mip, new_var <= set.upper, base_name = "sd_bound_le_$(sub.name)_$(JuMP.name(new_var))")
        return (sense = :le, coeffs = Dict(new_var => 1.0), rhs_constant = set.upper, rhs_x_terms = Tuple{Any,Float64}[])
    elseif set isa MOI.EqualTo{Float64}
        @constraint(mip, new_var == set.value, base_name = "sd_bound_eq_$(sub.name)_$(JuMP.name(new_var))")
        return (sense = :eq, coeffs = Dict(new_var => 1.0), rhs_constant = set.value, rhs_x_terms = Tuple{Any,Float64}[])
    elseif set isa MOI.ZeroOne || set isa MOI.Integer
        throw(ArgumentError("Subproblem $(sub.name) contains integer follower variables. The experimental strong-duality wrapper currently supports LP followers only."))
    else
        throw(ArgumentError("Unsupported variable constraint set $(typeof(set)) in subproblem $(sub.name)."))
    end
end

function _copy_affine_constraint!(mip, primal_var_map, obj)
    coeffs = Dict{VariableRef,Float64}()
    new_expr = AffExpr(JuMP.constant(obj.func))
    for (coef, var) in JuMP.linear_terms(obj.func)
        new_var = primal_var_map[var]
        add_to_expression!(new_expr, coef, new_var)
        coeffs[new_var] = get(coeffs, new_var, 0.0) + coef
    end

    set = obj.set
    constant = JuMP.constant(obj.func)
    if set isa MOI.LessThan{Float64}
        @constraint(mip, new_expr <= set.upper)
        return (sense = :le, coeffs = coeffs, rhs_constant = set.upper - constant, rhs_x_terms = Tuple{Any,Float64}[])
    elseif set isa MOI.GreaterThan{Float64}
        @constraint(mip, new_expr >= set.lower)
        return (sense = :ge, coeffs = coeffs, rhs_constant = set.lower - constant, rhs_x_terms = Tuple{Any,Float64}[])
    elseif set isa MOI.EqualTo{Float64}
        @constraint(mip, new_expr == set.value)
        return (sense = :eq, coeffs = coeffs, rhs_constant = set.value - constant, rhs_x_terms = Tuple{Any,Float64}[])
    end
    throw(ArgumentError("Unsupported affine constraint set $(typeof(set)) in experimental strong-duality wrapper."))
end

function _affexpr_involves_link_vars(expr::AffExpr, link_vars::Set{VariableRef})
    return any(var in link_vars for (_, var) in JuMP.linear_terms(expr))
end

function _validate_standard_link_constraint(obj, sub::SubSolverJuMP)
    obj.set isa MOI.LessThan{Float64} || throw(ArgumentError("Subproblem $(sub.name) contains a constraint involving internal linking copies that is not of the supported `y <= x` type."))
    iszero(obj.set.upper) || throw(ArgumentError("Subproblem $(sub.name) contains a linking constraint with nonzero right-hand side."))
    coeffs = Dict(var => coef for (coef, var) in JuMP.linear_terms(obj.func))
    for a in sub.A
        y = sub.y_vars[a]
        x = sub.link_varsC[a]
        if length(coeffs) == 2 && get(coeffs, y, 0.0) == 1.0 && get(coeffs, x, 0.0) == -1.0
            return true
        end
    end
    throw(ArgumentError("Subproblem $(sub.name) contains a constraint with internal linking copies that JuBiC cannot interpret as the standard `y <= x` link."))
end

function _copy_affexpr_to_model(expr::AffExpr, primal_var_map::Dict{VariableRef,VariableRef})
    copied = AffExpr(JuMP.constant(expr))
    for (coef, var) in JuMP.linear_terms(expr)
        haskey(primal_var_map, var) || throw(ArgumentError("The experimental strong-duality wrapper found an objective term depending on an unsupported variable $(JuMP.name(var))."))
        add_to_expression!(copied, coef, primal_var_map[var])
    end
    return copied
end

function _affexpr_coefficients(expr::AffExpr)
    coeffs = Dict{VariableRef,Float64}()
    for (coef, var) in JuMP.linear_terms(expr)
        coeffs[var] = get(coeffs, var, 0.0) + coef
    end
    return coeffs
end

function _invert_lookup(var_map::Dict{VariableRef,VariableRef}, copied_var::VariableRef)
    for (old_var, new_var) in var_map
        new_var == copied_var && return old_var
    end
    error("Could not invert the variable map in the experimental strong-duality wrapper.")
end

function _new_dual_variable(mip, sub_name, row_idx, sense::Symbol, rhs_x_terms)
    if sense == :le
        if isempty(rhs_x_terms)
            return @variable(mip, upper_bound = 0.0, base_name = "sd_pi_$(sub_name)_$(row_idx)")
        end
        return @variable(mip, upper_bound = 0.0, base_name = "sd_pi_$(sub_name)_$(row_idx)")
    elseif sense == :ge
        return @variable(mip, lower_bound = 0.0, base_name = "sd_pi_$(sub_name)_$(row_idx)")
    elseif sense == :eq
        return @variable(mip, base_name = "sd_pi_$(sub_name)_$(row_idx)")
    end
    error("Unknown row sense $(sense) in experimental strong-duality wrapper.")
end

function _linearize_binary_product!(mip, x::VariableRef, π::VariableRef, sense::Symbol, big_m::Float64, sub_name, row_idx, akey)
    big_m > 0 || throw(ArgumentError("The experimental strong-duality wrapper requires strictly positive big-M values."))
    lower, upper = if sense == :le
        (-big_m, 0.0)
    elseif sense == :ge
        (0.0, big_m)
    else
        (-big_m, big_m)
    end
    set_lower_bound(π, max(has_lower_bound(π) ? lower_bound(π) : -Inf, lower))
    set_upper_bound(π, min(has_upper_bound(π) ? upper_bound(π) : Inf, upper))

    z = @variable(mip, lower_bound = lower, upper_bound = upper, base_name = "sd_prod_$(sub_name)_$(row_idx)_$(hash(akey))")
    @constraint(mip, z >= lower * x)
    @constraint(mip, z <= upper * x)
    @constraint(mip, z >= π - upper * (1 - x))
    @constraint(mip, z <= π - lower * (1 - x))
    return z
end
