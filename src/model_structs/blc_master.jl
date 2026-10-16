# The master problem in the Benders-like Cuts model, i.e., the high point relaxation and a rule on how to generate big Ms.

"""
Master wrapper for the Benders-like Cuts (`BlC`) solver.

`BlCMaster` stores the high-point relaxation, binary linking variables,
subproblem objective expressions in the relaxation, and a user-provided
big-M function for BlC cut generation.
"""
struct BlCMaster{T}
    hpr::JuMP.Model  # The high-point-relaxation of the bilevel problem
    A::Vector{T}  # Iterable of common resources
    link_vars::Dict{T,VariableRef}  # Linking variables as dict. Keys=A
    big_m::Function  # A function that takes an a in A and sub_problem name and returns the corresponding big M
    sub_names::Vector{String}  # List of sub_problem names
    sub_objectives::Any  # Dict: sub name => sub_problem objective function (in hpr variables)
end

"""Evaluate a BlC big-M callback using the extended or legacy signature."""
function _evaluate_blc_big_m(big_m::Function, a, subproblem, current_cost, x_values)
    value = if applicable(big_m, a, subproblem, current_cost, x_values)
        big_m(a, subproblem, current_cost, x_values)
    else
        big_m(a, subproblem)
    end
    return _validate_blc_big_m(value)
end

function _validate_blc_big_m(value)
    value isa Real && isfinite(value) && value >= 0 ||
        throw(DomainError(value, "BlC big-M coefficients must be finite and nonnegative."))
    return value
end

function _apply_linking_branching_rule!(model::JuMP.Model, link_vars, params::SolverParam)
    rule = get(params.stats.data, "branching_rule", "default")
    if rule == "fixed_linking_order"
        _set_fixed_linking_branch_priorities!(model, link_vars)
    elseif rule == "linking_first"
        _set_linking_branch_priorities!(model, link_vars)
    elseif rule != "default"
        throw(ArgumentError("Unknown branching_rule '$rule'. Use default, linking_first, or fixed_linking_order."))
    end
    return nothing
end

"""Prioritize master linking variables for Gurobi branching.

All integer variables receive the default priority 1; linking variables
receive priority 2 and are therefore considered first by Gurobi.
"""
function _set_linking_branch_priorities!(model::JuMP.Model, link_vars)
    integer_variables = [
        variable for variable in all_variables(model) if is_binary(variable) || is_integer(variable)
    ]
    link_variables = collect(values(link_vars))

    for variable in integer_variables
        set_attribute(variable, Gurobi.VariableAttribute("BranchPriority"), Int64(1))
    end
    for variable in link_variables
        set_attribute(variable, Gurobi.VariableAttribute("BranchPriority"), Int64(2))
    end
    return nothing
end

"""Assign distinct deterministic priorities within linking/non-linking classes."""
function _set_fixed_linking_branch_priorities!(model::JuMP.Model, link_vars)
    integer_variables = [
        variable for variable in all_variables(model) if is_binary(variable) || is_integer(variable)
    ]
    link_indices = Set(index(variable) for variable in values(link_vars))
    link_variables = [variable for variable in integer_variables if index(variable) in link_indices]
    other_variables = [variable for variable in integer_variables if !(index(variable) in link_indices)]
    sort_key(variable) = (string(JuMP.name(variable)), string(index(variable)))
    sort!(link_variables; by=sort_key)
    sort!(other_variables; by=sort_key)

    # Higher priorities branch first. Linking variables remain above all others.
    for (rank, variable) in enumerate(other_variables)
        set_attribute(variable, Gurobi.VariableAttribute("BranchPriority"), Int64(length(other_variables) - rank + 1))
    end
    offset = length(other_variables)
    for (rank, variable) in enumerate(link_variables)
        set_attribute(variable, Gurobi.VariableAttribute("BranchPriority"), Int64(offset + length(link_variables) - rank + 1))
    end
    return nothing
end

function check(master::BlCMaster, params::SolverParam)
    # We need one linking variable for each resource
    if !(length(master.A) == length(master.link_vars))
        error(
            "You need to provide the same number of master linking variables as there are resources. |A| = $(length(master.A)) but |vars| = $(length(master.link_vars)).",
        )
    end

    # Check that linking variables are in the master problem
    for a in master.A
        if !is_valid(master.hpr, master.link_vars[a])
            error(
                "In the master problem, there is no linking variable for resource $a in the master JuMP model.",
            )
        end
    end

    # We need a big M for each resource a in A and sub_problem
    for sub in master.sub_names
        for a in master.A
            # Cost-aware callbacks need an actual follower solution. Validate
            # their values when the cut is generated, not using `nothing`.
            if applicable(master.big_m, a, sub)
                _validate_blc_big_m(master.big_m(a, sub))
            end
        end
    end
end
