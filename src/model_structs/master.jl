using JuMP

struct Master{T}
    model::JuMP.Model  # The master problem without second level
    A::Vector{T}  # Iterable of common resources
    link_vars::Dict{T,VariableRef}  # Linking variables as dict. Keys=A
    sub_names::Vector{String}  # List of sub_problem names

    # You can pass a function that takes 'm' and a list of master objective variables (Dict: a => var for a in A) and add some additional constraints/variables to 'm'.
    # This models a partial decomposition approach for the subproblems.
    partial_decomposition::Union{Function,Nothing} # TODO: We need a better documentation why this is needed. It is because you do not have the variables presenting the second-level objectives when initially building the model. This may not be oblivious to everyone
end

Master(model, A, link_vars, sub_names) = Master(model, A, link_vars, sub_names, nothing)


"""
    check(master::Master, params::SolverParam)

Check the master problem for consistency and correctness.
"""
function check(master::Master, params::SolverParam)
    # Check objective sense is minimization
    if !(objective_sense(master.model) == MIN_SENSE)
        error("The master must be formulated as a minimization problem.")
    end

    # We need one linking variable for each resource
    if !(length(master.A) == length(master.link_vars))
        error(
            "You need to provide the same number of master linking variables as there are resources. |A| = $(length(master.A)) but |vars| = $(length(master.link_vars)).",
        )
    end

    # Check that linking variables are in the master problem
    for a in master.A
        if !is_valid(master.model, master.link_vars[a])
            error(
                "In the master problem, there is no linking variable for resource $(a) in the master JuMP model.",
            )
        end
    end
end
