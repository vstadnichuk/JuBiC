# The master problem in the Benders-like Cuts model, i.e., the high point relaxation and a rule on how to generate big Ms.

struct BlCMaster{T}
    hpr::JuMP.Model  # The high-point-relaxation of the bilevel problem
    A::Vector{T}  # Iterable of common resources
    link_vars::Dict{T,VariableRef}  # Linking variables as dict. Keys=A
    big_m::Function  # A function that takes an a in A and sub_problem name and returns the corresponding big M
    sub_names::Vector{String}  # List of sub_problem names
    sub_objectives::Any  # Dict: sub name => sub_problem objective function (in hpr variables)
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
            big_m = master.big_m(a, sub)
            if big_m < 0
                error(
                    "For resource $a and sub_problem $sub, a negative big M $big_m was returned.",
                )
            end
        end
    end
end
