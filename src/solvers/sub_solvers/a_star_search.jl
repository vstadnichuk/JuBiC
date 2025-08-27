# This is a Subproblem solver implementation designed for subproblems that should be solved with a dynamic programm. 
#  It offers a highly flexible framework that allows you to implement 

struct AStarSolver{T} <: SubSolver
    # general solver settings
    name::String  # the name (unique identifier) of this sub_problem
    A::Vector{T}  # The set of resources
    structure::Any  # the structure of the problem. It will be passed to all wrapper functions
    link_constraints_capacities::Dict{T,<:Number}  # the capacity parameters in the interdiction linking constraints

    # astar specific settings
    max_cost::Number  # the maximal cost a path can have (before it is cutoff)
    check_cycles::Bool  # Sets the 'enable_closedset' parameter within the astar function of AStarSearch.jl. If you check for cycles in neighbour function, you can set it to false. It is recommended to set false here to avoid bugs
end

@enum AStarCostState begin
    MASTER_LEVEL  # Use the master problem (first level) objective
    SUB_PROBLEM_LEVEL  # Use the sub_problem (second level) objective
    CONNECTOR_BASED  # Use the objective induced by by passed ConnectorLP solution
end

"""
Information that is required to calculate the cost between states correctly. 
"""
struct CostStructure
    cost_state::AStarCostState
    sval::Any  # if non-applicable cost state -1 by default
    gval::Any  # if non-applicable cost state -1 by default
    kvals::Any  # empty dict if non-applicable cost state
end

CostStructure(cost_state::AStarCostState) =
    (cost_state == CONNECTOR_BASED) ?
    error("Cost structure lacking ConnectorLP parameters") :
    CostStructure(cost_state, -1, -1, Dict())


########## Wrapper functions for A*-search algorithm ##########

"""
    neighbours_w(state, A, xsol, structure)

This function should compute the set of neighbouring states. 
- 'state': The current state we search neighbours for. 
- 'xsol': Mapping resource => 0 (if forbidden) or 1 (if allowed).
- 'structure': Structure information stored in _AStarSolver_ object.
- 'params': the solver parameter.
"""
function neighbours_w(state, xsol, structure, params::SolverParam)
    error(
        "You have to provide an implementation of this neighbours function for your subtype.",
    )
end

"""
    heuristic_w(state, goal, xsol, structure, cs::CostStructure, params::SolverParam)

By reimplementing this function, you can provide a heuristic that estimates the remaining distance of the path. 
Note that you should provide an underestimation, i.e., an optimistic guess, of the remaining distance if you want the A*-search to produce correct results.
- 'state': The current state we search neighbours for. 
- 'goal' the goal state. 
- 'xsol': Mapping resource => 0 (if forbidden) or 1 (if allowed).
- 'structure': Structure information stored in _AStarSolver_ object.
- 'cs': The current cost structure to use as objective function.
- 'prarams': The solver parameters
"""
function heuristic_w(state, goal, xsol, structure, cs::CostStructure, params::SolverParam)
    return 0
end

"""
    cost_w(current, neighbour, structure, cs::CostStructure, params::SolverParam)

Should return the cost of transitioning form current state to the neighbour.
- 'current': The current state.
- 'neighbour': The next state we transition to.
- 'structure': Structure information stored in _AStarSolver_ object.
- 'cs': The cost structure to use for computing the transition cost.
- 'params': The parameters of the solver.
"""
function cost_w(current, neighbour, structure, cs::CostStructure, params::SolverParam)
    error("You have to provide an implementation of the cost calculation for your subtype.")
end

"""
    isgoal_w(state, goal, params::SolverParam)

Return true if 'state' is an end state for the A*-search. 'goal' is the state you passed with the 'end_state' function. 
"""
function isgoal_w(state, goal, params::SolverParam)
    error(
        "You have to provide an implementation of the 'isgoal' function for your subtype. There is no default implementation supported to avoid mistakes due to forgetting to implement this function.",
    )
end

"""
    hashfn_w(state)

Reimplement if you need custom stable hash function.
"""
function hashfn_w(state)
    return hash(state)
end


"""
    start_state(sol::AStarSolver)

Return the starting state that the A*-search starts with.
"""
function start_state(sol::AStarSolver)
    error("You have to provide the start state for your subtype.")
end

"""
    end_state(sol::AStarSolver)

Return the end state. The A*-search should terminate when reached. When you implement your own 'isgoal' function, please still provide a state here that then is passed to the 'isgoal' function.
"""
function end_state(sol::AStarSolver)
    error("You have to provide the end state for your subtype.")
end

"""
    used_resource(state, A, structure)

Return the list of resources used by this 'state'. If no resource is used, return an empty list. A 'state' can use multiple resources. Even if only 1 resource is used, return a list (with 1 element).
"""
function used_resource(state, A, structure)
    error(
        "You have to implement this function that returns the resources used by a state for your subtype.",
    )
end

"""
    risk(state, neighbour, structure)

Return the master objective value of transitioning form current state to the neighbour. 
- 'current': The current state.
- 'neighbour': The next state we transition to.
- 'structure': Structure information stored in _AStarSolver_ object.
"""
function risk(state, neighbour, structure)
    error(
        "You have to implement this function that returns the master objective function evaluated for a transition for your subtype.",
    )
end

"""
    risk(state, neighbour, structure)

Return the original cost, i.e. sub_problem objective value, of transitioning form current state to the neighbour. 
- 'current': The current state.
- 'neighbour': The next state we transition to.
- 'structure': Structure information stored in _AStarSolver_ object.
"""
function cost_subs(state, neighbour, structure)
    error(
        "You have to implement this function that returns the original cost, i.e. sub_problem objective value, for your subtype.",
    )
end

"""
    prepare(structure, cs::CostStructure, params::SolverParam)

This function is called before the A*-search is started to give you the oportunity to do some preprocessing steps (if needed).
- 'structure': Structure information stored in _AStarSolver_ object.
- 'cs': the current cost structure of the objective function.
- 'params': Solver parameter.
"""
function prepare(structure, cs::CostStructure, params::SolverParam)
    error(
        "Please provide a default implementation of the 'prepare' function. You can leave it empty, but a default implementation is not supported to avoi mistakes.",
    )
end

"""
    reconstract_path(goal)

Reconstructs the path from the passed label that reached the goal.
"""
function reconstract_path(goal)
    error(
        "You need to provide an implementation of this function 'reconstract_path' that computes the path based on the goal state.",
    )
end

"""
    dominate_w(me, a, b, structure, cs::CostStructure, params::SolverParam)

Check the dominance relationship between labels a and b.

# Arguments
- 'me': Passes an object of the type of the labels used within algorithms to specify which function to use. You can ignore the value but should specify the type here s.t. your dominance function is used. 
- 'a': A 'MyLabel' object from LabelingAlg.jl that contains not only the current label but also its cost and heuristic cost estimation. 
- 'b': A 'MyLabel' object from LabelingAlg.jl that contains not only the current label but also its cost and heuristic cost estimation. 
- 'structure': Structure information stored in _AStarSolver_ object.
- 'cs::CostStructure': the current cost structure of the objective function.
- 'params::SolverParam': Solver parameter.

# Returns
- 'dominated::Bool': True if a is dominated by b, i.e., b is strictly better that a.
"""
function dominate_w(me, a, b, structure, cs::CostStructure, params::SolverParam)
    error(
        "You have to provide an implementation of the dominance function for your subtype.",
    )
end

########## SubSolver functions ##########

function capacity_linking(sol::AStarSolver, a, params::SolverParam)
    return sol.link_constraints_capacities[a]
end

function check(sol::AStarSolver, params::SolverParam)
    # Do a very basic check if the necessary subfunctions were implemented
    xsold = Dict(a => 1 for a in sol.A)

    # create start and end label and test enighbours
    s = start_state(sol)
    e = end_state(sol)
    if !isgoal_w(e, e, params)
        error(
            "The #isgoal_w' function is not implemented correctly. end=end did result in false!",
        )
    end
    if hashfn_w(e) != hashfn_w(e)
        error(
            "The 'hashfn_w' function is implemented wrong. Two calls of hash produce different results!",
        )
    end
    nei = neighbours_w(s, xsold, sol.structure, params)
    if length(nei) == 0
        error(
            "Even if all resources are available, the starting state does not have any neighbours. Please check your implementation of the 'neighbours_w' function.",
        )
    end

    # it is not recommended to depend on the internal cycle avoiding strategie
    if sol.check_cycles
        printstyled(
            "It is not recommended to rely on the cycle avoiding algorithm in the A*-search implementation within sub_problem $(name(sol)). \n ";
            color=:cyan,
        )
    end

    # simple check that dominance wrapper returns a booleand and a label dominates itself (what should hold for any reasonable dominance function)
    if !dominate_w(
        s,
        MyLabel(s, 0, 0),
        MyLabel(s, 0, 0),
        sol.structure,
        CostStructure(MASTER_LEVEL),
        params,
    )
        error(
            "A label should dominate itself what it does not with the provided implementation of 'dominate_w' function",
        )
    end

    # test cost and resource functions
    c = cost_w(s, nei[1], sol.structure, CostStructure(SUB_PROBLEM_LEVEL), params)
    r = risk(s, nei[1], sol.structure)
    res = used_resource(s, sol.A, sol.structure)

    # basic tests passed
end

function compute_lower_bound_master_contribution(sol::AStarSolver, params::SolverParam, time_limit)
    cs = CostStructure(MASTER_LEVEL)
    xmapping = Dict(a => 1 for a in sol.A)

    result = run_astar(sol, xmapping, cs, params, time_limit)

    # check if solution was found
    if result.status == LS_NO_SOLUTION
        error(
            "We found no path for sub_problem $(sol.name) despite all resources being available!",
        )
    elseif result.status == LS_TIMEOUT
        @debug "We found no path for sub_problem $(sol.name) due to timeout!"
        throw(TimeoutException("We found no path for sub_problem $(sol.name) due to timeout!"))
    elseif result.status != LS_OPTIMAL
        error(
            "We have a status within returned 'LabelingSolution' that is neither of 'LS_OPTIMAL', 'LS_TIMEOUT', nor 'LS_NO_SOLUTION'.",
        )
    end
    @assert !isnothing(result.goal)

    # construct result path and return solution
    resultpath = reconstract_path(result.goal)
    rprisk = path_risk(resultpath, sol)
    @debug "Using L1 objective function with all resources available, we found for sub_problem $(sol.name) the solution value $(rprisk) with path $(resultpath)."
    return rprisk
end

function name(ss::AStarSolver)
    return ss.name
end

function separation!(sol::AStarSolver, sval, gvals, kvals::Dict, param::SolverParam, time_limit)
    # solve with A*-search
    cs = CostStructure(CONNECTOR_BASED, sval, gvals, kvals)
    xvalsd = Dict(a => 1 for a in sol.A)  # in separation, we have all resources but they may cost money
    result = run_astar(sol, xvalsd, cs, param, time_limit)

    # check if solution was found
    if result.status == LS_NO_SOLUTION
        error(
            "No path found for sub_problem $(sol.name) in separation using the new objective function based on g=$(gvals) and k=$(kvals).",
        )
    elseif result.status == LS_TIMEOUT
        @debug "No path found in 'separation' due to 'timeout' for sub_problem $(sol.name) using L2 objective function with x=$(xsol)."
        throw(TimeoutException("No path found in 'separation' due to 'timeout' for sub_problem $(sol.name) using L2 objective function with x=$(xsol)."))
    elseif result.status != LS_OPTIMAL
        error(
            "We have a status within returned 'LabelingSolution' that is neither of 'LS_OPTIMAL', 'LS_TIMEOUT', nor 'LS_NO_SOLUTION'.",
        )
    end

    # solution found. Construct path and return
    rpath = reconstract_path(result.goal)
    if length(rpath) == 0
        error(
            "Labeling algorithm returned that no path found but path computed by 'reconstract_path' function was empty for sub_problem $(sol.name) using L2 objective function with x=$(xvals).",
        )
    else
        rpvalue = path_cost_current(rpath, sol, cs, param)
        @debug "For sub_problem $(sol.name), we found the solution value $(rpvalue) (based on new merged cost function) with path $(rpath)."
        rip = resources_in_path(rpath, sol.A, sol.structure)
        is_violated = Bool(!(sval < rpvalue + 1e-6))
        risk = path_risk(rpath, sol)
        cost = path_cost_subs(rpath, sol)
        @debug "We found a violated=$(is_violated) constraint in sub_problem $(sol.name) A*-search because sval = $sval and rpvalue=$(rpvalue)"
        return SubSolution(is_violated, risk, cost, rip)
    end
end

function set_nthreads(sol::AStarSolver, n)
    printstyled(
        "Currently, multi-thread is not supported for A*-search algorithm within sub_problem $(name(sol)). \n ";
        color=:cyan,
    )
end

function solve_sub_for_x(sol::AStarSolver, xvals, params::SolverParam, time_limit)
    # solve with A*-search
    cs = CostStructure(SUB_PROBLEM_LEVEL)
    result = run_astar(sol, xvals, cs, params, time_limit)

    # check if solution was found
    if result.status == LS_NO_SOLUTION
        @debug "No path found for sub_problem $(sol.name) using L2 objective function with x=$(xvals)."
        return false, 0, Dict()
    elseif result.status == LS_TIMEOUT
        # we treat time out the same as no feasible path found
        @debug "No path found due to 'timeout' for sub_problem $(sol.name) using L2 objective function with x=$(xvals)."
        throw(TimeoutException("No path found due to 'timeout' for sub_problem $(sol.name) using L2 objective function and first-level solution x=$(xvals)."))
    elseif result.status != LS_OPTIMAL
        error(
            "We have a status within returned 'LabelingSolution' that is neither of 'LS_OPTIMAL', 'LS_TIMEOUT', nor 'LS_NO_SOLUTION'.",
        )
    end

    # solution found. Construct path and return
    rpath = reconstract_path(result.goal)
    if length(rpath) == 0
        error(
            "Labeling algorithm returned that no path found but path computed by 'reconstract_path' function was empty for sub_problem $(sol.name) using L2 objective function with x=$(xvals).",
        )
    else
        rpvalue = path_cost_subs(rpath, sol)
        @debug "Using L2 objective function with x=$(xvals), we found for sub_problem $(sol.name) the solution value $(rpvalue) with path $(rpath)."
        rip = resources_in_path(rpath, sol.A, sol.structure)
        rmap = Dict(a => (a in rip ? 1 : 0) for a in sol.A)
        return true, rpvalue, rmap
    end
end




########## Auxiliary functions ##########
"""
    path_risk(path, sol::AStarSolver)

Compute the path risk, i.e., its value according to master objective function.
"""
function path_risk(path, sol::AStarSolver)
    r = 0
    for e in zip(path[1:(end-1)], path[2:end])
        r += risk(e[1], e[2], sol.structure)
    end
    return r
end

"""
    path_cost_subs(path, sol::AStarSolver)

Compute path cost based on original subs objective function. 
"""
function path_cost_subs(path, sol::AStarSolver)
    c = 0
    for e in zip(path[1:(end-1)], path[2:end])
        c += cost_subs(e[1], e[2], sol.structure)
    end
    return c
end

"""
    path_cost_current(path, sol::AStarSolver, cs::CostStructure, params::SolverParam)

Compute path cost based on the current cost structure. 
"""
function path_cost_current(path, sol::AStarSolver, cs::CostStructure, params::SolverParam)
    c = 0
    for e in zip(path[1:(end-1)], path[2:end])
        cc = cost_w(e[1], e[2], sol.structure, cs, params)
        c += cc
    end
    return c
end

function resources_in_path(path, A, structure)
    resources = Set()
    for state in path
        union!(resources, used_resource(state, A, structure))
    end
    return collect(resources)
end


"""
    run_astar(sol::AStarSolver, xmapping, cs::CostStructure, params::SolverParam, time_limit)

Prepare and execute the labeling algorithm for the A*-search. 
**return:** The 'LabelingSolution' object that contains the found path (as custom label), the status of termination, and the runtime
"""
function run_astar(sol::AStarSolver, xmapping, cs::CostStructure, params::SolverParam, time_limit)
    prepare(sol.structure, cs, params)
    neighborhood(x) = neighbours_w(x, xmapping, sol.structure, params)
    heuristicf(x, g) = heuristic_w(x, g, xmapping, sol.structure, cs, params)
    costf(x, y) = cost_w(x, y, sol.structure, cs, params)
    isgoalf(x) = isgoal_w(x, end_state(sol), params)
    dominate(a, b) = dominate_w(start_state(sol), a, b, sol.structure, cs, params)
    result = shortest_path_labeling(
        neighborhood,
        start_state(sol),
        end_state(sol),
        heuristicf,
        costf,
        isgoalf,
        dominate,
        sol.max_cost,
        time_limit,
    )
    return result
end
