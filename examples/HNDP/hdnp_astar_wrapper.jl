# Transfer HNDP instance into subproblem within GBCSolver that is solved with A*-Search, i.e., implement necessary wrapper functions for subproblem
using JuBiC

import JuBiC: neighbours_w
import JuBiC: heuristic_w
import JuBiC: cost_w
import JuBiC: isgoal_w
import JuBiC: hashfn_w
import JuBiC: start_state
import JuBiC: end_state
import JuBiC: used_resource
import JuBiC: risk
import JuBiC: cost_subs
import JuBiC: prepare
import JuBiC: reconstract_path
import JuBiC: dominate_w



########## Implementation AStarSolver for HNDPwC ##########
struct Label
    mynode::Any
    myarc::Union{Tuple,Nothing}  # the arc used to obtain this node (as tuple). Or nothing if this is start label
    myweight::Any  # the accumulated weight of the path up to this label
    predecessor::Union{Label,Nothing}  # the predecessor (or nothing for start label)
end

"""
The structure that is contains information about a single user and is passed to the corresponding subsolver.
"""
mutable struct HNDPwC_Structure
    const graph  # the underlying network graph 
    const origin  # origin node
    const destination  # destination node
    const maxweight::Number  # maximum weight a path may accumulate

    const check_cycles::Bool  # if true, will check for cycles during the labeling algorithm (only needed if negative cost cycles may exist)

    const mrisk  # risk matrix
    const mcost  # cost matrix
    const mweight  # weight matrix

    const shortest_cost  # matrix with shortest path (cost) between nodes 
    const shortest_risk  # matrix with shortest path (risk) between nodes 
    const shortest_weight  # matrix with shortest path (weight) between nodes 
    shortest_adaptive::Any  # initialized as matrix with shortest path for each node pair acoridng to the current solution
end

"""
    isequal(t1::Label, t2::Label)

Two 'Label' are equal iff they have the same node. 
This prevents cycles while allowing us to store additional information within labels.
"""
function Base.:(==)(t1::Label, t2::Label)
    # TODO: double check if this makes sense
    return t1.mynode == t2.mynode
end

function Base.hash(t::Label)
    return Base.hash(t.mynode)
end

function neighbours_w(state::Label, xsol, structure::HNDPwC_Structure, params::SolverParam)
    nlist = []
    for neig in neighbors(structure.graph, state.mynode)
        if edge_allowed(xsol, state.mynode, neig)
            myarc = (state.mynode, neig)
            arcweight = structure.mweight[state.mynode, neig]
            newLabel = Label(neig, myarc, state.myweight + arcweight, state)
            if label_feasible(newLabel, structure)
                push!(nlist, newLabel)
            end
        end
    end
    return nlist
end

function heuristic_w(state::Label, goal::Label, xsol, structure::HNDPwC_Structure, cs::CostStructure, params::SolverParam)
    return structure.shortest_adaptive[state.mynode, goal.mynode]
end

function cost_w(current::Label, neighbour::Label, structure, cs::CostStructure, params::SolverParam)
    if cs.cost_state == MASTER_LEVEL
        return structure.mrisk[current.mynode, neighbour.mynode]
    elseif cs.cost_state == SUB_PROBLEM_LEVEL
        return structure.mcost[current.mynode, neighbour.mynode]
    else
        @assert cs.cost_state == CONNECTOR_BASED
        rc =
            structure.mcost[current.mynode, neighbour.mynode] * cs.gval +
            structure.mrisk[current.mynode, neighbour.mynode]
        rV = rc + cs.kvals[(current.mynode, neighbour.mynode)]
        return rV
    end
end

function isgoal_w(state::Label, goal, params::SolverParam)
    return state == goal
end

function start_state(sol::AStarSolver)
    return Label(sol.structure.origin, nothing, 0, nothing)
end

function end_state(sol::AStarSolver)
    return Label(sol.structure.destination, nothing, -1, nothing)  # we do not specify which arc we came from, what weight we should have, or which predecessor
end

function used_resource(state::Label, A, structure)
    if state.myarc in A
        return [state.myarc]
    else
        return []
    end
end

function risk(state::Label, neighbour::Label, structure::HNDPwC_Structure)
    return structure.mrisk[state.mynode, neighbour.mynode]
end

function cost_subs(state, neighbour, structure)
    return structure.mcost[state.mynode, neighbour.mynode]
end

function prepare(structure::HNDPwC_Structure, cs::CostStructure, params::SolverParam)
    # precompute the best distance matrix 
    if cs.cost_state == MASTER_LEVEL
        structure.shortest_adaptive = structure.shortest_risk
    elseif cs.cost_state == SUB_PROBLEM_LEVEL
        structure.shortest_adaptive = structure.shortest_cost
    else
        @assert cs.cost_state == CONNECTOR_BASED
        structure.shortest_adaptive = calculate_shortest_matrix(structure, cs)
    end
end

function reconstract_path(goal::Label)
    rV = [goal]
    if isnothing(goal.predecessor)
        return rV
    else
        prepend!(rV, reconstract_path(goal.predecessor))  # preserve order
        return rV
    end
end

function dominate_w(::Label, a, b, structure, cs::CostStructure, params::SolverParam)
    # label has to have same arc if it should dominate
    if a.me.mynode != b.me.mynode
        return false
    end

    # label dominates if cost and distance are lower
    if a.me.myweight >= b.me.myweight && a.cost >= b.cost
        return true
    end
    return false
end




########## Auxiliary functions ##########
function calculate_shortest_matrix(structure::HNDPwC_Structure, cs::CostStructure)
    @assert cs.cost_state == CONNECTOR_BASED
    mcscale = cs.gval .* structure.mcost
    madap = mcscale + structure.mrisk
    for a in keys(cs.kvals)
        madap[a[1], a[2]] += cs.kvals[a]
    end

    asp = floyd_warshall_shortest_paths(structure.graph, madap)
    return asp.dists
end

function edge_allowed(xsol, e1, e2)
    xvalue = xsol[e1, e2]
    return xvalue > 0.5  # xvalue in {0;1}
end

"""
    has_node(l::Label, n)

Check if 'l' or one of his predecessors contains the node 'n'.
"""
function has_node(l::Label, n)
    if l.mynode == n
        return true
    elseif isnothing(l.predecessor)
        return false
    else
        return has_node(l.predecessor, n)
    end
end


"""
    label_feasible(l::Label, structure::HNDPwC_Structure)

We check if the passed label is feasible. Currently, only check if goal is still reachable with remaining weight (l.myweight + min_remainin_gweight <= maxweight).
Return 'true' if all checks passed.
"""
function label_feasible(l::Label, structure::HNDPwC_Structure)
    # check if weight restrictions can still be satisfied
    min_remaining_weight = structure.shortest_weight[l.mynode, structure.destination]
    if l.myweight + min_remaining_weight > structure.maxweight
        return false
    end

    # check cycles 
    if structure.check_cycles
        if !isnothing(l.predecessor) && has_node(l.predecessor, l.mynode)
            # found cycle as mynode is contained in predecessors 
            return false
        end
    end

    return true
end
