# This is an own implementation of a labeling algorithm/A*-search. Compared to the 'AStarSearch.jl', it includes the option to prune labels.
## The labeling algorithm is implemented in a generic way that intends the users to manage subconstraints with the ''neighbours'' function. The labeling algorithms then handles the whole propagation of labels.
## Most of the time, you want to use custom labels that are passed to ''neighbours'' and ''dominance'' function to handle these checks. 
using DataStructures

########## We return a solution struct with solver information ##########
@enum LabelState begin
    LS_OPTIMAL = 0  # Terminated by finding best solution
    LS_NO_SOLUTION = 1  # All labels traversed but no solution found, i.e., nothing as solution is returned
    LS_TIMEOUT = 2  # Timeout reached before solution found
end

"""
The return value of the labeling algorithm. 
Returns the state in which algorithm terminated. If a solution was found. Also returns the goal label. 
"""
struct LabelingSolution
    goal::Any  # the final goal label (in user custom data structure). Nothing if no solution found
    status::LabelState  # the state in which the solver terminated
    costs::Any  # the cost of the solution (if found)
    runtime::Number # the time required to execute the labeling algorithm (in seconds)
end



########## The main labeling algorithms ##########
"""
An auxilliary structure that saves additional information for fast access on top of user defined state structure.
"""
struct MyLabel
    me::Any  # the user label structure 
    cost::Any  # the cost estimation for this label
    heu::Any  # heuristic underestimation of remaining cost to goal 
end

"""
    Base.isless(a::MyLabel, b::MyLabel)

Define a comparator function for labels.
This function compares labels by the estimated cost (cost + heuristic cost)
"""
function Base.isless(a::MyLabel, b::MyLabel)
    return (a.cost + a.heu <= b.cost + b.heu)
end


"""
    shortest_path_labeling(neighbours, starts, ends, heuristicf, costf, isgoal, dominance, max_cost)

Performs the labeling algorithm steps. Propagetes the start label till no label is left or first label reaching the goal is found. 
- 'neighbours': Mapping of label to its neighbours.
- 'starts': Start label.
- 'ends': Exemplary last node label. We call heuristic function for it and current label. 
- 'heuristicf': Mapping taking two labels and providing an underestimation for their distance. 
- 'costf': Mapping of two labels to the cost of moving from first to second. 
- 'isgoal': Function returning true if passed label is an end label, i.e., goal is reached. 
- 'dominance': Mapping taking two labels and returning true if second dominates the first. Note that as labels we pass our 'MyLabel' object that also contain cost of path till that point and heuristic cost estimation.
- 'max_cost': An upper bound on the cost of a label. If a label exceeds it, it is discarded. 
- 'time_limit': The time_limit in seconds.

"""
function shortest_path_labeling(neighbours, starts, ends, heuristicf, costf, isgoal, dominance, max_cost, time_limit)
    # TODO: One could definitly implement a faster labeling algorithms by doing smarted dominance checks that just going over the completle list of labels in each iteration. 
    # capture the starting time
    start_time = time()

    # open_list: labels to be extended, stored in a SortedSet ordered by cost.
    initial_label = MyLabel(starts, 0, heuristicf(starts, ends))
    open_list = SortedSet{typeof(initial_label)}()
    push!(open_list, initial_label)

    while !isempty(open_list)
        # Check if the elapsed time exceeds the time_limit.
        if time() - start_time > time_limit
            return LabelingSolution(nothing, LS_TIMEOUT, 0, max(0, time() - start_time))
        end

        # Retrieve and remove the label with the smallest cost.
        current_label = pop!(open_list)

        # If the current label is a goal state, add it to the solution set.
        if isgoal(current_label.me)
            return LabelingSolution(current_label.me, LS_OPTIMAL, current_label.cost, max(0, time() - start_time))
        end

        # Expand the current label.
        for new_label_o in neighbours(current_label.me)
            # build MyLabel structure
            new_label = MyLabel(new_label_o, current_label.cost + costf(current_label.me, new_label_o), heuristicf(new_label_o, ends))
            dominated = false

            # check if label satisfies max_cost restriction. If not, discard
            if new_label.cost + new_label.heu > max_cost
                continue
            end

            # Check if new_label is dominated by any label in open_list.
            for label in open_list
                if dominance(new_label, label)
                    dominated = true
                    break
                end
            end

            # If new_label is dominated, skip to the next one.
            if dominated
                continue
            end

            # Remove any labels dominated by new_label from open_list.
            to_delete = []
            for label in open_list
                if dominance(label, new_label)
                    push!(to_delete, label)
                end
            end
            for label in to_delete
                delete!(open_list, label)
            end

            # Add the new label to the open list.
            push!(open_list, new_label)
        end
    end

    # No solution found
    return LabelingSolution(nothing, LS_NO_SOLUTION, 0, max(0, time() - start_time))
end



"""
    enumerate_all_paths(neighbours, starts, ends, heuristicf, costf, isgoal, dominance, max_cost)

Enumerates all feasible paths. 
- 'neighbours': Mapping of label to its neighbours.
- 'starts': Start label.
- 'ends': Exemplary last node label. We call heuristic function for it and current label. 
- 'costf': Mapping of two labels to the cost of moving from first to second. 
- 'isgoal': Function returning true if passed label is an end label, i.e., goal is reached. 
- 'dominance': Mapping taking two labels and returning true if second dominates the first. Note that the concept of dominated paths must be defined for your enumeration algorithm. Otherwise, provide a function that always returns false.
- 'max_cost': An upper bound on the cost of a label. If a label exceeds it, it is discarded. Again, note that paths that correspond to discarded labels will not appear in the output. 
- 'time_limit': The time_limit in seconds. If the time_limit is reached, the algorithms returns alls paths found till that moment.

**return**: A ''LabelingSolution'' object where the ''goal'' element is the list of found paths. Note that the status ''LS_NO_SOLUTION'' is returned if no timeout occured, even if the list of paths found is empty.

"""
function enumerate_all_paths(neighbours, starts, ends, costf, isgoal, dominance, max_cost, time_limit)
    # Note: i'm aware that the code has high similarities to ''shortest_path_labeling'' function but I believe here it is easier to mentain two versions of code
    # TODO: Use more clear notation what is a passed label from input and what is our ''MyLabel'' 

    # capture the starting time
    start_time = time()

    # collection of all found paths represented by final label, i.e., final label 
    goal_labels = []
    goal_costs = []

    # open_list: labels to be extended, stored in a SortedSet ordered by cost.
    initial_label = MyLabel(starts, 0, 0)
    open_list = SortedSet{typeof(initial_label)}()
    push!(open_list, initial_label)

    while !isempty(open_list)
        # Check if the elapsed time exceeds the time_limit.
        if time() - start_time > time_limit
            return LabelingSolution(goal_labels, LS_TIMEOUT, goal_costs, max(0, time() - start_time))
        end

        # Retrieve and remove the label with the smallest cost.
        current_label = pop!(open_list)

        # If the current label is a goal state, add it to the solution set.
        if isgoal(current_label.me)
            push!(goal_labels, current_label.me)
            push!(goal_costs, current_label.cost)
        end

        # Expand the current label.
        for new_label_o in neighbours(current_label.me)
            # build MyLabel structure
            new_label = MyLabel(new_label_o, current_label.cost + costf(current_label.me, new_label_o), 0)
            dominated = false

            # check if label satisfies max_cost restriction. If not, discard
            if new_label.cost + new_label.heu > max_cost
                continue
            end

            # Check if new_label is dominated by any label in open_list.
            for label in open_list
                if dominance(new_label, label)
                    dominated = true
                    break
                end
            end

            # If new_label is dominated, skip to the next one.
            if dominated
                continue
            end

            # Remove any labels dominated by new_label from open_list.
            to_delete = []
            for label in open_list
                if dominance(label, new_label)
                    push!(to_delete, label)
                end
            end
            for label in to_delete
                delete!(open_list, label)
            end

            # Add the new label to the open list.
            push!(open_list, new_label)
        end
    end

    # Enumeration finished, return solution
    if length(goal_labels) == 0
        return LabelingSolution(goal_labels, LS_NO_SOLUTION, goal_costs, max(0, time() - start_time))
    end
    return LabelingSolution(goal_labels, LS_OPTIMAL, goal_costs, max(0, time() - start_time))
end
