using JuMP


"""
An defualt implementation of an SubSolver that uses an underlying JuMP model. 
Mainly intendendent for MIP as subsolver.
Note that the sub_problem is modified in the process. 
"""
struct SubSolverJuMP{T} <: SubSolver
    name::String  # the name (unique identifier) of this sub_problem
    mip_model::JuMP.Model
    A::Vector{T}  # Iterable of resources
    link_varsC::Any  # A list of linking variables (copy within sub_problem). Keys=A
    y_vars::Any  # Dict of sub_problem variables appearing in the interdiction constraints. Keys=A
    link_constraints_capacities::Dict{T,<:Number}  # the capacity parameters in the interdiction linking constraints. I would recommend to set it to 1 und use y_vars as indicator variables sub_problem formulation!
    r_objterm::AffExpr  # The objective function term within the master problem, i.e., the contribution to the cost in master (should be generated with sub_problem variables)
    c_objterm::AffExpr  # the objective function of the original sub_problem (should be generated with sub_problem variables)

    # extra functionalities
    extra_cuts::Function # Adds additional cuts in a branch and check manner. Called after 'mip_model' is solved to optimality. Should add the cut to the MIP itself. Input is only the time limit for this function. Then MIP is resolved. Return (b, time) where b is Bool if resolving should be used and time the runtime spend within the function
    # TODO: implement MIPGap for GBC solver (i.e., solve subproblem heuristically as long as it is not the coefficient of Benders-like cut)?
end

# TODO: It would maybe be much nicer to have a constructor who takes some struct representing the linking constraints that having them splitt between 3 different inputs

function (::Type{SubSolverJuMP{T}})(
    name, mip_model, A,
    link_varsC, y_vars, link_constraints_capacities,
    r_objterm, c_objterm
) where T
    SubSolverJuMP{T}(
      name, mip_model, A,
      link_varsC, y_vars, link_constraints_capacities,
      r_objterm, c_objterm,
      timelimit -> false, 0
    )
end

####### Generic public functions that support with implementing SubSolverJuMP #######
"""
    extra_cuts_benderslike_JuMP(jump::JuMP.Model, time_limit)

You can pass a JuMP model to this function that solves the second-level problem. It will act as an 'extra_cuts' function for the generation of additional Bender-like cuts for the Subsolver MIP. 
    
# Arguments
    - 'jump::JuMP.Model': The JuMP model that solves the second-level subproblem.
    - 'sub::JuMP.Model': The MIP of the SubSolverJuMP
    - 'A': The set of resources.
    - 'oL2': The second-level objective 
    - 'x_subproblem': The copies of the bilevel linking variables within the subproblem MIP
    - 'x_jump': The copies of the 'x_subproblem' variables within the 'jump' model.
    - 'y_jump': The variables in the 'jump' model that are indicted by the corresponding 'x_jump' variables.
    - 'bigMs::Function': The bigM function for generating coefficients for Benders-like cuts. Should take a as argument where a is the resource.
    - 'saved_cuts::Dict': A Dict in which we save all the generated cuts. Cuts are saved as tuples (rhs of blc with bigM terms, value of second level solution)
    - 'time_limit': The time limit that should not be exceeded by this function. Note that in case of time out an 'TimeoutException'
"""
function extra_cuts_benderslike_JuMP(jump::JuMP.Model, sub::JuMP.Model, A, oL2, x_subproblem, x_jump, y_jump, bigMs::Function, saved_cuts::Dict, time_limit)
    current_L2solution_sub = JuMP.value(oL2)
    x_vals = Dict(a => round(JuMP.value(x_subproblem[a])) for a in A)  # round to integer 0 or 1 values

    msolkey = JuBiC.key_master_sol(x_vals, A)
    if haskey(saved_cuts, msolkey)
    # if we already solved this sub_problem, we do not need to solve it again as we already added corresponding cut
    subopt = saved_cuts[msolkey][2]
    need_cut = (current_L2solution_sub > subopt)
        if need_cut
            @debug "Need to recover cut with rhs $subopt while current master solution is $current_L2solution_sub"
            cutopt = @constraint(sub, oL2 <= saved_cuts[msolkey][1]) 
            @debug "Recovered to Subsolver the additional constraint $cutopt"
        else
            @debug "We recovered cut but no need to add it as cut rhs is $subopt while current master solution is $current_L2solution_sub"
            return false, 0  # we ignore runtime for the check for internal runtime tracking
        end
    end

    @debug "Starting Benders-like cuts within Subsolver MIP by solving helper MIP for sub solution $x_vals."
    # set time limit for solving
    set_time_limit_sec(jump, time_limit)
    # fix values for linking variables
    @constraint(jump, fixc[a=A], x_jump[a] == round(x_vals[a]))
    # now, solve the helper MIP
    optimize!(jump)

    # handle solution 
    solve_time = JuMP.solve_time(jump)
    status = JuMP.termination_status(jump)
    @debug "Helper MIP was solved with status $status"
    if status == MOI.OPTIMAL
        y_vals = JuMP.value.(y_jump)
        subopt = JuMP.objective_value(jump)
    elseif status == MOI.TIME_LIMIT
        throw(TimeoutException("We hit the time limit within 'extra_cuts_benderslike_JuMP' function when trying to find the new cut to separate."))
    else
        error("Helper MIP terminated with status $status")
    end

    # generate cut coefficients
    bigMterms = 0
    for a in A
        bigMterms +=
            bigMs(a) *
            y_vals[a] *
            (1 - x_subproblem[a])
    end

    # save found cut
    saved_cuts[msolkey] = (subopt + bigMterms, subopt)

    # clean up for next run by removing constraints enforcing model properties 
    @debug "clean up helper model for next callback call"
    for a in A
        delete(jump, fixc[a])
    end
    unregister(jump, :fixc)

    # check if cut is needed to be added
    @debug "Solved helper MIP and found solution with value $subopt while current master solution is $current_L2solution_sub"
    need_cut = (current_L2solution_sub > subopt)
    if need_cut
        cutopt = @constraint(sub, oL2 <= subopt + bigMterms) 
        @debug "Added to Subsolver the additional constraint $cutopt"
                
        # finished adding cut. Check if new cut implies resolving and then return
        return true, solve_time # always resolve as we added a constraint. Time to build cut neglected in internal runtime computations
    else
        @debug "No need to add cut as current Subsolver solution is valid"
        return false, solve_time
    end
end


####### Implementation of Solver required functions #######
function capacity_linking(sol::SubSolverJuMP, a, params::SolverParam)
    return sol.link_constraints_capacities[a]
end

function check(sol::SubSolverJuMP, params::SolverParam)
    # check objective sense is minimization
    if !(objective_sense(sol.mip_model) == MIN_SENSE)
        error("The subsolver $(name(sol)) must be formulated as an minimization problem.")
    end

    # if output is turned on, write a warning that it will slow down the solver
    try
        if get_optimizer_attribute(sol.mip_model, "OutputFlag") == 1
            printstyled(
                "Currently, the subsolver $(name(sol)) will create output to console or file. 
    It is highly recommended to set model to silent to improve computational performance. \n ";
                color=:orange,
            )
        end
    catch e
        printstyled(
            "We encured an error when checking output status of subsolver $(name(sol)). Please do it manually! \n ";
            color=:red,
        )
        println(e)
    end

    # check if containers for linking vars/constraints have same size
    if !(length(sol.link_constraints_capacities) == length(sol.link_varsC))
        error(
            "In subsolver $(name(sol)), the capacities in linking constraints do not match the (copied) master linking variables:
       |capa| = $(length(sol.link_constraints_capacities)) and |link_vars| = $(length(sol.link_varsC)).",
        )
    elseif !(length(sol.link_varsC) == length(sol.y_vars))
        error(
            "In subsolver $(name(sol)), the number of linking variables does not match the number of second level variables in linking constraints! 
      |link_vars| = $(length(sol.link_varsC)) and |y_vars| = $(length(sol.y_vars)).",
        )
    elseif !(length(sol.y_vars) == length(sol.A))
        error(
            "In subsolver $(name(sol)), the number of second level variables in linking constraints does not match the number of resources. 
      Please remeber that for each resource you need to have a linking constraint and corresponding variables! |A| = $(length(sol.A)) and |y_vars| = $(length(sol.y_vars)).",
        )
    end

    # check if linking variables and second level variables in linking constraints are contained in model
    for a in sol.A
        if !is_valid(sol.mip_model, sol.link_varsC[a])
            error(
                "Linking variable for resource $(a) is not contained in subsolver $(name(sol))!",
            )
        elseif !is_valid(sol.mip_model, sol.y_vars[a])
            error(
                "Second level variable (y) for resource $(a) is not contained in subsolver $(name(sol))!",
            )
        end
    end

    # checks finished
end

function compute_lower_bound_master_contribution(sol::SubSolverJuMP, params::SolverParam, time_limit)
    # set objective function corresponding to master problem contribution
    @objective(sol.mip_model, Min, sol.r_objterm)

    # print the MIP to solve when debbug mode
    try 
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/ubbound_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    # solve adjusted sub 
    solve_mip(sol, time_limit)

    # check if optimal solution found and otherwise handle exceptions
    check_solution_status(sol)

    # in debbug mode, print solution before return
    @debug "The found lower bound solution for $(sol.name) is $(objective_value(sol.mip_model))."
    return objective_value(sol.mip_model)
end

function name(ss::SubSolverJuMP)
    return ss.name
end

function separation!(sol::SubSolverJuMP, sval, gvals, kvals::Dict, params::SolverParam, time_limit)
    # TODO: we use hard coded numeric tolerances here...
    obj_tolerance = 4 # 10e-4
    var_non_zero_tolerance = 10e-4

    # set new objective function
    k_term = sum([kvals[a] * sol.link_varsC[a] for a in sol.A])  # TODO: This is only correct if linking variables on sub_problem side are binary? I think, it should be correct like this, but it is at least not very intuitive. 
    new_obj = sol.c_objterm * gvals + sol.r_objterm + k_term
    @objective(sol.mip_model, Min, new_obj)

    # Solve sub_problem (and print it in debbug mode)
    try 
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/sub_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end
    solve_mip(sol, time_limit)

    # check if optimal solution found and otherwise handle exceptions
    check_solution_status(sol)

    # construct return value 
    opt_obj = objective_value(sol.mip_model)
    @debug "Optimal solution of sub " * sol.name * " is " * string(opt_obj)

    is_violate = Bool(!(sval < opt_obj + var_non_zero_tolerance))
    @debug "The violated status of sub $(sol.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    # TODO: this rounding can be dangerous (but without code also breaks)...
    r = round(value(sol.r_objterm); digits=obj_tolerance)
    c = round(value(sol.c_objterm); digits=obj_tolerance)
    as = [a for a in sol.A if value(sol.y_vars[a]) > var_non_zero_tolerance]  
    return SubSolution(is_violate, r, c, as)
end

function set_nthreads(sol::SubSolverJuMP, n)
    set_attribute(sol.mip_model, MOI.NumberOfThreads(), n)
end

function solve_sub_for_x(sol::SubSolverJuMP, xvals, params::SolverParam, time_limit)
    # set original objective function (just to make sure)
    @objective(sol.mip_model, Min, sol.c_objterm)

    # fix values for linking variables
    @constraint(sol.mip_model, fixc[a=sol.A], sol.link_varsC[a] == round(xvals[a]))  # TODO: The rounding seems to help with numerics
    #@constraint(sol.mip_model, fixc[a=sol.A], sol.link_varsC[a] == xvals[a])

    # solve adjusted sub (and write to file in debbug mode)
    try
        if should_debbug_print(params)
            write_to_file(
                sol.mip_model,
                "$(params.output_folder_path)/subfix_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    # set adjusted time limit
    @debug "Subproblem $(sol.name) MIP was adjusted. Start MIP solver to solve it. "
    solve_mip(sol, time_limit)

    try
        # here, we can have infeeasible solutions due to wrong first-level decision. Catch this case
        status = termination_status(sol.mip_model)
        if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
            @debug "The Subproblem $(sol.name) was infeasible"
            return false, 0, Dict()
        end

        # check if optimal solution found and otherwise handle exceptions
        check_solution_status(sol)

        # evaluate solution
        y_vals = value.(sol.y_vars)
        osol = objective_value(sol.mip_model)
        @debug "We found an solution for sub_problem $(sol.name) with value $(osol)."

        # if debbug mode, print model after clean up
        try
            if should_debbug_print(params)
                write_to_file(
                    sol.mip_model,
                    "$(params.output_folder_path)/subfixCLEAN_$(sol.name).$(params.file_format_output)",
                )
            end
        catch err
            @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
        end

        # return found solution
        return true, osol, y_vals
    finally
        # cleanup
        for a in sol.A # TODO: can we implement this faster?
            delete(sol.mip_model, fixc[a])
        end
        unregister(sol.mip_model, :fixc)
    end
end



####### Auxiliary functions #######
"""
    check_solution_status()

Check if model terminated due to finding an optimal solution. If not, throw an error. For timeouts, throw 'TimeoutException'.
"""
function check_solution_status(sol::SubSolverJuMP)
    # check for time limit
    status = termination_status(sol.mip_model)
    if status == MOI.TIME_LIMIT
        @debug "We hit the time limit while solving MIP of subsolver $(sol.name)"
        throw(TimeoutException("We hit the time limit while solving MIP of subsolver $(sol.name)"))
    end

    # check for infeasibility (or other non desired status)
    status = termination_status(sol.mip_model)
    if status == MOI.INFEASIBLE
        @debug "The sub_problem $(sol.name) is infeasible for the current master solution. Starting computations of IIS"
        #unset_silent(sol.mip_model)  # do not do this as it breaks the IIS computations somehow
        compute_conflict!(sol.mip_model)
        @debug "IIS computed, now printing it"
        iis_model, _ = copy_conflict(sol.mip_model)
        if get_attribute(sol.mip_model, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
            @debug iis_model # printing to file just causes error...
        else
            @debug "The IIS was not computed succesfully??? Maybe numerics??"
        end
        error("Submodel $(sol.name) was infeasible. Computed IIS but stopping solution process (as it is cleary a bug)")
    elseif status == MOI.DUAL_INFEASIBLE
        error("The second level was unbounded what violates assumption on algorithm.")
    end

    if !(status == MOI.OPTIMAL)
        @warn "The following error is most likely due to the user given function adding a constraint to MIP but forcing the solver to resolve the problem afterwards."
        error("Submodel $(sol.name) could not be solver to optimality but terminated with status $status")
    end
end


"""
    solve_mip(sol::SubSolverJuMP, time_limit)

Solve the underlying MIP. Mainly calls 'optimize!' and handles Branch&Check for additional cuts. It also sets the time limit and adjusts it if necessary.
"""
function solve_mip(sol::SubSolverJuMP, time_limit)
    inner_time_limit = time_limit
    need_solving = true

    while need_solving
        need_solving = false

        if inner_time_limit <= 0
            @debug "We hit the time limit while solving MIP of subsolver $(sol.name) in the branch-and-cut framework."
            throw(TimeoutException("We hit the time limit while solving MIP of subsolver $(sol.name) in the branch-and-cut framework."))
        end

        # set adjusted time limit
        @debug "Setting time limit of subsolver $(sol.name) to $inner_time_limit"
        set_time_limit_sec(sol.mip_model, inner_time_limit)

        # solve MIP 
        optimize!(sol.mip_model)
        inner_time_limit = inner_time_limit - solve_time(sol.mip_model)

        # if solved to optimality, apply branch and check
        status = termination_status(sol.mip_model)
        @debug "Subproblem $(sol.name) MIP was solved with status $status"
        if status == MOI.OPTIMAL
            need_solving, time_spend = sol.extra_cuts(time_limit) 
            inner_time_limit = max(inner_time_limit - time_spend, 0) # deduce time spend in subsolver from time limit
            if need_solving  
                @debug "User provided additional cut for branch-and-check of Subsolver MIP $(name(sol)). Initializing resolving of MIP."
            else
                @debug "User function checked model in branch-and-cut and decided that no reoptimization is needed."
            end
        end
    end
end

