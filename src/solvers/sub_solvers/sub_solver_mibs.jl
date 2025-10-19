# It is essentially the MIP subsolver, but it uses MiBS as the solver. The new implementation makes it easier to debug because we need to manage two models within the BilevelJuMP instance.
## It is intended for use when solving a bilevel setting within the Lagrangian subproblem. 
## NOTE that MiBS is still an experimental feature, and does not support many of the desired features (e.g, time limit setting)

using JuMP, BilevelJuMP, MibS_jll


"""
An defualt implementation of an SubSolver that uses an underlying JuMP model. 
Mainly intendendent for MIP as subsolver.
Note that the sub_problem is modified in the process. 
"""
struct SubSolverMiBS{T} <: SubSolver
    name::String  # the name (unique identifier) of this sub_problem
    bi_model::BilevelModel
    A::Vector{T}  # Iterable of resources
    y_vars::Any  # Dict of sub_problem (lower-level) variables appearing in the interdiction constraints. Keys=A. NOTE: We essentially assume that y-vars are binary, and we do not make any copies by default
    r_objterm::GenericAffExpr  # The objective function term within the master problem, i.e., the contribution to the cost in master (should be generated with sub_problem variables)
    c_objterm::GenericAffExpr  # the objective function of the original sub_problem (should be generated with sub_problem variables)
end

####### Implementation of Solver required functions #######
function capacity_linking(sol::SubSolverMiBS, a, params::SolverParam)
    return 1
end

function check(sol::SubSolverMiBS, params::SolverParam)
    # check if containers for linking vars/constraints have same size
    if !(length(sol.y_vars) == length(sol.A))
        error(
            "In subsolver $(name(sol)), the number of second level variables in linking constraints does not match the number of resources. 
      Please remeber that for each resource you need to have a linking constraint and corresponding variables! |A| = $(length(sol.A)) and |y_vars| = $(length(sol.y_vars)).",
        )
    end

    # check if linking variables and second level variables in linking constraints are contained in model
    for a in sol.A
        if !is_valid(Lower(sol.bi_model), sol.y_vars[a])
            error(
                "Second level variable (y) for resource $(a) is not contained in subsolver lower level $(name(sol))!",
            )
        end
    end

    # checks finished
end

function compute_lower_bound_master_contribution(sol::SubSolverMiBS, params::SolverParam, time_limit)
    # set objective function corresponding to master problem contribution
    @objective(Upper(sol.bi_model), Min, sol.r_objterm)

    # print the MIP to solve when debbug mode
    try 
        if should_debbug_print(params)
            write_to_file(
                sol.bi_model,
                "$(params.output_folder_path)/ubbound_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    # solve adjusted sub 
    solution = _solve_mip(sol, params.debbug_out)

    # in debbug mode, print solution before return
    @debug "The found lower bound solution for $(sol.name) is $(solution.objective)."
    return solution.objective
end

function name(ss::SubSolverMiBS)
    return ss.name
end

function separation!(sol::SubSolverMiBS, sval, gvals, kvals::Dict, params::SolverParam, time_limit)
    # TODO: we use hard coded numeric tolerances here...
    obj_tolerance = 4 # 10e-4
    var_non_zero_tolerance = 10e-4

    # set new objective function
    k_term = sum([kvals[a] * sol.link_varsC[a] for a in sol.A])  # TODO: This is only correct if linking variables on sub_problem side are binary? I think, it should be correct like this, but it is at least not very intuitive. 
    new_obj = sol.c_objterm * gvals + sol.r_objterm + k_term
    @objective(Upper(sol.bi_model), Min, new_obj)

    # Solve sub_problem (and print it in debbug mode)
    try 
        if should_debbug_print(params)
            write_to_file(
                sol.bi_model,
                "$(params.output_folder_path)/sub_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end
    solution = _solve_mip(sol, params.debbug_out)

    # check solution status
    if !solution.status 
        @error "When solving subproblem $(sub_solver.name) for separating GBC cuts, the Lagrangian subproblem was reported to be infeasible but MiBS."
    end

    # construct return value 
    opt_obj = solution.objective
    @debug "Optimal solution of sub " * sol.name * " is " * string(opt_obj)

    is_violate = Bool(!(sval < opt_obj + var_non_zero_tolerance))
    @debug "The violated status of sub $(sol.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    r = eval_affexpr_by_name(sol.r_objterm, sol.bi_model, solution.all_lower)
    c = eval_affexpr_by_name(sol.c_objterm, sol.bi_model, solution.all_lower)
    as = indices_where_coeff_pos(solution.all_lower, sol.y_vars)
    return SubSolution(is_violate, r, c, as)
end


function separation_BlC!(sub_solver::SubSolverMiBS, sval, kvals::Dict, params::SolverParam, time_limit)
    # TODO: we use hard coded numeric tolerances here...
    obj_tolerance = 4 # 10e-4
    var_non_zero_tolerance = 10e-4

    # set new objective function
    k_term = sum([kvals[a] * (1-sub_solver.y_vars[a]) for a in sub_solver.A]) 
    new_obj = sub_solver.c_objterm - k_term
    @objective(Upper(sub_solver.bi_model), Max, new_obj)

    # Solve sub_problem (and print it in debbug mode)
    try 
        if should_debbug_print(params)
            write_to_file(
                sub_solver.bi_model,
                "$(params.output_folder_path)/sub_$(sub_solver.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sub_solver.name) to file. error message $err"
    end
    solution = _solve_mip(sub_solver, params.debbug_out)

    # check solution status
    if !solution.status 
        @error "When solving subproblem $(sub_solver.name) for separating BlC cuts, the Lagrangian subproblem was reported to be infeasible but MiBS."
    end

    # construct return value 
    opt_obj =solution.objective
    @debug "Optimal solution of sub " * sub_solver.name * " is " * string(opt_obj)

    is_violate = Bool(!(sval > opt_obj - var_non_zero_tolerance))
    @debug "The violated status of sub $(sub_solver.name) is $(is_violate) as sval=$(sval) and found obj=$(opt_obj)."
    # TODO: this rounding can be dangerous (but without code also breaks)...
    r = eval_affexpr_by_name(sol.r_objterm, sol.bi_model, solution.all_lower)
    c = eval_affexpr_by_name(sol.c_objterm, sol.bi_model, solution.all_lower)
    as = indices_where_coeff_pos(solution.all_lower, sol.y_vars; non_used = true)  # Note that we need the ressources that were NOT used
    #@debug "The x-solution of sub $(sub_solver.name) is $([(a, value(sub_solver.y_vars[a])) for a in sub_solver.A])." 
    return SubSolution(is_violate, r, c, as)
end

function set_nthreads(sol::SubSolverMiBS, n)
    #set_attribute(sol.bi_model, MOI.NumberOfThreads(), n)
    @warn "We currently not diretly support setting the number of threads for MiBS solver as MiBS is strangly accessed over a separate function and not a solver interface."
end

function solve_sub_for_x(sol::SubSolverMiBS, xvals, params::SolverParam, time_limit)
    # set original objective function (just to make sure)
    @objective(Upper(sol.bi_model), Min, sol.c_objterm)

    # fix values for linking variables
    ## we must set the constraints to lower level as MiBS does not support 
    @constraint(Lower(sol.bi_model), fixc[a=sol.A], sol.y_vars[a] == round(xvals[a]))  # The rounding seems to help with numerics

    # solve adjusted sub (and write to file in debbug mode)
    try
        if should_debbug_print(params)
            write_to_file(
                sol.bi_model,
                "$(params.output_folder_path)/subfix_$(sol.name).$(params.file_format_output)",
            )
        end
    catch err
        @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
    end

    # set adjusted time limit
    @debug "Subproblem $(sol.name) MIP was adjusted. Start MIP solver to solve it. "
    solution = _solve_mip(sol, params.debbug_out)

    try
        # here, we can have infeeasible solutions due to wrong first-level decision. Catch this case
        
        if !solution.status
            @debug "The Subproblem $(sol.name) was infeasible"
            return false, 0, Dict()
        end

        # evaluate solution
        y_vals = index_value_map(y_vals, solution.all_lower)
        osol = solution.objective
        @debug "We found an solution for sub_problem $(sol.name) with value $(osol)."

        # if debbug mode, print model after clean up
        try
            if should_debbug_print(params)
                write_to_file(
                    sol.bi_model,
                    "$(params.output_folder_path)/subfixCLEAN_$(sol.name).$(params.file_format_output)",
                )
            end
        catch err
            @error "Could not print Submodel MIP $(sol.name) to file. error message $err"
        end

        # evaluate first-level objective
        osol_L1 = eval_affexpr_by_name(sol.r_objterm, sol.bi_model, solution.all_lower)

        # return found solution
        return true, osol, osol_L1, y_vals
    finally
        # cleanup
        for a in sol.A # TODO: can we implement this faster?
            delete(Lower(sol.bi_model), fixc[a])
        end
        unregister(Lower(sol.bi_model), :fixc)
    end
end



###            Auxiliary Functions              ###

"""
Evaluate the passed AffExpr based on the values in the passed Dict. Dict keys must be variable names.
For expressions with variables not present in your dict, decide whether to:
    - Treat them as zero (set require_all = false), or
    - Throw an error to ensure completeness (require_all = true).
"""
function eval_affexpr_by_name(A::JuMP.AffExpr, model::JuMP.Model, values::Dict{String,<:Real}; require_all::Bool = true)
    s = JuMP.constant(A)
    # TODO: we use hard coded numeric tolerances here...
    obj_tolerance = 4 # 10e-4
    # Add contributions from the variables for which you provided values
    for (nm, val) in values
        vref = JuMP.variable_by_name(model, nm)
        vref === nothing && error("No variable named '$nm' in the model.")
        s += JuMP.coefficient(A, vref) * round(float(val); digits=obj_tolerance)
    end
    if require_all
        # Optional: ensure all variables in A have values
        for v in JuMP.all_variables(model)
            c = JuMP.coefficient(A, v)
            if c != 0.0 && !haskey(values, JuMP.name(v))
                error("Missing value for variable $(JuMP.name(v)) used in the expression.")
            end
        end
    end
    return s
end

"""
Extract the ressources for which the solution Dict reports values greater that 0 (i.e., we expect the binary variables to be 1).
    if 'nonused' is true, instead filter all elements that were set to 0.
"""
function indices_where_coeff_pos(values::Dict{String,<:Real}, B::AbstractDict{I,JuMP.VariableRef}; non_used=false) where {I}
    # TODO: we use hard coded numeric tolerances here...
    varnonzero_tolerance = 0.0004 # 10e-4
    result = I[]
    for (i, v) in B
        nm = JuMP.name(v)
        if haskey(values, nm)
            val = values[nm]
            if val < -varnonzero_tolerance @error " The value of variable $nm was negative in MiBS solution, but variable should be binary!" end
            if non_used
                if val < 1 - varnonzero_tolerance
                    push!(result, i)
                end
            else
                if val > varnonzero_tolerance
                    push!(result, i)
                end
            end
        else
            if missing === :error
                error("No value provided for variable '$nm'.")
            elseif missing === :zero
                if 0.0 > c
                    push!(result, i)
                end
            elseif missing === :skip
                # do nothing

            else
                error("Unknown missing policy: $missing")
            end
        end
    end
    return result
end

"""
Return a Dict mapping index => value for the variables, using a Dict of values keyed by VariableRef.
"""
function index_value_map(B::AbstractDict{I,JuMP.VariableRef}, values::AbstractDict{String,<:Real}) where {I}
    # TODO: we use hard coded numeric tolerances here...
    obj_tolerance = 4 # 10e-4
    result = Dict{I, Float64}()
    for (i, v) in B
        nm = JuMP.name(v)
        if haskey(values, nm)
            result[i] = round(float(values[nm]); digits=obj_tolerance)
        else
            if missing === :skip
                continue
            elseif missing === :error
                error("No value provided for variable '$nm'.")
            elseif missing === :default
                result[i] = float(default)
            else
                error("Unknown missing policy: $missing")
            end
        end
    end
    return result
end

"""
    _solve_mip(sol::SubSolverMiBS, printlog::Bool)

Solve the underlying Bilevel JuMP Model using MiBS. Returns the 'solution' struct of the 'solve_with_MibS' function.
"""
function _solve_mip(sol::SubSolverMiBS, printlog::Bool)
    # set adjusted time limit
    @debug "Starting MiBS subsolver $(sol.name)"
    #set_time_limit_sec(sol.bi_model, time_limit)

    # solve Bilevel JuMP model with MiBS 
    solution = BilevelJuMP.solve_with_MibS(sol.bi_model, MibS_jll.mibs; verbose_results=printlog, keep_files=printlog)
    return solution
end






