using JuMP

### Constructor and basic functions ###
struct ConSubsolCut
    res::AbstractVector  # List of ressources used by solution
    objL2::Number  # L2 solution value, i.e, cost
    objL1::Number  # L1 solution value, i.e., risk
end

function Base.:(==)(scs1::ConSubsolCut, scs2::ConSubsolCut)
    if scs1.objL2 == scs2.objL2 && scs1.objL1 == scs2.objL1
        for r1 in scs1.res
            if !(r1 in scs2.res)
                return false
            end
        end
        for r2 in scs2.res
            if !(r2 in scs1.res)
                return false
            end
        end
        return true
    end
    return false
end

"""
The LP Benders solver for a single sub_problem. 
"""
struct ConnectorLP{T}
    lp::Model  # The lp model. Variables should be registered as :s, :k[A], :g
    A::Vector{T}  # Iterable of resources
    link_vars::Dict{T,VariableRef}  # The linking variables (from master) as dict. Keys=A
    sub_solver::SubSolver  # The sub_solver that is used to solve the sub_problem
    lower_bound_obj_contribution::Number  # The lower bound on the contribution to the first level objective function

    blc_cut_generator::Union{ConnectorLP_BlC, Nothing}  # subroutine for generating BlC cut coefficients required for better GBC cut coefficients
    my_subsolutions::AbstractVector{ConSubsolCut}  # the storage of found solutions of subproblem used to generate cuts. Should be init with empty list

    g_round_digits::Int  # we round obtained g-value to this digits (see 'ceil' docu of digits). It ensures numerical stability while weakening cuts.
    numeric_state::Dict{Symbol,Any}  # stores temporary numerical diagnostics/overrides for the current connector solve
end

function name(clp::ConnectorLP)
    return name(clp.sub_solver)
end

function _reset_numeric_state!(subLP::ConnectorLP)
    empty!(subLP.numeric_state)
    return nothing
end

function _connector_solution_snapshot(subLP::ConnectorLP)
    return Dict{Symbol,Any}(
        :s => Float64(value(subLP.lp[:s])),
        :g => Float64(value(subLP.lp[:g])),
        :k => Dict(a => Float64(value(subLP.lp[:k][a])) for a in subLP.A),
    )
end

function _snapshot_s(subLP::ConnectorLP)
    snapshot = get(subLP.numeric_state, :cut_solution_snapshot, nothing)
    return isnothing(snapshot) ? value(subLP.lp[:s]) : snapshot[:s]
end

function _snapshot_g(subLP::ConnectorLP)
    if haskey(subLP.numeric_state, :g_override)
        return subLP.numeric_state[:g_override]
    end
    snapshot = get(subLP.numeric_state, :cut_solution_snapshot, nothing)
    return isnothing(snapshot) ? value(subLP.lp[:g]) : snapshot[:g]
end

function _snapshot_k(subLP::ConnectorLP)
    snapshot = get(subLP.numeric_state, :cut_solution_snapshot, nothing)
    return isnothing(snapshot) ? value.(subLP.lp[:k]) : snapshot[:k]
end

function _pareto_row_scale(subLP::ConnectorLP, current_obj::Real, g_obj_coef::Real)
    cssol = Float64(value(subLP.lp[:s]))
    cgsol = Float64(value(subLP.lp[:g]))
    kval_sum = sum(abs(Float64(value(subLP.lp[:k][a]))) for a in subLP.A)
    return max(
        1.0,
        abs(cssol),
        abs(g_obj_coef * cgsol),
        kval_sum,
        abs(current_obj),
        abs(subLP.lower_bound_obj_contribution),
    )
end

function _pareto_band_tolerance(
    subLP::ConnectorLP,
    current_obj::Real,
    g_obj_coef::Real,
    params::GBCparam,
)
    return max(0.0, Float64(params.pareto_band_tolerance))
end

function _mip_gap_tolerances(model::JuMP.Model)
    rel_gap = try
        Float64(get_optimizer_attribute(model, "MIPGap"))
    catch
        1e-4
    end
    abs_gap = try
        Float64(get_optimizer_attribute(model, "MIPGapAbs"))
    catch
        1e-10
    end
    return rel_gap, abs_gap
end

function _subsolver_gap_tolerances(sub_solver)
    model = try
        getproperty(sub_solver, :mip_model)
    catch
        nothing
    end

    if isnothing(model) || !(model isa JuMP.Model)
        return 1e-4, 1e-10
    end
    return _mip_gap_tolerances(model)
end

function _estimate_g_rounding_tolerance(
    subLP::ConnectorLP,
    sub_solver::SubSolution,
)
    alpha = Float64(sub_solver.obj_second_level)
    iszero(alpha) && return 0.0, 0.0, 0.0

    sval = Float64(value(subLP.lp[:s]))
    gval = Float64(value(subLP.lp[:g]))
    kvals = Dict(a => Float64(value(subLP.lp[:k][a])) for a in subLP.A)
    rhs = Float64(sub_solver.obj_first_level)
    opt_obj = Float64(sub_solver.obj_compare)

    row_scale = max(
        1.0,
        abs(sval),
        sum(abs(kvals[a]) for a in subLP.A),
        abs(alpha * gval),
        abs(rhs),
    )
    delta_g_connector = 1e-6 * row_scale / abs(alpha)

    rel_gap, abs_gap = _subsolver_gap_tolerances(subLP.sub_solver)
    sub_obj_scale = max(1.0, abs(opt_obj))
    delta_g_sub = max(abs_gap, rel_gap * sub_obj_scale) / abs(alpha)

    delta_g_total = max(delta_g_connector, delta_g_sub)
    return delta_g_connector, delta_g_sub, delta_g_total
end

function _required_g_value(
    subLP::ConnectorLP,
    sub_solver::SubSolution,
)
    alpha = Float64(sub_solver.obj_second_level)
    iszero(alpha) && return Inf
    sval = Float64(value(subLP.lp[:s]))
    rhs = Float64(sub_solver.obj_first_level)
    ksum = sum(Float64(value(subLP.lp[:k][a])) for a in sub_solver.A_sub)
    return (sval - ksum - rhs) / alpha
end

function _duplicate_cut_warning_message(
    subLP::ConnectorLP,
    csc::ConSubsolCut,
    sval::Real,
    opt_obj::Real,
    g_current::Real,
    g_required::Real,
    g_rounded::Real,
    delta_g_connector::Real,
    delta_g_sub::Real,
    delta_g_total::Real,
)
    resource_text = isempty(csc.res) ? "no resources" : join(string.(csc.res), ", ")
    return "JuBiC detected a potentially numerically unstable repeated cut while solving ConnectorLP $(name(subLP.sub_solver)). " *
           "The same cut for resource pattern [$(resource_text)] was generated again. " *
           "Current connector value: s=$(sval). Current follower value: $(opt_obj). " *
           "Current g=$(g_current), required g to satisfy the cut=$(g_required), numerically rounded g=$(g_rounded). " *
           "Estimated g tolerances: connector=$(delta_g_connector), subproblem=$(delta_g_sub), combined=$(delta_g_total). " *
           "Under these tolerances, JuBiC will not add the same cut again and treats the current connector solution as numerically optimal."
end

function _warn_cut_overlap!(
    subLP::ConnectorLP,
    params::GBCparam,
    a,
    overlap_kind::AbstractString,
    overlap_value,
    kval,
)
    @warn "ConnectorLP $(name(subLP.sub_solver)) generated a GBC cut where the same resource appears in both the k-term and the $(overlap_kind)-term. JuBiC continues, but this indicates numerically delicate or unintuitive behavior and the resulting solution should be double-checked. Resource=$(a), $(overlap_kind)=$(overlap_value), k=$(kval)."
    return nothing
end

function _record_opt_cut_validation!(
    subLP::ConnectorLP,
    params::GBCparam,
    incumbent_cut_value::Real,
    incumbent_reference_value::Real,
    cut_rhs::Real,
)
    if params.integer_obj
        cut_integer = round(Float64(incumbent_cut_value))
        expected_integer = round(Float64(incumbent_reference_value))
        cut_is_integer = abs(Float64(incumbent_cut_value) - cut_integer) <= 1e-8
        expected_is_integer = abs(Float64(incumbent_reference_value) - expected_integer) <= 1e-8
        if !(cut_is_integer && expected_is_integer && cut_integer == expected_integer)
            throw(NumericalIssueException(
                "ConnectorLP $(name(subLP.sub_solver)) generated an integer-object optimality cut that does not evaluate exactly at the incumbent x-solution. cut(x)=$(incumbent_cut_value), incumbent_reference=$(incumbent_reference_value), cut_rhs=$(cut_rhs), rounded_cut=$(cut_integer), rounded_reference=$(expected_integer).",
                "Terminate_Numerics",
            ))
        end
        return nothing
    end

    diff = Float64(incumbent_cut_value - incumbent_reference_value)
    tol = 1e-4
    if abs(diff) > tol
        @warn "ConnectorLP $(name(subLP.sub_solver)) generated an optimality cut that does not reproduce the incumbent reference value implied by the cut structure. cut(x)=$(incumbent_cut_value), incumbent_reference=$(incumbent_reference_value), cut_rhs=$(cut_rhs), diff=$(diff)."
        if haskey(params.stats.data, "NOptCutValidationWarnings")
            add_stat!(params.stats, "NOptCutValidationWarnings", 1)
        else
            new_stat!(params.stats, "NOptCutValidationWarnings", 1)
        end
        offenders = get!(params.stats.data, "OptCutValidationUsers", String[])
        if !(name(subLP.sub_solver) in offenders)
            push!(offenders, name(subLP.sub_solver))
        end
    end
    return nothing
end

function _round_cut_component(value::Real, params::GBCparam)
    if params.integer_obj
        return Float64(ceil(value))
    end
    return Float64(value)
end

function _sanitize_nonnegative_opt_cut_coefficient(
    value::Real,
    coeff_name::AbstractString,
    subLP::ConnectorLP;
    tol::Float64=1e-3,
)
    numeric_value = Float64(value)
    if numeric_value >= 0.0
        return numeric_value
    end
    if numeric_value >= -tol
        @debug "ConnectorLP $(name(subLP.sub_solver)) rounds the $(coeff_name) coefficient $(numeric_value) up to 0.0 within numerical tolerance $(tol)."
        return 0.0
    end
    throw(ArgumentError("The $(coeff_name) coefficient $(numeric_value) used in ConnectorLP $(name(subLP.sub_solver)) is negative!"))
end

function _compute_opt_cut_bound_coefficient(subLP::ConnectorLP, sval, optL2, gval, kvals, x_vals)
    bound_value = sval - optL2 * gval - subLP.lower_bound_obj_contribution
    for a in subLP.A
        if Float64(kvals[a]) > 0.0 && Float64(x_vals[a]) > 0.0
            @debug "ConnectorLP $(name(subLP.sub_solver)) detected k[$(a)]=$(kvals[a]) with x[$(a)]=$(x_vals[a]) > 0 while computing xi=$(bound_value). This indicates an edge case where the cut could potentially be strengthened further by exploiting the simultaneous activity of k and x on arc $(a)."
        end
    end
    @debug "ConnectorLP $(name(subLP.sub_solver)) uses bound coefficient xi=$(bound_value) with s=$(sval), optL2=$(optL2), g=$(gval), and lower_bound=$(subLP.lower_bound_obj_contribution)."
    return _sanitize_nonnegative_opt_cut_coefficient(bound_value, "bound / xi", subLP)
end

function _compute_opt_cut_k_coefficients(
    subLP::ConnectorLP,
    kvals,
    bound_value::Real,
    params::GBCparam,
)
    coeffs = Dict{Any,Float64}()
    for a in subLP.A
        raw_value = Float64(kvals[a])
        if params.trim_coeff
            raw_value = min(raw_value, Float64(bound_value))
        end
        coeffs[a] = _round_cut_component(
            _sanitize_nonnegative_opt_cut_coefficient(raw_value, "k[$(a)]", subLP),
            params,
        )
    end
    return coeffs
end

function _compute_opt_cut_blc_g_coefficients(
    subLP::ConnectorLP,
    cutcoeff_BlC,
    gval,
    bound_value::Real,
    params::GBCparam,
)
    coeffs = Dict{Any,Float64}()
    adjusted_g = _adjust_g(subLP, gval)
    for a in keys(cutcoeff_BlC)
        raw_value = adjusted_g * Float64(cutcoeff_BlC[a])
        if params.trim_coeff
            raw_value = min(raw_value, Float64(bound_value))
        end
        coeffs[a] = _round_cut_component(
            _sanitize_nonnegative_opt_cut_coefficient(raw_value, "g[$(a)]", subLP),
            params,
        )
    end
    return coeffs
end

function _rounded_binary_y_values(subLP::ConnectorLP, y_vals)
    rounded = Dict{Any,Float64}()
    warned = false
    for a in subLP.A
        raw = Float64(y_vals[a])
        if abs(raw) <= 1e-8
            rounded[a] = 0.0
        elseif abs(raw - 1.0) <= 1e-8
            rounded[a] = 1.0
        else
            rounded[a] = min(max(Float64(ceil(raw)), 0.0), 1.0)
            if !warned
                @warn "ConnectorLP $(name(subLP.sub_solver)) received non-binary follower y-values when generating an optimality cut. JuBiC rounded them up to binary values for numerical stabilization."
                warned = true
            end
        end
    end
    return rounded
end

function _rounded_binary_x_values(subLP::ConnectorLP, x_vals, params::GBCparam)
    if !params.integer_obj
        return x_vals
    end

    rounded = Dict{Any,Float64}()
    for a in subLP.A
        raw = Float64(x_vals[a])
        rounded[a] = min(max(Float64(round(raw)), 0.0), 1.0)
    end
    return rounded
end



"""
    genBenders_cut!(subLP::ConnectorLP, link_vals::Dict{Tuple{Int64, Int64}, Float64}, params::GBCparam, time_limit)

Generate an general Benders (feasibility or optimality) cut.

# Arguments
- 'subLP::ConnectorLP{T}': The LP that must be solved for the sub_problem. Note that the LP is modified during the solving process.
- 'link_vals::Dict{T, Float64}': The mapping of resource to values of the linking variables from the master.
- 'params::GBCparam': Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns

- 'feas::Bool': True iff we need a feasibility cut.
- 'cut': The left-hand (i.e., non-trivial) side of the cut.
- 'bigMcut': If big M coefficients for BlC were computed as subroutine, return the right-hand (i.e., non-trivial) side of the BlC. Otherwise, return 'nothing'. 
- 'pobj': The contribution of the cut to the objective for the current solution (0 if feasibility cut).
"""
function genBenders_cut!(subLP::ConnectorLP{T}, link_vals::Dict{T,Float64}, params::GBCparam, time_limit) where T
    _reset_numeric_state!(subLP)
    # adjust sub_problem by setting new objective
    @debug "We now solve for the found optimal master solution the ConnectorLP $(name(subLP.sub_solver))."

    foundfeas, optL2, optL2_risk, y_vals = solve_sub_for_x(subLP.sub_solver, link_vals, params, time_limit)
    new_obj = subLP.lp[:s] - optL2 * subLP.lp[:g]
    new_obj +=
        -sum([
            round(capacity_linking(subLP.sub_solver, a, params) * link_vals[a]) * subLP.lp[:k][a]  # TODO: lets see if rounding helps numerics
            for a in subLP.A
        ])
    @objective(subLP.lp, Max, new_obj)

    # solve sub LP for new master
    @debug "Start iterative solution procedure for connector $(name(subLP.sub_solver))."

    # based if we found solution for current sub_problem, we need optimality or feasibility cut
    feas = !foundfeas  # a bit irreating, as we need feasibility cut iff sub is infeasible
    if feas
        # if infeasible, we do not have unbounded LP because of variable bounds added when initialized. This has some advantages:
        # 1. We can do the pareto-optimal cuts approach on this adjusted problem
        # 2. It seems to be numerical more stable for this LP that relying on dual information
        # 3. There seems to be only very little support for retrieving unbounded ray in JuMP

        # fix g to make sure that we obtain a ray with desired properties
        fix_g_constraint = @constraint(subLP.lp, subLP.lp[:g] == 0)

        @debug "Solving ConnectorLP $(name(subLP.sub_solver)) for feasibility cut generation."
        time_iterate = iterate_subsolver(subLP, params, time_limit)  # solve with new constraint

        
        # pareto optimality step for feasibility cuts
        if params.pareto == PARETO_OPTIMALITY_AND_FEASIBILITY
            time_limit_pareto = time_limit - time_iterate
            @debug "Start pareto-optimal Benders cut generation procedure for connector $(name(subLP.sub_solver)) for feasibility cut construction with remaining time limit $time_limit_pareto."
            pareto_optimal_decomposition(subLP, new_obj, optL2, params, time_limit_pareto)
        end

        # build feas cut
        cut = build_feas_cut(subLP)
        pobj = 0
        bigMcut = nothing
    else
        # solve sub_problem 
        @debug "Solving ConnectorLP $(name(subLP.sub_solver)) for optimality cut generation."
        time_iterate = iterate_subsolver(subLP, params, time_limit)
        time_limit_build_cut = time_limit - time_iterate
        pobj = value(new_obj)

        if get(subLP.numeric_state, :accepted_numerically, false)
            @debug "ConnectorLP $(name(subLP.sub_solver)) was accepted numerically after detecting a repeated cut. We skip adding the duplicate cut to the ConnectorLP itself and now build the usual optimality cut using the numerically stabilized g-value."
        elseif params.pareto == PARETO_OPTIMALITY_AND_FEASIBILITY || params.pareto == PARETO_OPTIMALITY_ONLY
            # pareto optimality step for optimality cuts
            time_limit_pareto = time_limit - time_iterate
            pareto_snapshot = _connector_solution_snapshot(subLP)
            try
                if !haskey(params.stats.data, "ConnectorLPTimePareto")
                    new_stat!(params.stats, "ConnectorLPTimePareto", 0.0)
                end
                @debug "Start pareto-optimal Benders cut generation procedure for connector $(name(subLP.sub_solver)) for optimality cut construction with remaining time limit $time_limit_pareto."
                pareto_time = @elapsed pareto_optimal_decomposition(subLP, new_obj, optL2, params, time_limit_pareto)
                add_stat!(params.stats, "ConnectorLPTimePareto", pareto_time)
                time_limit_build_cut = time_limit_pareto
                pobj = value(new_obj)
            catch err
                @warn "Pareto-optimal cut generation failed for connector $(name(subLP.sub_solver)). JuBiC falls back to the pre-pareto connector solution and continues with the standard cut. Error: $(sprint(showerror, err))"
                params.stats.data["Opt_status_override"] = "Numerics"
                params.stats.data["GBCStatus"] = "Numerics"
                subLP.numeric_state[:cut_solution_snapshot] = pareto_snapshot
            end
        end

        # build the usual optimality cut. In the numerical fallback path this uses
        # the stabilized g-value stored in subLP.numeric_state[:g_override].
        cut, bigMcut = build_opt_cut(subLP, optL2, optL2_risk, y_vals, link_vals, params, time_limit_build_cut)
    end

    # clean up and return
    if feas
        delete(subLP.lp, fix_g_constraint)
    end
    pareto_optimal_decomposition_cleanup(subLP)
    if !params.warmstart
        # okay, I understand that these tests are kind of interesting, BUT I never felt so stupid implementing a feature
        @debug "As requested cleaning up connectorLP $(name(subLP.sub_solver)) by removing all generated constraints."
        for cref in all_constraints(subLP.lp; include_variable_in_set_constraints=false)
            @debug "Now deleting constraint $cref from ConnectorLP $(name(subLP.sub_solver))."
            delete(subLP.lp, cref)
        end
        empty!(subLP.my_subsolutions)
    end
    @debug "Finished solving connectLP $(name(subLP.sub_solver))."
    return feas, cut, bigMcut, pobj
end


"""
    build_opt_cut(subLP::ConnectorLP, optL2, y_vals, x_vals, params::GBCparam, timelimit)

# Arguments
- 'subLP::ConnectorLP{T}': The LP that must be solved for the sub_problem. Note that the LP is modified during the solving process.
- 'optL2': The value of optimal L2 solution.
- 'optL2_risk': The value of the optimal L2 solution evaluated with L1 objective.
- 'y_vals': The values of the optimal L2 solution for the current master solution. 
- 'x_vals': The current master solution.
- 'params::GBCparam': Parameters passed down from the main solver.
- 'time_limit': The time limit for this subroutine. If exceeded, throws a 'TimeoutException'.

# Returns

- 'cut': The left-hand (i.e., non-trivial) side of the BlC constraint.
- 'bigMcut': If big M coefficients of BlC were computed, return here the right-hand (i.e., non trivial) side of BlC constraint. Otherwise, return 'nothing'
"""
function build_opt_cut(subLP::ConnectorLP, optL2, optL2_risk, y_vals, x_vals, params::GBCparam, timelimit)
    master_vars = subLP.link_vars
    bigMcut = nothing
    rounded_x_vals = _rounded_binary_x_values(subLP, x_vals, params)

    # g term of the cut
    gval = _snapshot_g(subLP)

    # s term of the cut
    sval = _snapshot_s(subLP)
    cut_rhs = _adjust_optcut_constant(sval - optL2 * gval, optL2_risk, subLP, params)
    cut = cut_rhs
    incumbent_cut_value = Float64(cut)
    incumbent_reference_value = Float64(cut_rhs)

    kvals = _snapshot_k(subLP)
    bound_value = _compute_opt_cut_bound_coefficient(subLP, sval, optL2, gval, kvals, x_vals)
    k_coeffs = _compute_opt_cut_k_coefficients(subLP, kvals, bound_value, params)

    cut -= sum(k_coeffs[a] * master_vars[a] for a in subLP.A)
    incumbent_cut_value -= sum(k_coeffs[a] * rounded_x_vals[a] for a in subLP.A)
    incumbent_reference_value -= sum(k_coeffs[a] * rounded_x_vals[a] for a in subLP.A)

    @debug "ConnectorLP $(name(subLP.sub_solver)) uses k-coefficients $(k_coeffs) and base g-coefficient xi=$(bound_value) for this optimality cut."

    ## generate cut from bound or by solving Lagrangian dual for BlC
    blc_subroutine = params.bigMwithLC && gval > 0 
    if blc_subroutine
        @debug "The conditions are met s.t. we generate a BlC for obtaining better big M coef. It is g=$gval "
        # If g=0, investing computational efford into generating a Lagrangian cut is just waste of computational ressources
        # solve subroutine approximating 
        synchronize_blc(subLP, subLP.blc_cut_generator)  # update ConnectorLP_BlC s.t. we preserve the subsolver solutions we found till now
        blc_cut, _, cutcoeff_BlC = genBenderslike_cut!(subLP.blc_cut_generator, x_vals, params, timelimit)  # if we get a timeout error, we just let it through
        
        # build GBC cut, i.e., Bilevel Lagrangian cut, coefficients for theta terms
        @debug "The big M values computed from Lagrangian dual now used in BlC are: $cutcoeff_BlC and optL2 * gval=$(ceil(Int, optL2 * gval))"
        blc_g_coeffs = _compute_opt_cut_blc_g_coefficients(subLP, cutcoeff_BlC, gval, bound_value, params)
        for a in keys(blc_g_coeffs) 
            rounded_coef = blc_g_coeffs[a]
            # we do not multiply the big M with yvals as there can exist multiple equivalent solutions of L2 with same objective, leeding to different BlC 
            cut -= rounded_coef * (1-master_vars[a]) #* y_vals[a] # TODO: rounding to avoid numeric trouble
            incumbent_cut_value -= rounded_coef * (1 - rounded_x_vals[a])
        end
        add_stat!(params.stats, "NBigMlagCuts", 1)

        # Overlap between k and BlC-based big-M coefficients can occur on
        # numerically delicate instances. We continue, but flag the run.
        for a in subLP.A
            blc_coef = get(cutcoeff_BlC, a, 0.0)
            if kvals[a] > 0 && blc_coef > 0
                _warn_cut_overlap!(subLP, params, a, "blc_coef", blc_coef, kvals[a])
            end
        end

        # generate bigMcut
        bigMcut = blc_cut
    else
        if gval > 0
            rounded_y_vals = _rounded_binary_y_values(subLP, y_vals)
            theta_xXg = sum(bound_value * (1 - master_vars[a]) * rounded_y_vals[a] for a in subLP.A)
            cut -= theta_xXg
            incumbent_cut_value -= sum(bound_value * (1 - rounded_x_vals[a]) * rounded_y_vals[a] for a in subLP.A)

            # Overlap between k and the y-pattern can occur on numerically
            # delicate instances. We continue, but flag the run.
            for a in subLP.A
                if kvals[a] > 0 && rounded_y_vals[a] > 0
                    _warn_cut_overlap!(subLP, params, a, "y", rounded_y_vals[a], kvals[a])
                end
            end
        else
            @debug "ConnectorLP $(name(subLP.sub_solver)) omits g-based optimality-cut coefficients because g=0."
        end
    end

    _record_opt_cut_validation!(subLP, params, incumbent_cut_value, incumbent_reference_value, cut_rhs)

    
    # update gbc solver with cuts generated by connectorLP_BlC. We need to do it at the end as we otherwise would need to reoptimize our model
    if blc_subroutine
        synchronize_gbc(subLP.blc_cut_generator, subLP)  # update this ConnectorLP s.t. we preserve the subsolver solutions we found till now

        # if debug mode, print again adjusted LP to file
        if params.debbug_out
            write_to_file(
                subLP.lp,
                "$(params.output_folder_path)/conLP_$(name(subLP.sub_solver)).$(params.file_format_output)",
            )
        end
    end
    _reset_numeric_state!(subLP)
    return cut, bigMcut
end

function build_feas_cut(subLP::ConnectorLP)
    master_vars = subLP.link_vars
    sval = _snapshot_s(subLP)
    kvals = _snapshot_k(subLP)
    fcut = sum(
        (kvals[a] / sval > 0.9 ? 1.0 : kvals[a] / sval) * master_vars[a]
        for a in subLP.A
    )
    return fcut
end

"""
    check_solution_status_LP(lp::JuMP.Model)

Check if our connector LP model terminated with an optimal solution
"""
function check_solution_status_LP(me::ConnectorLP)
    status = termination_status(me.lp) 
    if status == MOI.TIME_LIMIT
        @debug "We hit the time limit while solving the ConnectorLP of user $(name(me))"
        throw(TimeoutException( "We hit the time limit while ConnectorLP of user $(name(me))" ))
    end


    if status == MOI.INFEASIBLE
        @error "The ConnectorLP $(name(me)) is infeasible. Starting computations of IIS. But most of the time this is implied by numerical issues."
        compute_conflict!(me.lp)
        @debug "IIS computed, now printing it"
        iis_model, _ = copy_conflict(me.lp)
        if get_attribute(me.lp, MOI.ConflictStatus()) == MOI.CONFLICT_FOUND
            @error iis_model # printing to file just causes error...
        else
            @error "The IIS was not computed succesfully??? Most likely numerics??"
        end
        error("ConnectorLP $(name(me)) was infeasible. Computed IIS but stopping solution process (as it is clearly a bug). Most likely, it was caused by numerical issues.")
    elseif status == MOI.DUAL_INFEASIBLE
        error("The ConnectorLP $(name(me)) was unbounded what violates the way we handle its. Because we add significantly large bounds, we avoid unboundnes; so, this should not have happended.")
    end

    if !(status == MOI.OPTIMAL)
        error("ConnectorLP $(name(me)) could not be solver to optimality but terminated with status $(status). Connot continue as this is undefined behavior.")
    end
end

function _duplicate_connector_cut_error(
    subLP::ConnectorLP,
    csc::ConSubsolCut,
    sval::Real,
    opt_obj::Real,
    tolerance::Real,
)
    target = opt_obj + tolerance
    resource_text = isempty(csc.res) ? "no resources" : join(string.(csc.res), ", ")
    return NumericalIssueException(
        "JuBiC detected a possible numerical issue while solving ConnectorLP $(name(subLP.sub_solver)). " *
        "The same cut for resource pattern [$(resource_text)] was generated again even though it was already added earlier. " *
        "The current connector-LP value is s=$(sval), while the follower still certifies a value of $(opt_obj). " *
        "To avoid cycling on numerically unstable cuts, JuBiC stops here. " *
        "For this cut to be considered non-violated, the connector-LP would need to bring s below approximately $(target) (using tolerance $(tolerance)).",
        "NumericalIssue_DuplicateCut",
    )
end


"""
    iterate_subsolver(subLP::ConnectorLP, params::GBCparam)

Solve the LP by a separation procedure. 
Uses an iterative while loop in the implementation, avoiding stack overflow exceptions compared to this function recursive version.

In case the separation takes longer that the set time limit, throw an exception. 

# Return
    The time it required to execute this function
"""
function iterate_subsolver(subLP::ConnectorLP, params::GBCparam, time_limit)
    # this part of the solver can run into quite some nasty infinity loops. To prevent the software just freezing (or seeming to freeze for the user),
    ## we stop the separation process in case we reach the timelimit set in the parameters (what still can take long, but at least the user is expected to wait this long) 
    current_time = time()

    if !haskey(params.stats.data, "ConnectorLPIterations")
        new_stat!(params.stats, "ConnectorLPIterations", 0)
    end
    if !haskey(params.stats.data, "ConnectorLPTimeLP")
        new_stat!(params.stats, "ConnectorLPTimeLP", 0.0)
    end
    if !haskey(params.stats.data, "ConnectorLPTimePricing")
        new_stat!(params.stats, "ConnectorLPTimePricing", 0.0)
    end

    violated_cut = true # true as long as violated constraint could exist in LP
    while violated_cut
        add_stat!(params.stats, "ConnectorLPIterations", 1)
        # solve the sub_problem iteratively (but first debbug output)
        if params.debbug_out
            write_to_file(
                subLP.lp,
                "$(params.output_folder_path)/conLP_$(name(subLP.sub_solver)).$(params.file_format_output)",
            )
        end

        # first, solve the LP
        #@debug subLP.lp # this debug output is not helpfull
        lp_time = @elapsed optimize!(subLP.lp) # Assumption: Solving the LP consumes neglectable time
        add_stat!(params.stats, "ConnectorLPTimeLP", lp_time)
        check_solution_status_LP(subLP)  

        # solve sub_problem for found solution
        kvals = Dict(a => value(subLP.lp[:k][a]) for a in subLP.A)
        @debug "The found sub_problem ConnectorLP solution is s=$(value(subLP.lp[:s])), g=$(value(subLP.lp[:g])), and non-zero k=$(Dict(key => k for (key, k) in kvals if k != 0)). "
        pricing_time = @elapsed begin
            sub_solver = separation!(
                subLP.sub_solver,
                value(subLP.lp[:s]),
                value(subLP.lp[:g]),
                kvals,
                params,
                time_limit
            )
            subLP.numeric_state[:last_sub_solver_solution] = sub_solver
        end
        add_stat!(params.stats, "ConnectorLPTimePricing", pricing_time)
        sub_solver = subLP.numeric_state[:last_sub_solver_solution]
        delete!(subLP.numeric_state, :last_sub_solver_solution)

        # if we found a violated constraint, add it and resolve
        if sub_solver.vio
            # add constraints
            @debug "For connector $(name(subLP.sub_solver)), we found a sub_problem solution that violates current LP solution and uses resources $(sub_solver.A_sub)."
            csc = ConSubsolCut(sub_solver.A_sub, sub_solver.obj_second_level, sub_solver.obj_first_level)
            if csc in subLP.my_subsolutions
                sval = value(subLP.lp[:s])
                opt_obj = sub_solver.obj_compare
                g_current = value(subLP.lp[:g])
                delta_g_connector, delta_g_sub, delta_g_total = _estimate_g_rounding_tolerance(subLP, sub_solver)
                g_required = _required_g_value(subLP, sub_solver)
                g_rounded = g_current + delta_g_total

                if g_rounded >= g_required
                    warning_msg = _duplicate_cut_warning_message(
                        subLP,
                        csc,
                        sval,
                        opt_obj,
                        g_current,
                        g_required,
                        g_rounded,
                        delta_g_connector,
                        delta_g_sub,
                        delta_g_total,
                    )
                    @warn warning_msg
                    subLP.numeric_state[:accepted_numerically] = true
                    subLP.numeric_state[:g_override] = g_rounded
                    params.stats.data["Opt_status_override"] = "Numerics"
                    params.stats.data["GBCStatus"] = "Numerics"
                    violated_cut = false
                    continue
                else
                    tolerance = 10e-4
                    err = _duplicate_connector_cut_error(subLP, csc, sval, opt_obj, tolerance)
                    @error err.message
                    throw(err)
                end
            end
            new_const_left =
                subLP.lp[:s] - sum(subLP.lp[:k][a] for a in sub_solver.A_sub; init=0) -
                sub_solver.obj_second_level * subLP.lp[:g] 
            c = @constraint(subLP.lp, new_const_left <= sub_solver.obj_first_level)
            @debug "We added an violated constraint $(c) for connector $(name(subLP.sub_solver)). Continue separation."

            # save found solution to our internal storage
            push!(subLP.my_subsolutions, csc)

            # check for time limit
            if current_time + params.runtime < time() 
                @error "We reached the time limit when solving ConnectorLP $(name(subLP.sub_solver)). Terminating cut generation process."
                throw(TimeoutException("We reached the time limit when solving ConnectorLP $(name(subLP.sub_solver)). Terminating GBC solution procedure."))
            end

            # resolve
            violated_cut = true
        else
            # no violated constraint many more
            violated_cut = false
        end
    end

    return time() - current_time # return how long the process took
end


"""
    iterate_subsolver_recursive(subLP::ConnectorLP, params::GBCparam)

@Deprecated !

    Solves the current version of the LP using a recursive call of this function. 
Compared to the non-recursive 'iterate_subsolver' version of this function, it should be in my understanding of Julia.
However, it can easily lead to an 'StackOverflowError' if multiple constraints need to be added within one iteration. 
Therefore, this function DEPRECATED and should be avoided. 
"""
function iterate_subsolver_recursive(subLP::ConnectorLP, params::GBCparam)
    # solve the sub_problem iteratively (but first debbug output)
    if params.debbug_out
        write_to_file(
            subLP.lp,
            "$(params.output_folder_path)/conLP_$(name(subLP.sub_solver)).$(params.file_format_output)",
        )
    end

    # first, solve the LP
    #@debug subLP.lp # this debug output is not helpfull
    optimize!(subLP.lp)
    check_solution_status_LP(subLP)

    # solve sub_problem for found solution
    kvals = Dict(a => value(subLP.lp[:k][a]) for a in subLP.A)
    @debug "The found sub_problem ConnectorLP solution is s=$(value(subLP.lp[:s])), g=$(value(subLP.lp[:g])), and k=$(kvals)"
    sub_solver = separation!(
        subLP.sub_solver,
        value(subLP.lp[:s]),
        value(subLP.lp[:g]),
        kvals,
        params
    )

    # if we found a violated constraint, add it and resolve
    if sub_solver.vio
        # add constraints
        @debug "For connector $(name(subLP.sub_solver)), we found a sub_problem solution that violates current LP solution and uses resources $(sub_solver.A_sub)."
        new_const_left =
            subLP.lp[:s] - sum(subLP.lp[:k][a] for a in sub_solver.A_sub; init=0) -
            sub_solver.obj_second_level * subLP.lp[:g]
        c = @constraint(subLP.lp, new_const_left <= sub_solver.obj_first_level)
        @debug "We added an violated constraint $(c) for connector $(name(subLP.sub_solver)). Continue separation."

        # resolve
        iterate_subsolver_recursive(subLP, params)
    end
end



"""
    pareto_optimal_decomposition(subLP::ConnectorLP, lp_obj, g_obj_coef, params::GBCparam)

When called after LP was solved to optimality, transforms the LP to generate a pareto-optimal cut and then resolves it.

# Arguments
- 'subLP::ConnectorLP': The LP
- 'lp_obj': The LP objective function (as JuMP expression)
- 'g_obj_coef': The coefficient of g in the objective function, i.e., the current solution value of the follower problem for selected master solution, or 0 if infeasible.
- 'params::GBCparam': The parameters
- 'time_limit': The runtime limit for this subroutine. In case of time out, throw an 'TimeoutException'
"""
function pareto_optimal_decomposition(subLP::ConnectorLP, lp_obj, g_obj_coef, params::GBCparam, time_limit)
    # resolve the subLP with steps necessary for pareto-optimal cuts
    cssol = value(subLP.lp[:s])
    cgsol = value(subLP.lp[:g])  # TODO: no rounding of g value for pareto optimal cut good idea 

    # The 0.1 term -> you still want to minimize the g
    bigMterm = (cssol - g_obj_coef * cgsol - subLP.lower_bound_obj_contribution) + 0.1 # this is just an educated guess, i.e., heuristic pareto-optimal cut generation
    if bigMterm < 0
        error(
            "We got a negative big_m in ConnectorLP $(name(subLP)) because s=$(cssol), g=$(cgsol), gcoef=$(g_obj_coef), and bound=$(subLP.lower_bound_obj_contribution).",
        )
    end

    # build adjusted model
    current_obj = objective_value(subLP.lp)
    pBd_obj = sum(subLP.lp[:k]) + bigMterm * subLP.lp[:g]
    pBd_tol = _pareto_band_tolerance(subLP, current_obj, g_obj_coef, params)
    @objective(subLP.lp, Min, pBd_obj)
    @constraint(subLP.lp, pBd_const_lb, lp_obj >= current_obj - pBd_tol)
    @constraint(subLP.lp, pBd_const_ub, lp_obj <= current_obj + pBd_tol)
    @debug "Pareto band for connector $(name(subLP.sub_solver)) fixes the original objective within +/-$(pBd_tol) around $(current_obj)."

    # print pareto adjusted model to file in debbug mode 
    if params.debbug_out
        write_to_file(
            subLP.lp,
            "$(params.output_folder_path)/conLP_ParetoInit_$(name(subLP.sub_solver)).$(params.file_format_output)",
        )
    end

    # resolve optimization problem
    iterate_subsolver(subLP, params, time_limit)
end


function pareto_optimal_decomposition_cleanup(subLP::ConnectorLP)
    # do the cleanup required for pareto-optimal Benders cuts
    if haskey(object_dictionary(subLP.lp), :pBd_const_lb)
        delete(subLP.lp, subLP.lp[:pBd_const_lb])
        unregister(subLP.lp, :pBd_const_lb)
    end
    if haskey(object_dictionary(subLP.lp), :pBd_const_ub)
        delete(subLP.lp, subLP.lp[:pBd_const_ub])
        unregister(subLP.lp, :pBd_const_ub)
    end
    if haskey(object_dictionary(subLP.lp), :pBd_const)
        delete(subLP.lp, subLP.lp[:pBd_const])
        unregister(subLP.lp, :pBd_const)
    end
    # no need to set new obj. This will happend when new master solution passed for next run
end

"""
    synchronize_blc(con::ConnectorLP{T}, con_blc::ConnectorLP_BlC{T}) where T

Synchronize the subsolver solutions and constraints from ConnectorLP to its BlC counterpart. 
"""
function synchronize_blc(con::ConnectorLP{T}, con_blc::ConnectorLP_BlC{T}) where T
    return #TODO: We solve very different Lagrangian subproblems. Transferring one solution to another is not straightforward. However, we leave the function to give room a possible extension in the future.
    """
    for scs in con.my_subsolutions
        if !(scs in con_blc.my_subsolutions_blc)
            # add solution to list 
            #push!(con_blc.my_subsolutions_blc, scs)

            add_cut_from_scs_blc(con_blc, scs; synchro=true)
            @debug "Added solution $scs and corresponding constraint from ConnectorLP to its BlC counterpart."
        end
    end
    """
end


"""
    synchronize_gbc(con_blc::ConnectorLP_BlC{T}, con::ConnectorLP{T}) where T
    
Synchronize the subsolver solutions and constraints from ConnectorLP to its BlC counterpart. 
"""
function synchronize_gbc(con_blc::ConnectorLP_BlC{T}, con::ConnectorLP{T}) where T 
    return # TODO: We solve very different Lagrangian subproblems. Transferring one solution to another is not straightforward. However, we leave the function to give room a possible extension in the future.
    """
    throw(ErrorException("Does not work after rework of BlC"))
    for scs in con_blc.my_subsolutions_blc
        if !(scs in con.my_subsolutions)
            # add solution to list 
            push!(con.my_subsolutions, scs)

            # add constraint
            new_const_left =
                con.lp[:s] - sum(con.lp[:k][a] for a in scs.res; init=0) - scs.objL2 * con.lp[:g]
            @constraint(con.lp, new_const_left <= scs.objL1, base_name="synchro")
            
            @debug "Added solution $scs and corresponding constraint from ConnectorLP_BlC to its GBC counterpart."
        end
    end
    """
end


"""
    _adjust_g(subLP::ConnectorLP, g_value)

Round value of g to avoid numeric trouble. Only intended use case is to round 'g_value' before multiplying it with other coef of x-variables from master. 
This should ensure higher numric stability while preserving correctness of cut.
"""
function _adjust_g(subLP::ConnectorLP, g_value)
    return ceil(g_value, digits=subLP.g_round_digits) 
end


"""
    _adjust_optcut_constant(constvalue, optL2_risk, subLP::ConnectorLP, params::GBCparam)

Numerically savely adjust the constant term of the GBC optimality cut.  
"""
function _adjust_optcut_constant(constvalue, optL2_risk, subLP::ConnectorLP, params::GBCparam)
    num_tol = 10e-2  # TODO_: Again, hard coded numerics
    if params.integer_obj
        rounded_constvalue = round(Float64(constvalue))
        rounded_optL2_risk = round(Float64(optL2_risk))
        if rounded_constvalue == rounded_optL2_risk
            @debug "Set constant term of GBC to rounded integer value $(rounded_constvalue)"
            return rounded_constvalue
        end
        @debug "In GBCSolver $(name(subLP)), the next cut to be generated has constant term $(constvalue) but the current solution should have value $(optL2_risk).
            After rounding, the constant term is $(rounded_constvalue) while the supplied current solution value rounds to $(rounded_optL2_risk).
            This can occur if there are multiple solutions with the same L2 objective, but the subsolver chooses one with a non-minimal objective value according to the L1 objective.
            Since we cannot guess the correct value for the cut constant, we keep the rounded constant term."
        return rounded_constvalue
    end

    target_rhs = optL2_risk
    if !(constvalue >= target_rhs - num_tol && constvalue <= target_rhs + num_tol)
        @debug "In GBCSolver $(name(subLP)), the next cut to be generated has constant term $(constvalue) but the current solution should have value $(optL2_risk).
            This can occur if there are multiple solutions with the same L2 objective, but the subsolver chooses one with a non-minimal objective value according to the L1 objective.
            Since we cannot guess the correct value for the cut constant, we omit rounding here."
        return constvalue
    end
    @debug "Set constant term of GBC to $(target_rhs)"
    return target_rhs
end


Base.show(io::IO, conLP::ConnectorLP) = print(io, "ConnectorLP_$(name(conLP.sub_solver))")
