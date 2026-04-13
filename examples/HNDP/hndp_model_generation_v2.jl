using JuBiC
using JuMP
using Graphs

include("hndp_network_generation_v2.jl")

const HNDP_BIGM_FIXED_NETWORK_PATH = :fixed_network_path
const HNDP_BIGM_N_MINUS_ONE = :n_minus_one_most_expensive

"""
    build_hndp_blc_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH)

Build a JuBiC `BlCMaster` instance for the passed HNDP network.
"""
function build_hndp_blc_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
)
    all_arcs = _hndp_all_arcs(hndp)
    user_names = [string(user.uname) for user in hndp.users]

    hpr = Model(() -> get_next_optimizer(solver))
    @variable(hpr, x[all_arcs], Bin)
    for arc in all_arcs
        if !(arc in hndp.edgeA)
            @constraint(hpr, x[arc] == 1, base_name = "fix_$(arc)")
        end
    end

    master_terms = Dict{String,Any}()
    sub_terms = Dict{String,Any}()
    for user in hndp.users
        leader_obj, follower_obj, _ = _add_user_flow_constraints!(
            hpr,
            user,
            hndp,
            all_arcs;
            xvars=x,
            binary_flow=true,
            add_strong_duality=false,
        )
        user_name = string(user.uname)
        opt_l1 = @variable(hpr, base_name = "optL1_$(user_name)")
        opt_l2 = @variable(hpr, base_name = "optL2_$(user_name)")
        @constraint(hpr, opt_l1 == leader_obj, base_name = "leader_obj_$(user_name)")
        @constraint(hpr, opt_l2 == follower_obj, base_name = "follower_obj_$(user_name)")
        master_terms[user_name] = opt_l1
        sub_terms[user_name] = opt_l2
    end

    constructioncost = sum(hndp.edge_price[a] * x[a] for a in hndp.edgeA)
    @variable(hpr, construction_cost_var)
    @constraint(hpr, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(hpr, Min, construction_cost_var + sum(values(master_terms)))

    big_m_by_user = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
    big_m_function(a, uname) = big_m_by_user[string(uname)]
    xdec = Dict(a => x[a] for a in hndp.edgeA)
    master = BlCMaster(hpr, hndp.edgeA, xdec, big_m_function, user_names, sub_terms)

    subs = Any[]
    for user in hndp.users
        sub_model = Model(() -> get_next_optimizer(solver))
        leader_obj, follower_obj, flow = _add_user_flow_constraints!(
            sub_model,
            user,
            hndp,
            all_arcs;
            xvars=nothing,
            binary_flow=true,
            add_strong_duality=false,
        )
        @objective(sub_model, Min, follower_obj)
        set_silent(sub_model)
        y_vars = Dict(a => flow[a] for a in hndp.edgeA)
        push!(
            subs,
            SubSolverJuMP(
                string(user.uname),
                sub_model,
                hndp.edgeA,
                y_vars,
                leader_obj,
                follower_obj,
            ),
        )
    end

    return Instance(master, subs)
end

"""
    build_hndp_sd_instance(hndp, solver; big_m_mode=HNDP_BIGM_FIXED_NETWORK_PATH, indicator_constraints=false, bound_duals=true)

Build a strong-duality reformulation of the HNDP. This formulation is only
valid for instances without weight bounds, because then each follower problem is
an LP shortest-path model.
"""
function build_hndp_sd_instance(
    hndp::HNDPwC,
    solver::SolverWrapper;
    big_m_mode::Symbol=HNDP_BIGM_FIXED_NETWORK_PATH,
    indicator_constraints::Bool=false,
    bound_duals::Bool=true,
)
    any(!isnothing(user.weighlimit) for user in hndp.users) &&
        throw(ArgumentError("The strong-duality HNDP formulation is currently only supported for instances without weight bounds."))

    all_arcs = _hndp_all_arcs(hndp)
    mip = Model(() -> get_next_optimizer(solver))
    @variable(mip, x[all_arcs], Bin)
    for arc in all_arcs
        if !(arc in hndp.edgeA)
            @constraint(mip, x[arc] == 1, base_name = "fix_$(arc)")
        end
    end

    big_m_by_user = _derive_user_big_m_values(hndp, all_arcs, solver, big_m_mode)
    big_m_function(a, uname) = big_m_by_user[string(uname)]

    leader_terms = Dict{String,Any}()
    for user in hndp.users
        leader_obj, follower_obj, _ = _add_user_flow_constraints!(
            mip,
            user,
            hndp,
            all_arcs;
            xvars=x,
            binary_flow=false,
            add_strong_duality=true,
            big_m_function=big_m_function,
            indicator_constraints=indicator_constraints,
            bound_duals=bound_duals,
        )
        user_name = string(user.uname)
        opt_l2 = @variable(mip, base_name = "optL2_$(user_name)")
        @constraint(mip, opt_l2 == follower_obj, base_name = "follower_obj_$(user_name)")
        leader_terms[user_name] = leader_obj
    end

    constructioncost = sum(hndp.edge_price[a] * x[a] for a in hndp.edgeA)
    @variable(mip, construction_cost_var)
    @constraint(mip, construction_cost_var == constructioncost, base_name = "construction_cost")
    @objective(mip, Min, construction_cost_var + sum(values(leader_terms)))

    return Instance(MIPMaster(mip), nothing)
end

function _hndp_all_arcs(hndp::HNDPwC)
    return [(src(e), dst(e)) for e in edges(hndp.mygraph)]
end

function _add_user_flow_constraints!(
    model::JuMP.Model,
    user::User,
    hndp::HNDPwC,
    all_arcs;
    xvars,
    binary_flow::Bool,
    add_strong_duality::Bool,
    big_m_function=nothing,
    indicator_constraints::Bool=false,
    bound_duals::Bool=true,
)
    if binary_flow
        flow = @variable(model, [a in all_arcs], Bin, base_name = "f_$(user.uname)")
    else
        flow = @variable(model, [a in all_arcs], lower_bound = 0, upper_bound = 1, base_name = "f_$(user.uname)")
    end

    if !isnothing(user.weighlimit)
        @constraint(
            model,
            sum(user.mweight[a...] * flow[a] for a in all_arcs) <= user.weighlimit,
            base_name = "weight_$(user.uname)",
        )
    end

    for node in vertices(hndp.mygraph)
        rhs = node == user.origin ? 1 : (node == user.destination ? -1 : 0)
        outgoing = sum((flow[(node, neigh)] for neigh in outneighbors(hndp.mygraph, node)); init=0.0)
        incoming = sum((flow[(neigh, node)] for neigh in inneighbors(hndp.mygraph, node)); init=0.0)
        @constraint(model, outgoing - incoming == rhs, base_name = "flow_$(user.uname)_$(node)")
    end

    if !isnothing(xvars)
        for a in all_arcs
            @constraint(model, flow[a] <= xvars[a], base_name = "link_$(user.uname)_$(a)")
        end
    end

    leader_obj = @expression(model, sum(user.mrisk[a...] * flow[a] for a in all_arcs))
    follower_obj = @expression(model, sum(user.mcost[a...] * flow[a] for a in all_arcs))

    if add_strong_duality
        isnothing(big_m_function) && throw(ArgumentError("A big-M function is required when adding strong-duality constraints."))
        isnothing(xvars) && throw(ArgumentError("Strong-duality constraints require first-level x-variables."))

        potentials = @variable(model, [node in vertices(hndp.mygraph)], lower_bound = 0, base_name = "pi_$(user.uname)")
        if bound_duals
            max_big_m = maximum(big_m_function(a, string(user.uname)) for a in all_arcs)
            for node in vertices(hndp.mygraph)
                set_upper_bound(potentials[node], max_big_m)
            end
        end

        for a in all_arcs
            if indicator_constraints
                @constraint(
                    model,
                    xvars[a] --> {potentials[a[1]] - potentials[a[2]] <= user.mcost[a...]},
                    base_name = "dual_ind_$(user.uname)_$(a)",
                )
            else
                @constraint(
                    model,
                    potentials[a[1]] - potentials[a[2]] <= user.mcost[a...] + big_m_function(a, string(user.uname)) * (1 - xvars[a]),
                    base_name = "dual_$(user.uname)_$(a)",
                )
            end
        end

        @constraint(
            model,
            follower_obj == potentials[user.origin] - potentials[user.destination],
            base_name = "strong_duality_$(user.uname)",
        )
    end

    return leader_obj, follower_obj, flow
end

function _derive_user_big_m_values(hndp::HNDPwC, all_arcs, solver::SolverWrapper, big_m_mode::Symbol)
    big_m_mode in (HNDP_BIGM_FIXED_NETWORK_PATH, HNDP_BIGM_N_MINUS_ONE) ||
        throw(ArgumentError("Unknown HNDP big-M mode $(big_m_mode)."))

    values = Dict{String,Float64}()
    for user in hndp.users
        uname = string(user.uname)
        if big_m_mode == HNDP_BIGM_FIXED_NETWORK_PATH
            path_cost = _fixed_network_path_cost(user, hndp, all_arcs, solver)
            isnothing(path_cost) &&
                throw(ArgumentError("No feasible path exists in the fixed network for user $(uname), so the fixed-network big-M mode is not available."))
            values[uname] = path_cost + 1.0
        else
            costs = sort([user.mcost[a...] for a in all_arcs], rev=true)
            n_take = min(length(costs), max(nv(hndp.mygraph) - 1, 1))
            values[uname] = sum(costs[1:n_take])
        end
    end
    return values
end

function _fixed_network_path_cost(user::User, hndp::HNDPwC, all_arcs, solver::SolverWrapper)
    model = Model(() -> get_next_optimizer(solver))
    @variable(model, flow[a in all_arcs], Bin)

    for a in hndp.edgeA
        @constraint(model, flow[a] == 0, base_name = "forbid_$(user.uname)_$(a)")
    end

    if !isnothing(user.weighlimit)
        @constraint(
            model,
            sum(user.mweight[a...] * flow[a] for a in all_arcs) <= user.weighlimit,
            base_name = "weight_$(user.uname)",
        )
    end

    for node in vertices(hndp.mygraph)
        rhs = node == user.origin ? 1 : (node == user.destination ? -1 : 0)
        outgoing = sum((flow[(node, neigh)] for neigh in outneighbors(hndp.mygraph, node)); init=0.0)
        incoming = sum((flow[(neigh, node)] for neigh in inneighbors(hndp.mygraph, node)); init=0.0)
        @constraint(model, outgoing - incoming == rhs, base_name = "flow_$(user.uname)_$(node)")
    end

    @objective(model, Min, sum(user.mcost[a...] * flow[a] for a in all_arcs))
    set_silent(model)
    optimize!(model)

    status = termination_status(model)
    if status == MOI.OPTIMAL || status == MOI.LOCALLY_SOLVED
        return objective_value(model)
    elseif status == MOI.INFEASIBLE
        return nothing
    end
    throw(ArgumentError("Unexpected status $(status) while computing fixed-network big-M for user $(user.uname)."))
end
