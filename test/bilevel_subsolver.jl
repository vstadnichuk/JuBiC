using JuBiC, JuMP, Gurobi, Test
import MathOptInterface as MOI

function _simple_follower_oracle_factory()
    return function (x_vals, params, time_limit)
        oracle = Model(optimizer)
        set_silent(oracle)
        set_time_limit_sec(oracle, time_limit)
        @variable(oracle, y[1:2], Bin)
        @constraint(oracle, y[1] <= x_vals[1])
        @constraint(oracle, y[2] <= x_vals[2])
        @constraint(oracle, y[1] == 1)
        @objective(oracle, Min, -y[1] - y[2])
        optimize!(oracle)

        status = termination_status(oracle)
        if status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
            return false, 0, 0, Dict(1 => 0.0, 2 => 0.0)
        elseif status != MOI.OPTIMAL
            error("Oracle lower-level solve terminated with status $status")
        end

        y_vals = Dict(1 => value(y[1]), 2 => value(y[2]))
        return true, objective_value(oracle), 10 * y_vals[2], y_vals
    end
end

function generate_gbc_simple_bilevel_instance_blc_subsolver()
    A = [1, 2]
    nsub = "SubBlCJuMP"

    mm = Model(optimizer)
    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1] - x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])

    sub = Model(optimizer)
    set_silent(sub)
    @variable(sub, y[1:2], Bin)
    @constraint(sub, cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    @objective(sub, Min, sub_obj)
    master_sub_obj = 10 * y[2]
    subS = SubSolverBlCJuMP(
        nsub,
        sub,
        A,
        y,
        master_sub_obj,
        sub_obj,
        _simple_follower_oracle_factory(),
        a -> 100.0,
    )

    return Instance(master, [subS])
end

function generate_gbc_simple_bilevel_instance_blc_subsolver_auto_oracle()
    A = [1, 2]
    nsub = "SubBlCJuMPAuto"

    mm = Model(optimizer)
    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1] - x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])

    sub = Model(optimizer)
    set_silent(sub)
    @variable(sub, y[1:2], Bin)
    @constraint(sub, cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    @objective(sub, Min, sub_obj)
    master_sub_obj = 10 * y[2]
    subS = SubSolverBlCJuMP(
        nsub,
        sub,
        A,
        y,
        master_sub_obj,
        sub_obj,
        a -> 100.0,
    )

    return Instance(master, [subS])
end

function generate_blclag_simple_bilevel_instance_blc_subsolver()
    A = [1, 2]
    nsub = "SubBlCJuMP"

    hpr = Model(optimizer)
    @variable(hpr, x[A], Bin)
    @variable(hpr, yh[A], Bin)
    @constraint(hpr, hpr_link_1, yh[1] <= x[1])
    @constraint(hpr, hpr_link_2, yh[2] <= x[2])
    @constraint(hpr, hpr_fix, yh[1] == 1)
    @objective(hpr, Min, x[1] - x[2] + 10 * yh[2])
    xdict = Dict(a => x[a] for a in A)
    sub_obj_hpr = @expression(hpr, -sum(yh[A]))
    master = BlCLagMaster(hpr, A, xdict, [nsub], Dict(nsub => sub_obj_hpr))

    sub = Model(optimizer)
    set_silent(sub)
    @variable(sub, y[1:2], Bin)
    @constraint(sub, cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    @objective(sub, Min, sub_obj)
    master_sub_obj = 10 * y[2]
    subS = SubSolverBlCJuMP(
        nsub,
        sub,
        A,
        y,
        master_sub_obj,
        sub_obj,
        _simple_follower_oracle_factory(),
        a -> 100.0,
    )

    return Instance(master, [subS])
end

function generate_blclag_simple_bilevel_instance_blc_subsolver_auto_oracle()
    A = [1, 2]
    nsub = "SubBlCJuMPAuto"

    hpr = Model(optimizer)
    @variable(hpr, x[A], Bin)
    @variable(hpr, yh[A], Bin)
    @constraint(hpr, hpr_link_1, yh[1] <= x[1])
    @constraint(hpr, hpr_link_2, yh[2] <= x[2])
    @constraint(hpr, hpr_fix, yh[1] == 1)
    @objective(hpr, Min, x[1] - x[2] + 10 * yh[2])
    xdict = Dict(a => x[a] for a in A)
    sub_obj_hpr = @expression(hpr, -sum(yh[A]))
    master = BlCLagMaster(hpr, A, xdict, [nsub], Dict(nsub => sub_obj_hpr))

    sub = Model(optimizer)
    set_silent(sub)
    @variable(sub, y[1:2], Bin)
    @constraint(sub, cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    @objective(sub, Min, sub_obj)
    master_sub_obj = 10 * y[2]
    subS = SubSolverBlCJuMP(
        nsub,
        sub,
        A,
        y,
        master_sub_obj,
        sub_obj,
        a -> 100.0,
    )

    return Instance(master, [subS])
end

function test_gbc_simple_bilevel_with_blc_subsolver()
    model = generate_gbc_simple_bilevel_instance_blc_subsolver()
    outdir = logging_folder * "/gbc_simple_bilevel_blc_subsolver"
    mkpath(outdir)
    parameter = GBCparam(
        GurobiSolver(),
        true,
        outdir,
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    stats = solve_instance!(model, parameter)
    @test haskey(stats.data, "Opt")
    @test stats.data["Opt"] ≈ 1
end

function test_blclag_simple_bilevel_with_blc_subsolver()
    model = generate_blclag_simple_bilevel_instance_blc_subsolver()
    outdir = logging_folder * "/blclag_simple_bilevel_blc_subsolver"
    mkpath(outdir)
    parameter = BlCLagparam(
        GurobiSolver(),
        true,
        outdir,
        "lp",
        PARETO_OPTIMALITY_ONLY,
        true,
        3600,
    )

    stats = solve_instance!(model, parameter)
    @test haskey(stats.data, "Opt")
    @test stats.data["Opt"] ≈ 1
end

function test_gbc_simple_bilevel_with_blc_subsolver_auto_oracle()
    model = generate_gbc_simple_bilevel_instance_blc_subsolver_auto_oracle()
    outdir = logging_folder * "/gbc_simple_bilevel_blc_subsolver_auto_oracle"
    mkpath(outdir)
    parameter = GBCparam(
        GurobiSolver(),
        true,
        outdir,
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    stats = solve_instance!(model, parameter)
    @test haskey(stats.data, "Opt")
    @test stats.data["Opt"] ≈ 1
end

function test_blclag_simple_bilevel_with_blc_subsolver_auto_oracle()
    model = generate_blclag_simple_bilevel_instance_blc_subsolver_auto_oracle()
    outdir = logging_folder * "/blclag_simple_bilevel_blc_subsolver_auto_oracle"
    mkpath(outdir)
    parameter = BlCLagparam(
        GurobiSolver(),
        true,
        outdir,
        "lp",
        PARETO_OPTIMALITY_ONLY,
        true,
        3600,
    )

    stats = solve_instance!(model, parameter)
    @test haskey(stats.data, "Opt")
    @test stats.data["Opt"] ≈ 1
end

test_gbc_simple_bilevel_with_blc_subsolver()
test_blclag_simple_bilevel_with_blc_subsolver()
test_gbc_simple_bilevel_with_blc_subsolver_auto_oracle()
test_blclag_simple_bilevel_with_blc_subsolver_auto_oracle()
