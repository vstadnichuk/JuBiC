using JuBiC, JuMP, Gurobi, BilevelJuMP
import JuBiC: ConnectorLP_BlC, SubSolution, SubSolver, genBenderslike_cut!

_blclag_test_tempdir() = JuBiC.repo_local_tempdir("tests", "blclag"; prefix="blclag_test")
const mktempdir = _blclag_test_tempdir

mutable struct MockBlCParetoFallbackSubSolver <: SubSolver
    name::String
    mip_model::Model
    A::Vector{Int}
    separation_blc_calls::Int
end

function JuBiC.name(sub_solver::MockBlCParetoFallbackSubSolver)
    return sub_solver.name
end

function JuBiC.check(sub_solver::MockBlCParetoFallbackSubSolver, params::SolverParam)
    return nothing
end

function JuBiC.capacity_linking(sub_solver::MockBlCParetoFallbackSubSolver, a, params::SolverParam)
    return 1
end

function JuBiC.solve_sub_for_x(
    sub_solver::MockBlCParetoFallbackSubSolver,
    xvals,
    params::SolverParam,
    time_limit,
)
    return true, 1.0, 0.0, Dict(a => 0.0 for a in sub_solver.A)
end

function JuBiC.separation_BlC!(
    sub_solver::MockBlCParetoFallbackSubSolver,
    sval,
    kvals::Dict,
    params::SolverParam,
    time_limit,
)
    sub_solver.separation_blc_calls += 1
    if sub_solver.separation_blc_calls == 1
        return SubSolution(false, 0.0, 1.0, 1.0, Int[])
    end
    error("forced pareto fallback failure")
end

function JuBiC.set_nthreads(sub_solver::MockBlCParetoFallbackSubSolver, n)
    return nothing
end

"""
    generate_blclag_simple_bilevel_instance()

Generate the corresponding high-point-relaxation instance for the BlCLag solver
on the same simple bilevel problem as
`generate_gbc_simple_bilevel_instance`.

As for the classical BlC formulation, the master is the HPR with copied
second-level variables and constraints. In contrast to `BlCMaster`, no big-M
callback is required because the BlCLag solver generates the Benders-like cuts
automatically from the Lagrangian dual.
"""
function generate_blclag_simple_bilevel_instance()
    A = [1, 2]
    nsub = "Sub0"

    hpr = Model(optimizer)
    @variable(hpr, x[A], Bin)
    @variable(hpr, yh[A], Bin)
    @constraint(hpr, hpr_link_1, yh[1] <= x[1])
    @constraint(hpr, hpr_link_2, yh[2] <= x[2])
    @constraint(hpr, hpr_fix, yh[1] == 1)
    @objective(hpr, Min, x[1] - x[2] + 10 * yh[2])
    xdict = Dict(a => x[a] for a in A)
    sub_obj_hpr = @expression(hpr, -sum(yh[A]))
    master = BlCLagMaster(
        hpr,
        A,
        xdict,
        [nsub],
        Dict(nsub => sub_obj_hpr),
    )

    sub = BilevelModel()
    @variable(Upper(sub), x[A], Bin)
    @variable(Lower(sub), y[A], Bin)
    @constraint(Lower(sub), c_link[a in A], y[a] <= x[a])
    @constraint(Lower(sub), cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    master_sub_obj = @expression(sub, 10 * y[2])
    @objective(Upper(sub), Min, sub_obj)
    @objective(Lower(sub), Min, sub_obj)
    subS = SubSolverMiBS(
        nsub,
        sub,
        A,
        x,
        y,
        master_sub_obj,
        sub_obj,
    )
    return Instance(master, [subS])
end

"""
    test_blclag_requires_bilevel_subsolver()

Check that BlCLag rejects subsolvers that cannot solve bilevel follower
problems. This protects the requirement that BlCLag only runs with subsolvers
that explicitly support bilevel follower solves.
"""
function test_blclag_requires_bilevel_subsolver()
    A = [1, 2]
    nsub = "SubInvalidBlCLag"

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
    @variable(sub, y[A], Bin)
    @constraint(sub, cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    @objective(sub, Min, sub_obj)
    master_sub_obj = 10 * y[2]
    subS = SubSolverJuMP(
        nsub,
        sub,
        A,
        y,
        master_sub_obj,
        sub_obj,
        timelimit -> (false, 0),
    )
    set_silent(sub)

    model = Instance(master, [subS])
    outdir = logging_folder * "/blclag_invalid_subsolver"
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

    err = try
        solve_instance!(model, parameter)
        nothing
    catch e
        e
    end
    @test err isa ErrorException
    @test occursin("requires a bilevel-capable subsolver", sprint(showerror, err))
end

function test_blclag_pareto_fallback_numerical_guard()
    A = [1]

    master = Model(optimizer)
    @variable(master, x[A], Bin)
    link_vars = Dict(a => x[a] for a in A)

    sub_model = Model(optimizer)
    set_silent(sub_model)
    sub_solver = MockBlCParetoFallbackSubSolver("SubBlCParetoFallback", sub_model, A, 0)

    connector_lp = Model(optimizer)
    set_silent(connector_lp)
    @variable(connector_lp, -10.0 <= s <= 10.0)
    @variable(connector_lp, k[A] >= 0)
    @constraint(connector_lp, k[1] == 0.0)

    connector = ConnectorLP_BlC(
        connector_lp,
        A,
        link_vars,
        sub_solver,
        JuBiC.ConSubsolCut_BlC[],
        Dict{Symbol,Any}(),
    )

    parameter = BlCLagparam(
        GurobiSolver(),
        false,
        mktempdir(),
        "lp",
        PARETO_OPTIMALITY_ONLY,
        true,
        60.0,
    )

    cut, pobj, cutcoeff = genBenderslike_cut!(connector, Dict(1 => 0.0), parameter, 60.0)

    @test pobj ≈ 1.0
    @test cut.constant ≈ 1.0
    @test cutcoeff[1] == 0
    @test haskey(parameter.stats.data, "Opt_status_override")
    @test parameter.stats.data["Opt_status_override"] == "Numerics"
    @test haskey(parameter.stats.data, "BlCLagStatus")
    @test parameter.stats.data["BlCLagStatus"] == "Numerics"
    @test isempty(connector.numeric_state)
end

test_blclag_requires_bilevel_subsolver()
test_blclag_pareto_fallback_numerical_guard()
