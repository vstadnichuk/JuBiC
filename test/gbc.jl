using JuBiC, JuMP, Gurobi, BilevelJuMP
import JuBiC: ConSubsolCut, ConnectorLP, SubSolution, SubSolver, genBenders_cut!

mutable struct MockNumericDuplicateSubSolver <: SubSolver
    name::String
    mip_model::Model
    A::Vector{Int}
end

function JuBiC.name(sub_solver::MockNumericDuplicateSubSolver)
    return sub_solver.name
end

function JuBiC.check(sub_solver::MockNumericDuplicateSubSolver, params::SolverParam)
    return nothing
end

function JuBiC.capacity_linking(sub_solver::MockNumericDuplicateSubSolver, a, params::SolverParam)
    return 1
end

function JuBiC.compute_lower_bound_master_contribution(
    sub_solver::MockNumericDuplicateSubSolver,
    params::SolverParam,
    time_limit,
)
    return 0.0
end

function JuBiC.solve_sub_for_x(
    sub_solver::MockNumericDuplicateSubSolver,
    xvals,
    params::SolverParam,
    time_limit,
)
    return true, 0.0, 1.0e9, Dict(a => 0.0 for a in sub_solver.A)
end

function JuBiC.separation!(
    sub_solver::MockNumericDuplicateSubSolver,
    sval,
    gvals,
    kvals::Dict,
    params::SolverParam,
    time_limit,
)
    return SubSolution(true, 0.0, 1000.0, sval - 1.0, [1])
end

function JuBiC.set_nthreads(sub_solver::MockNumericDuplicateSubSolver, n)
    return nothing
end

"""
    generate_gbc_simple_bilevel_instance()

Generate the simple bilevel test instance used for the GBC solver tests.

The upper-level problem is:
    min x1 - x2 + 10 y2

The lower-level problem is:
    min -y1 - y2
    s.t. y1 <= x1, y2 <= x2, y1 = 1

The optimal first-level solution is x1 = 1, x2 = 0, with overall objective 1.
This function returns the decomposed instance with a `Master` and one
`SubSolverJuMP`.
"""
function generate_gbc_simple_bilevel_instance()
    A = [1, 2]
    nsub = "Sub0"

    mm = Model(optimizer)
    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1] - x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])

    sub = Model(optimizer)
    @variable(sub, y[1:2], Bin)
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
    return Instance(master, [subS])
end

"""
    generate_blc_simple_bilevel_instance()

Generate the corresponding high-point-relaxation instance for the classical BlC
solver on the same simple bilevel problem as
`generate_gbc_simple_bilevel_instance`.

The master is the high-point relaxation, i.e., it contains copies of the
second-level variables and constraints, and the subproblem is the follower model
used for cut generation. The instance has one subproblem and uses a constant big
M value of 100 in the tests.
"""
function generate_blc_simple_bilevel_instance()
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
    master = BlCMaster(
        hpr,
        A,
        xdict,
        (a, sub_name) -> 100.0,
        [nsub],
        Dict(nsub => sub_obj_hpr),
    )

    sub = Model(optimizer)
    @variable(sub, y[1:2], Bin)
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
    return Instance(master, [subS])
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
    simple_bilevel_example()

Generates a very simple MIP bilevel problem and solves it using the GBC solver.

# Description

The upper-level problem is:

    min x1 -x2 + 10y2
    s.t.
        x1, x2 in {0,1}

The lower-level problem is:

    min -y1 -y2
    s.t.
        y1 <= x1
        y2 <= x2
        y1 = 1
        y1, y2 in {0,1}

The optimal solution is: x1 = 1 and x2 = 0.
"""
function test_gbc_simple_bilevel()
    model = generate_gbc_simple_bilevel_instance()

    parameter = GBCparam(
        GurobiSolver(),
        true,
        logging_folder * "/gbc_simple_bilevel",
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    mkpath(logging_folder * "/gbc_simple_bilevel")
    stats = solve_instance!(model, parameter)

    @test haskey(stats.data, "Opt") && stats.data["Opt"] ≈ 1
end

function test_gbc_solver_instance_io_roundtrip()
    instance = generate_gbc_simple_bilevel_instance()
    export_dir = mktempdir()
    output_GBC_solver_instance(instance, export_dir)

    @test isfile(joinpath(export_dir, "master.lp"))
    @test isfile(joinpath(export_dir, "sub_Sub0.lp"))
    @test isfile(joinpath(export_dir, "metadata.json"))

    imported = read_GBC_solver_instance(export_dir, GurobiSolver())
    @test imported.master isa Master
    @test length(imported.subproblems) == 1
    @test imported.subproblems[1] isa SubSolverJuMP

    solve_dir = mktempdir()
    parameter = GBCparam(
        GurobiSolver(),
        true,
        solve_dir,
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    stats = solve_instance!(imported, parameter)
    @test haskey(stats.data, "Opt")
    @test stats.data["Opt"] ≈ 1
end


"""
    feasibility_cuts_test()

Simple test for generation of feasibility cuts. the solver is required to solve a bilevel problem with first level.

# Description

The upper-level problem is:

    min x1 + x2
    s.t.
        x1, x2 in {0,1}

The lower-level problem is:

    min 100
    s.t.
        y1 <= x1
        y2 <= x2
        y1 + y2 = 2
        y1, y2 in {0,1}

The two feasibility cuts x1 >= 1 and x2 >= 1 should be generated.
"""
function test_gbc_feasibility_cuts()
    # common ressources
    A = [1, 2]

    # names of subs
    nsub = "SubFeas0"

    # master MIP
    mm = Model(optimizer)

    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1] + x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])

    # subproblem 
    sub = Model(optimizer)

    @variable(sub, y[1:2], Bin)
    @constraint(sub, cB, y[1] + y[2] == 2)
    sub_obj = @expression(sub, AffExpr(100))
    @objective(sub, Min, sub_obj)
    master_sub_obj = AffExpr(0)
    subS = SubSolverJuMP(
        nsub,
        sub,
        A,
        y,
        master_sub_obj,
        sub_obj,
        timelimit -> (false, 0)
    )

    # set model to silent and continue
    set_silent(sub)
    model = Instance(master, [subS])

    gbc_logging_folder = logging_folder * "/gbc_feasibility_cuts"
    parameter = GBCparam(
        GurobiSolver(),
        true,
        gbc_logging_folder,
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    mkpath(gbc_logging_folder)
    stats = solve_instance!(model, parameter)


    # check whether the cuts x1 >= 1 and x2 >= 1 were generated
    cuts_collection = gbc_logging_folder * "/mastercuts_collection.txt"
    cuts_content = if isfile(cuts_collection) read(cuts_collection, String) else "" end
    @test occursin("ScalarConstraint{AffExpr, MathOptInterface.GreaterThan{Float64}}(x[1], MathOptInterface.GreaterThan{Float64}(1.0))", cuts_content)
    @test occursin("ScalarConstraint{AffExpr, MathOptInterface.GreaterThan{Float64}}(x[2], MathOptInterface.GreaterThan{Float64}(1.0))", cuts_content)
end


"""
    two_follower_test()

A simple tests for generating optimality cuts for a bilevel problem with 2 independent followern.
    
# Description

The upper-level problem is:

    min -x1 + y1 + y2
    s.t.
        x1 in {0,1}

The lower-level has 2 followers. The first is:

    min -y1
    s.t.
        y1 <= x1
        y1 in {0,1}
The second is:

    min -y2
    s.t.
        y2 <= 3x1
        y2 in {0,1}

The optimal solution is x1=0. 
"""
function test_gbc_two_follower()
    # common ressources
    A = [1]

    # names of subs
    nsub1 = "SubTwo1"
    nsub2 = "SubTwo2"

    # master MIP
    mm = Model(optimizer)

    @variable(mm, x[A], Bin)
    @objective(mm, Min, -x[1])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub1, nsub2])

    # subproblem 1
    sub1 = Model(optimizer)

    @variable(sub1, y1[[1]], Bin)
    sub_obj1 = @expression(sub1, -y1[1])
    @objective(sub1, Min, sub_obj1)
    master_sub_obj1 = 1*y1[1]
    subS1 = SubSolverJuMP(
        nsub1,
        sub1,
        A,
        y1,
        master_sub_obj1,
        sub_obj1,
        timelimit -> (false, 0)
    )

    # subproblem 2
    sub2 = Model(optimizer)

    @variable(sub2, y2[[1]], Bin)
    sub_obj2 = @expression(sub2, -y2[1])
    @objective(sub2, Min, sub_obj2)
    master_sub_obj2 = 1*y2[1]
    subS2 = SubSolverJuMP(
        nsub2,
        sub2,
        A,
        y2,
        master_sub_obj2,
        sub_obj2,
        timelimit -> (false, 0)
    )

    # disable output subs
    set_silent(sub1)
    set_silent(sub2)

    model = Instance(master, [subS1, subS2])

    parameter = GBCparam(
        GurobiSolver(),
        true,
        logging_folder * "/gbc_two_follower",
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    mkpath(logging_folder * "/gbc_two_follower")
    stats = solve_instance!(model, parameter)

    @test haskey(stats.data, "Opt") && stats.data["Opt"] ≈ 0
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

function test_gbc_duplicate_cut_numerical_guard()
    A = [1]

    master = Model(optimizer)
    @variable(master, x[A], Bin)
    link_vars = Dict(a => x[a] for a in A)

    sub_model = Model(optimizer)
    set_silent(sub_model)
    sub_solver = MockNumericDuplicateSubSolver("SubNumericDuplicate", sub_model, A)

    connector_lp = Model(optimizer)
    set_silent(connector_lp)
    @variable(connector_lp, s <= 1.0e9)
    @variable(connector_lp, k[A] >= 0)
    @variable(connector_lp, 0 <= g <= 1.0e9)
    @constraint(connector_lp, s == 1.0e9)
    @constraint(connector_lp, g == 1.0e6)
    @constraint(connector_lp, k[1] == 0.0)

    duplicate_cut = ConSubsolCut([1], 1000.0, 0.0)
    connector = ConnectorLP(
        connector_lp,
        A,
        link_vars,
        sub_solver,
        0.0,
        nothing,
        ConSubsolCut[duplicate_cut],
        0,
        Dict{Symbol,Any}(),
    )

    parameter = GBCparam(
        GurobiSolver(),
        false,
        mktempdir(),
        "lp",
        PARETO_NONE,
        60.0,
    )

    feas, cut, bigMcut, pobj = genBenders_cut!(connector, Dict(1 => 0.0), parameter, 60.0)

    @test !feas
    @test bigMcut === nothing
    @test pobj ≈ 1.0e9
    @test cut.constant > 0.0
    @test haskey(parameter.stats.data, "Opt_status_override")
    @test parameter.stats.data["Opt_status_override"] == "Opt_Numerics"
    @test haskey(parameter.stats.data, "GBCStatus")
    @test parameter.stats.data["GBCStatus"] == "Opt_Numerics"
    @test length(connector.my_subsolutions) == 1
    @test isempty(connector.numeric_state)
end


test_gbc_simple_bilevel()
test_gbc_solver_instance_io_roundtrip()
test_gbc_feasibility_cuts()
test_gbc_two_follower()
test_blclag_requires_bilevel_subsolver()
test_gbc_duplicate_cut_numerical_guard()
