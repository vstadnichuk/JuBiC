using JuBiC, JuMP, Gurobi


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
    # common ressources
    A = [1, 2]

    # names of subs
    nsub = "Sub0"

    # master MIP
    mm = Model(optimizer)

    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1] - x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])

    # subproblem 
    sub = Model(optimizer)

    @variable(sub, y[1:2], Bin)
    @variable(sub, xc[A], Bin)
    @constraint(sub, c1, y[1] <= xc[1])
    @constraint(sub, c2, y[2] <= xc[2])
    @constraint(sub, cB, y[1] == 1)
    sub_obj = @expression(sub, -sum(y[A]))
    @objective(sub, Min, 1 * sub_obj)
    master_sub_obj = 10 * y[2]
    link_constraints_capacities = Dict(1 => 1.0, 2 => 1.0)
    subS = SubSolverJuMP(
        nsub,
        sub,
        A,
        xc,
        y,
        link_constraints_capacities,
        master_sub_obj,
        sub_obj,
        timelimit -> (false, 0)
    )

    # set model to silent and continue
    set_silent(sub)
    model = Instance(master, [subS])

    parameter = GBCparam(
        GurobiSolver(Gurobi.Env()),
        true,
        logging_folder * "/gbc_simple_bilevel",
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    mkpath(logging_folder * "/gbc_simple_bilevel")
    stats = solve_instance!(model, parameter)

    @test haskey(stats.data, "Opt") && stats.data["Opt"] ≈ 1
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
    @variable(sub, xc[A], Bin)
    @constraint(sub, c1, y[1] <= xc[1])
    @constraint(sub, c2, y[2] <= xc[2])
    @constraint(sub, cB, y[1] + y[2] == 2)
    sub_obj = @expression(sub, AffExpr(100))
    @objective(sub, Min, sub_obj)
    master_sub_obj = AffExpr(0)
    link_constraints_capacities = Dict(1 => 1.0, 2 => 1.0)
    subS = SubSolverJuMP(
        nsub,
        sub,
        A,
        xc,
        y,
        link_constraints_capacities,
        master_sub_obj,
        sub_obj,
        timelimit -> (false, 0)
    )

    # set model to silent and continue
    set_silent(sub)
    model = Instance(master, [subS])

    gbc_logging_folder = logging_folder * "/gbc_feasibility_cuts"
    parameter = GBCparam(
        GurobiSolver(Gurobi.Env()),
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
    @variable(sub1, x1c[A], Bin)
    @constraint(sub1, c1, y1[1] <= x1c[1])
    sub_obj1 = @expression(sub1, -y1[1])
    @objective(sub1, Min, sub_obj1)
    master_sub_obj1 = 1*y1[1]
    link_constraints_capacities1 = Dict(1 => 1.0)
    subS1 = SubSolverJuMP(
        nsub1,
        sub1,
        A,
        x1c,
        y1,
        link_constraints_capacities1,
        master_sub_obj1,
        sub_obj1,
        timelimit -> (false, 0)
    )

    # subproblem 2
    sub2 = Model(optimizer)

    @variable(sub2, y2[[1]], Bin)
    @variable(sub2, x2c[A], Bin)
    @constraint(sub2, c1, y2[1] <= x2c[1])
    sub_obj2 = @expression(sub2, -y2[1])
    @objective(sub2, Min, sub_obj2)
    master_sub_obj2 = 1*y2[1]
    link_constraints_capacities2 = Dict(1 => 3.0)
    subS2 = SubSolverJuMP(
        nsub2,
        sub2,
        A,
        x2c,
        y2,
        link_constraints_capacities2,
        master_sub_obj2,
        sub_obj2,
        timelimit -> (false, 0)
    )

    # disable output subs
    set_silent(sub1)
    set_silent(sub2)

    model = Instance(master, [subS1, subS2])

    parameter = GBCparam(
        GurobiSolver(Gurobi.Env()),
        true,
        logging_folder * "/gbc_two_follower",
        "lp",
        PARETO_OPTIMALITY_ONLY,
    )

    mkpath(logging_folder * "/gbc_two_follower")
    stats = solve_instance!(model, parameter)

    @test haskey(stats.data, "Opt") && stats.data["Opt"] ≈ 0
end


test_gbc_simple_bilevel()
test_gbc_feasibility_cuts()
test_gbc_two_follower()
