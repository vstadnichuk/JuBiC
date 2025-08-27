"""
    simple_bilevel_example()

Generates a very simple MIP bilevel problem.

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
function simple_bilevel_example()
    # common ressources
    A = [1, 2]

    # names of subs
    nsub = "Sub0"

    # master MIP
    mm = Model()

    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1] - x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])
    set_optimizer(master.model, optimizer)  # test to set optimizer after model construction

    # subproblem 
    sub = Model()
    set_optimizer(sub, optimizer)

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
    )

    # set model to silent and continue
    set_silent(sub)
    rV = Instance(master, [subS])
    return rV
end


"""
    feasibility_cuts_test()

Simple test for generation of feasibility cuts. the solver is required to solve a bilevel problem with first level \n 
min x1 + x2 \n 
s.t. x in {0,1} \n 
and the second level \n 
min 100 \n
s.t. y1 <= x1 \n
y2 <= x2 \n
y1 + y2 = 2 \n 
y in {0,1} \n 
The two feasibility cuts x1 >= 1 and x2 >= 1 should be generated.
"""
function feasibility_cuts_test()
    # common ressources
    A = [1, 2]

    # names of subs
    nsub = "SubFeas0"

    # master MIP
    mm = Model()
    set_optimizer(mm, optimizer) 

    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1]+x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub])

    # subproblem 
    sub = Model()
    set_optimizer(sub, optimizer)


    @variable(sub, y[1:2], Bin)
    @variable(sub, xc[A], Bin)
    @constraint(sub, c1, y[1] <= xc[1])
    @constraint(sub, c2, y[2] <= xc[2])
    @constraint(sub, cB, y[1] + y[2] == 2)
    sub_obj = @expression(sub, 100)
    @objective(sub, Min, sub_obj)
    master_sub_obj = 0
    link_constraints_capacities = Dict(1 => 1.0 , 2 => 1.0)
    subS = SubSolverJuMP(nsub, sub, A, xc, y, link_constraints_capacities, master_sub_obj, sub_obj)

    # set model to silent and continue
    set_silent(sub)
    return Instance(master, [subS])
end


"""
    two_follower_test()

A simple tests for generating optimality cuts for a bilevel problem with 2 independent followern. The first level is \n 
min -x1 + y1 + y2 \n 
s.t. x1 in {0,1} \n 
The second level has 2 follower. The first is \n 
min -y1 \n
s.t. y1 <= x1 \n
y1 in {0,1} \n 
The second is \n
min -y2 \n
s.t. y2 <= 3*x1 \n
y2 in {0,1} \n 
The optimal solution is x1=0. 
"""
function two_follower_test()
    # common ressources
    A = [1]

    # names of subs
    nsub1 = "SubTwo1"
    nsub2 = "SubTwo2"

    # master MIP
    mm = Model()
    set_optimizer(mm, optimizer) 

    @variable(mm, x[A], Bin)
    @objective(mm, Min, x[1])
    xdict = Dict(a => x[a] for a in A)
    master = Master(mm, A, xdict, [nsub1, nsub2])

    # subproblem 1
    sub1 = Model()
    set_optimizer(sub1, optimizer)

    @variable(sub1, y1[[1]], Bin)
    @variable(sub1, xc[A], Bin)
    @constraint(sub1, c1, y1[1] <= xc[1])
    sub_obj1 = @expression(sub1, -y1[1])
    @objective(sub1, Min, sub_obj1)
    master_sub_obj1 = -y1[1]
    link_constraints_capacities1 = Dict(1 => 1.0)
    subS1 = SubSolverJuMP(nsub1, sub1, A, xc, y1, link_constraints_capacities1, master_sub_obj1, sub_obj1)

    # subproblem 2
    sub2 = Model()
    set_optimizer(sub2, optimizer)

    @variable(sub2, y2[[1]], Bin)
    @variable(sub2, xc[A], Bin)
    @constraint(sub2, c1, y2[1] <= xc[1])
    sub_obj2 = @expression(sub2, -y2[1])
    @objective(sub2, Min, sub_obj2)
    master_sub_obj2 = -3*y2[1]
    link_constraints_capacities2 = Dict(1 => 1.0)
    subS2 = SubSolverJuMP(nsub2, sub2, A, xc, y2, link_constraints_capacities2, master_sub_obj2, sub_obj2)

    # disable output subs
    set_silent(sub1)
    set_silent(sub2)

    return Instance(master, [subS1, subS2])
end
















