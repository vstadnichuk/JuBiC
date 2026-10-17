using JuBiC, JuMP, Gurobi

function build_simple_bilevel_instance()
    solver = GurobiSolver()
    optimizer = () -> JuBiC.get_next_optimizer(solver)

    A = [1, 2]
    sub_name = "Sub0"

    master_model = Model(optimizer)
    @variable(master_model, x[A], Bin)
    @objective(master_model, Min, x[1] - x[2])
    xdict = Dict(a => x[a] for a in A)
    master = Master(master_model, A, xdict, [sub_name])

    follower_model = Model(optimizer)
    @variable(follower_model, y[1:2], Bin)
    @constraint(follower_model, y[1] == 1)
    sub_obj = @expression(follower_model, -y[1] - y[2])
    @objective(follower_model, Min, sub_obj)

    first_level_follower_term = 10 * y[2]
    follower = SubSolverJuMP(
        sub_name,
        follower_model,
        A,
        y,
        first_level_follower_term,
        sub_obj,
        timelimit -> (false, 0),
    )

    set_silent(follower_model)
    return Instance(master, [follower]), solver
end

instance, solver = build_simple_bilevel_instance()

gbc_out = JuBiC.repo_local_tempdir("docs", "getting_started"; prefix = "gbc")
gbc_param = GBCparam(
    solver,
    false,
    gbc_out,
    "lp",
    PARETO_OPTIMALITY_ONLY,
    60.0,
)
gbc_param.stats.data["enable_output_logs"] = false

gbc_stats = solve_instance!(instance, gbc_param)

println("GBC status: ", gbc_stats.data["Opt_status"])
println("GBC objective: ", gbc_stats.data["Opt"])

result = gbc_stats
