# Getting Started

This page shows one small bilevel model from the test suite, first in mathematical notation and then in JuBiC syntax. The example is then solved in two ways:

- with JuBiC's generalized Benders decomposition (`GBC`),
- and with the direct `MiBS` wrapper on an explicit bilevel model.

## The Example Problem

We use a single-follower bilevel model with two binary linking variables.

```math
\begin{aligned}
\min_{x_1, x_2, y_1, y_2}\quad & x_1 - x_2 + 10 y_2 \\
\text{s.t.}\quad & x_1, x_2 \in \{0,1\}, \\
& (y_1, y_2) \in \arg\min_{\bar y_1, \bar y_2} \left\{-\bar y_1 - \bar y_2 :
\bar y_1 \le x_1,\; \bar y_2 \le x_2,\; \bar y_1 = 1,\; \bar y_1,\bar y_2 \in \{0,1\}\right\}.
\end{aligned}
```

The optimal first-level solution is:

```math
x_1 = 1,\qquad x_2 = 0,
```

with total objective value `1`.

## JuBiC Code

The following code builds the model in JuBiC syntax.

```julia
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
```

This constructs a JuBiC decomposition instance:

- `master_model` contains the first-level variables and first-level objective terms,
- `follower_model` contains the follower problem,
- `Master(...)` wraps the first-level JuMP model,
- `SubSolverJuMP(...)` wraps the follower JuMP model,
- and `Instance(master, [follower])` combines both into the object expected by `solve_instance!`.

## Structural Modeling Requirements

The following structural requirements should be taken into account when modeling with JuBiC:

- The linking variables must be binary. Here those are the `x[a]` variables in the master and the `y[a]` variables passed into the follower wrapper.
- The shared resource set `A` is the common indexing set that connects both levels.
- The line

```julia
xdict = Dict(a => x[a] for a in A)
```

defines which first-level variables are the linking variables seen by JuBiC.

- The line

```julia
SubSolverJuMP(sub_name, follower_model, A, y, first_level_follower_term, sub_obj, ...)
```

defines the follower-side linking variables. They must match the shared resource set `A`.


## Solving with `GBC`

JuBiC provides multiple solver interfaces. Here, we exemplary use the native JuBiC solver `GBC`.

```julia
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
```

The parameter object specifies how JuBiC should run this solver. In this example:

- `solver` selects the optimizer backend,
- `gbc_out` is the output directory,
- `"lp"` is the model export format if output files are enabled,
- `PARETO_OPTIMALITY_ONLY` selects the GBC cut strategy,
- and `60.0` is the runtime limit in seconds.

The output reports that the model solved to optimality and that the objective value is `1` (up to numerical precision).

```text
GBC status: Optimal
GBC objective: 1.0000999834060664
```

## Solving the Same Problem with `MiBS`

JuBiC also supports accessing existing bilevel solvers such as `MiBS`. Different solver routes rely on different structural assumptions. In particular, the direct `MiBS` route requires a `MIP-MIP` bilevel model. In JuBiC this route is built on top of the [BilevelJuMP](https://joaquimg.github.io/BilevelJuMP.jl/stable/) wrapper, so we formulate the same example a bit differently here and build it directly as a `BilevelJuMP.BilevelModel`.

```julia
using BilevelJuMP

mibs_model = BilevelModel()
@variable(Upper(mibs_model), x[1:2], Bin)
@variable(Lower(mibs_model), y[1:2], Bin)

@constraints(Lower(mibs_model), begin
    y[1] <= x[1]
    y[2] <= x[2]
    y[1] >= 1
    y[1] <= 1
end)

@objective(Upper(mibs_model), Min, x[1] - x[2] + 10 * y[2])
@objective(Lower(mibs_model), Min, -y[1] - y[2])

mibs_instance = Instance(MibSMaster(mibs_model), nothing)
mibs_out = JuBiC.repo_local_tempdir("docs", "getting_started"; prefix = "mibs")
mibs_param = MibSparam(false, mibs_out, 60.0)
mibs_param.stats.data["enable_output_logs"] = false

mibs_stats = solve_instance!(mibs_instance, mibs_param)

println("MiBS status: ", mibs_stats.data["MibSStatus"])
println("MiBS objective: ", mibs_stats.data["Opt"])
```

The same pattern applies here:

- `MibSparam(...)` provides the output directory and runtime limit,
- `enable_output_logs = false` disables solver-side debug artifacts,
- and `solve_instance!` returns the run statistics.

The output reports that `MiBS` solved the model to optimality and found objective value `1`.

```text
MiBS status: Optimal
MiBS objective: 1.0
```

## Next Steps

After this first example, the most relevant follow-up pages are:

- [Using JuBiC](model_objects.md)
- [Core Solver API](solver_api.md)
- [Generalized Benders Cuts (GBC)](solvers/gbc.md)
- [Direct MiBS Wrapper](solvers/mibs.md)

