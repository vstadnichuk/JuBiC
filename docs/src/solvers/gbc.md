# Generalized Benders Cuts (`GBC`)

`GBC` is JuBiC's native generalized Benders decomposition solver. Technical
details of the method are described in
[Stadnichuk and Koster](https://optimization-online.org/?p=28877).

## Mathematical Structure

For a first-level solution `x`, follower `k` solves its own lower-level problem

```math
\min\{ f_k(y) : y \in Y_k,\; y_a \le x_a \ \forall a \in A \}.
```

The contribution of follower `k` to the first-level objective is then evaluated
only over lower-level optimal solutions:

```math
\varphi_k(x) =
\min\{ r_k(y) :
    y \in \arg\min\{ f_k(z) : z \in Y_k(x) \}
\}.
```

The first-level problem solved by `GBC` can therefore be written as

```math
\begin{aligned}
\min_x\quad & F(x) + \sum_{k \in \mathcal{K}} \varphi_k(x) \\
\text{s.t.}\quad & x \in X \subseteq \{0,1\}^{|A|}.
\end{aligned}
```

## Lazy Cut Structure

The core idea of `GBC` is to approximate each unknown value function ``\varphi_k``
from below by linear constraints in the first-level variables. The master
contains one auxiliary variable `subObj[k]` per follower and lazy cuts enforce
increasingly tight lower bounds on these variables.

If the current first-level solution is feasible for follower `k`, `GBC`
generates an optimality cut. These cuts generalize
[Lagrangian cuts](https://doi.org/10.1007/s10107-018-1249-5) to the bilevel
setting and are added as linear constraints of the form

```math
\alpha_k + \sum_{a \in A} \beta_{ka} x_a \le \text{subObj}_k,
```

where the coefficients define a generated linear expression in the first-level
linking variables.

If the current first-level solution is infeasible for a follower, `GBC`
generates a feasibility cut to cut off this master solution. These are
automatically generated
[combinatorial Benders cuts](https://doi.org/10.1007/978-3-540-25960-2_14)
in the first-level variables.
They are added as linear constraints of the form

```math
\gamma + \sum_{a \in A} \delta_a x_a \ge 1.
```

## Master Representation Used by `GBC`

The first level is passed as a JuMP model through `Master`.

Important inputs are:

- `model`: the JuMP model containing the first-level variables, constraints, and the direct first-level objective term `F(x)`.
- `A`: the set over which the binary linking variables are defined.
- `link_vars`: a dictionary mapping each `a in A` to the corresponding first-level variable `x_a`.
- `sub_names`: names of the follower subproblems; these names are used to match master-side objects with subsolver objects.
- `objL2`: optional expressions for the lower-level objective contribution in the master representation, when this information is available and useful to the solver.
- `partial_decomposition`: optional callback for adding additional master-side variables or constraints when only part of the follower structure is decomposed.

The master model should contain the first-level structure only. The follower
value functions are represented through the generated lazy cuts.

## Compatible Subsolvers

`GBC` is decomposition-based, so one
[subsolver](../sub_solvers.md) is required for each follower. The currently
implemented supported subsolvers are:

- `SubSolverJuMP`
- `SubSolverMiBS`
- `AStarSolver`

## Solver Parameters

`GBC` is configured with `GBCparam`. The most relevant inputs are:

- `solver`: the MIP solver wrapper used for the master problem and auxiliary JuMP models.
- `debbug_out`: whether instance-level debug artifacts should be written.
- `output_folder_path`: output directory used when output logs are enabled.
- `file_format_output`: file format used for exported JuMP models, for example `"lp"` or `"mps"`.
- `pareto`: Pareto-cut mode, for example `PARETO_OPTIMALITY_ONLY`.
- `runtime`: runtime limit in seconds.
- `warmstart`: whether connector state is reused across callback solves.
- `bigMwithLC`: whether BlC-style subroutines are used to strengthen some GBC coefficients.
- `trim_coeff`: whether generated coefficients are trimmed by available bounds.
- `parallel_separation`: whether follower-side connector separation is parallelized across multiple Julia workers.
- `threads_master` and `threads_sub_con`: thread limits for the master MIP and follower-side solves.
- `connector_add_current_solution_cut`: if enabled, the follower solution computed for the current first-level point is inserted into the corresponding `ConnectorLP` before the connector separation loop starts, unless the same row already exists.
- `subsolver_numerical_preprocessing`: if enabled, compatible subsolvers may simplify numerically extreme connector objectives before solving the pricing problem.

Additional constructor variants expose seed, thread, and numerical-tolerance
settings; see [Core API Reference](../solver_api.md).

If `parallel_separation = true`, JuBiC requires `threads_sub_con = 1`. In that
mode the parallelism comes from solving several follower-side models
concurrently, not from letting each worker model use multiple Gurobi threads.

The last two options are intended for numerically difficult `GBC` runs. They
can change the internal sequence of generated connector rows, but not the
mathematical problem being solved. The first option can give the connector LP
more information early. The second option is currently implemented for
[`SubSolverJuMP`](../subsolvers/subsolver_jump.md).

## Minimal Working Example

The following example is the simple unit-test bilevel problem used in
`test/gbc.jl`.

```math
\begin{aligned}
\min_{x_1,x_2,y_1,y_2}\quad & x_1 - x_2 + 10 y_2 \\
\text{s.t.}\quad & x_1, x_2 \in \{0,1\}, \\
& (y_1, y_2) \in \arg\min \Big\{-y_1 - y_2 : \\
& \qquad y_1 \le x_1, \\
& \qquad y_2 \le x_2, \\
& \qquad y_1 = 1, \\
& \qquad y_1, y_2 \in \{0,1\}\Big\}.
\end{aligned}
```

```julia
using JuBiC
using JuMP

function build_gbc_example()
    solver = GurobiSolver()
    optimizer = () -> get_next_optimizer(solver)

    A = [1, 2]
    sub_name = "Sub0"

    master_model = Model(optimizer)
    @variable(master_model, x[A], Bin)
    @objective(master_model, Min, x[1] - x[2])
    master = Master(master_model, A, Dict(a => x[a] for a in A), [sub_name])

    sub_model = Model(optimizer)
    set_silent(sub_model)
    @variable(sub_model, y[1:2], Bin)
    @constraint(sub_model, y[1] == 1)
    follower_obj = @expression(sub_model, -y[1] - y[2])
    @objective(sub_model, Min, follower_obj)
    master_obj_term = @expression(sub_model, 10 * y[2])

    subsolver = SubSolverJuMP(
        sub_name,
        sub_model,
        A,
        y,
        master_obj_term,
        follower_obj,
    )

    return Instance(master, [subsolver]), solver
end

instance, solver = build_gbc_example()
params = GBCparam(
    solver,
    false,
    mktempdir(),
    "lp",
    PARETO_OPTIMALITY_ONLY,
)

stats = solve_instance!(instance, params)
println(stats.data["Opt"])
```
