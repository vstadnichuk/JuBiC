# Benders-like Cuts (`BlC`)

`BlC` is JuBiC's basic implementation of Benders-like decomposition where the
user provides the big-M coefficients used in the Benders-like cuts. The method
follows the Benders-like cuts approach of Israeli and Wood.

## Mathematical Structure

`BlC` starts from the high-point relaxation (HPR) of the full bilevel problem.
The HPR contains the upper-level constraints together with the feasibility
constraints of all follower problems, while follower optimality is enforced
later through lazy cuts. For JuBiC's multi-follower setting, the HPR therefore
contains the first-level variables `x` and copied follower variables `\hat y^k`
for all followers `k`.

The relaxation has the structure

```math
\begin{aligned}
\min_{x,\hat y}\quad
    & F(x) + \sum_{k \in \mathcal{K}} r_k(\hat y^k) \\
\text{s.t.}\quad
    & x \in X, \\
    & \hat y^k \in Y_k, && k \in \mathcal{K}, \\
    & \hat y^k_a \le x_a, && a \in A,\; k \in \mathcal{K}.
\end{aligned}
```

The HPR is one relaxation of the full bilevel model; it is not a separate
relaxation per follower. It omits the conditions that the copied follower
solutions `\hat y^k` are optimal for their follower objectives. `BlC` restores
these conditions lazily by approximating the follower value functions with
linear constraints.

For an incumbent first-level solution `\bar x`, JuBiC solves the corresponding
follower problem and obtains a follower optimum `\bar y^k`. If the HPR solution
violates follower optimality, `BlC` adds a cut of the form

```math
f_k(\hat y^k)
\le
f_k(\bar y^k) + \sum_{a \in A} M_{ka} \bar y^k_a (1 - x_a).
```

The user must provide valid coefficients `M_{ka}`.

## JuBiC Object Structure

`solve_instance!(inst, params::BLCparam)` expects:

- `inst.master isa BlCMaster`
- `inst.subproblems` to be a collection of compatible follower subsolvers

The master model passed to `BlCMaster` is the HPR model. The nontrivial input is
the big-M callback:

```julia
big_m(a, sub_name)
```

This function must return the big-M coefficient for resource `a` and follower
`sub_name`. JuBiC calls it while generating BlC cuts. The callback must be
defined for every `a in A` and every name in `sub_names`.

## Compatible Subsolvers

`BlC` uses follower subsolvers to evaluate follower optima at incumbent
first-level solutions. The currently implemented supported subsolvers are:

- `SubSolverJuMP`
- `SubSolverMiBS`
- `AStarSolver`

## Solver Parameters

`BlC` is configured with `BLCparam`. The most relevant inputs are:

- `solver`: the MIP solver wrapper used for the HPR master and follower subsolvers.
- `debbug_out`: whether instance-level debug artifacts should be written.
- `output_folder_path`: output directory used when output logs are enabled.
- `file_format_output`: file format used for exported JuMP models, for example `"lp"` or `"mps"`.
- `runtime`: runtime limit in seconds.
- `parallel_separation`: whether follower separations are solved in parallel for one incumbent master solution.
- `threads_master` and `threads_sub_con`: thread limits for the master MIP and subproblem-side work.

The big-M coefficients are not part of `BLCparam`; they are passed through the
`big_m(a, sub_name)` function stored in `BlCMaster`.

## Minimal Working Example

The example solves

```math
\begin{aligned}
\min_{x_1,x_2,y_1,y_2}\quad & x_1 - x_2 + 10 y_2 \\
\text{s.t.}\quad & x_1, x_2 \in \{0,1\}, \\
& (y_1, y_2) \in \arg\min \Big\{-y_1 - y_2 : \\
& \qquad y_1 \le x_1,\; y_2 \le x_2,\; y_1 = 1,\;
          y_1,y_2 \in \{0,1\}\Big\}.
\end{aligned}
```

In the `BlC` interface, the master model is the high-point relaxation.

```julia
using JuBiC
using JuMP

function build_blc_example()
    solver = GurobiSolver()
    optimizer = () -> get_next_optimizer(solver)

    A = [1, 2]
    sub_name = "Sub0"

    hpr = Model(optimizer)
    @variable(hpr, x[A], Bin)
    @variable(hpr, yh[A], Bin)
    @constraint(hpr, yh[1] <= x[1])
    @constraint(hpr, yh[2] <= x[2])
    @constraint(hpr, yh[1] == 1)
    @objective(hpr, Min, x[1] - x[2] + 10 * yh[2])
    second_level_obj = @expression(hpr, -yh[1] - yh[2])

    master = BlCMaster(
        hpr,
        A,
        Dict(a => x[a] for a in A),
        (a, sub_name) -> 100.0,
        [sub_name],
        Dict(sub_name => second_level_obj),
    )

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

instance, solver = build_blc_example()
params = BLCparam(
    solver,
    false,
    mktempdir(),
    "lp",
    600.0,
)

stats = solve_instance!(instance, params)
println(stats.data["Opt"])
```
