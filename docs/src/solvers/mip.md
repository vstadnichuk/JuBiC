# Compact MIP Wrapper (`MIP`)

`MIP` is JuBiC's direct wrapper around a JuMP MIP model. It is useful when an
instance is already available as a compact single-level model and should be
solved through JuBiC's common `solve_instance!` interface.

## Required Instance Shape

`solve_instance!(inst, params::MIPparam)` expects:

- `inst.master isa MIPMaster`
- `inst.subproblems == nothing`

`MIPMaster` stores the JuMP model in `mymip`. No decomposition-specific master,
connector, or follower subsolver objects are used.

## Solver Parameters

`MIP` is configured with `MIPparam`. The most relevant inputs are:

- `solver`: the MIP solver wrapper used to optimize the stored JuMP model.
- `debbug_out`: whether instance-level debug artifacts should be written.
- `output_folder_path`: output directory used when output logs are enabled.
- `file_format_output`: file format used if the model is exported.
- `runtime`: runtime limit in seconds.
- `threads_master`: thread limit passed to the MIP solve.

## Minimal Working Example

The following example solves a small binary MIP through JuBiC's common
`solve_instance!` interface.

```math
\begin{aligned}
\min_x\quad & x_1 + 2x_2 \\
\text{s.t.}\quad & x_1 + x_2 \ge 1, \\
& x_1,x_2 \in \{0,1\}.
\end{aligned}
```

```julia
using JuBiC
using JuMP

solver = GurobiSolver()
optimizer = () -> get_next_optimizer(solver)

model = Model(optimizer)
@variable(model, x[1:2], Bin)
@constraint(model, x[1] + x[2] >= 1)
@objective(model, Min, x[1] + 2 * x[2])

instance = Instance(MIPMaster(model), nothing)
params = MIPparam(
    solver,
    false,
    mktempdir(),
    "lp",
    60.0,
)

stats = solve_instance!(instance, params)
println(stats.data["Opt"])
```
