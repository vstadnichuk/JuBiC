# Direct `MiBS` Wrapper

JuBiC contains a direct wrapper for the external
[`MiBS`](https://github.com/coin-or/MibS) solver.

## Required Instance Shape

`solve_instance!(inst, params::MibSparam)` expects:

- `inst.master isa MibSMaster`
- `inst.subproblems == nothing`

`MibSMaster` stores a
[`BilevelJuMP.BilevelModel`](https://joaquimg.github.io/BilevelJuMP.jl/stable/).
JuBiC uses the BilevelJuMP-based interface to pass the model to MiBS, so this
wrapper operates on a bilevel model directly rather than on a decomposed
master/subsolver representation.

## Solver Parameters

`MiBS` is configured with `MibSparam`. The most relevant inputs are:

- `debbug_out`: whether instance-level debug artifacts should be written.
- `output_folder_path`: output directory used when output logs are enabled.
- `runtime`: requested runtime limit in seconds.

## Minimal Working Example

The following example builds a small MIP-MIP bilevel model directly as a
`BilevelJuMP.BilevelModel` and solves it through the direct MiBS wrapper.

```math
\begin{aligned}
\min_{x_1,x_2,y_1,y_2}\quad & x_1 - x_2 + 10 y_2 \\
\text{s.t.}\quad & x_1, x_2 \in \{0,1\}, \\
& (y_1, y_2) \in \arg\min \Big\{-y_1 - y_2 : \\
& \qquad y_1 \le x_1,\; y_2 \le x_2,\; y_1 = 1,\;
          y_1,y_2 \in \{0,1\}\Big\}.
\end{aligned}
```

```julia
using JuBiC
using BilevelJuMP

model = BilevelModel()
@variable(Upper(model), x[1:2], Bin)
@variable(Lower(model), y[1:2], Bin)

@constraints(Lower(model), begin
    y[1] <= x[1]
    y[2] <= x[2]
    y[1] >= 1
    y[1] <= 1
end)

@objective(Upper(model), Min, x[1] - x[2] + 10 * y[2])
@objective(Lower(model), Min, -y[1] - y[2])

instance = Instance(MibSMaster(model), nothing)
params = MibSparam(false, mktempdir(), 60.0)

stats = solve_instance!(instance, params)
println(stats.data["Opt"])
```
