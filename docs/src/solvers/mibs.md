# Mixed Integer Bilevel Solver (MibS)

The MibS solver serves as a wrapper around the bilevel optimization solver [MibS_jll.jl](https://github.com/JuliaBinaryWrappers/MibS_jll.jl/), which itself wraps the underlying [MibS solver](https://github.com/coin-or/MibS) developed by COIN-OR.


## Master Type

The `MibSMaster` type has the following fields:

- `model`: The bilevel model expressed using BilevelJuMP

When creating a MibS instance, set the `master` to the `MibSMaster` and set the `subproblems` to nothing.


## Solver Parameters

The MibS solver accepts the following parameters:

- `output_folder_path`: The path to the output directory
- `file_format_output`: The format of the output files (e.g., "lp", "mps")
- `stats`: The runtime statistics (default: `JuBiC.RunStats()`)

The parameters can be set using one of the following constructors:

```julia
MibSparam(output_folder_path, file_format_output)
MibSparam(output_folder_path, file_format_output, stats)
```


## Solver Output Statistics

- `Solver`: 'MibSSolver' the name of the solver used
- `Opt`: The objective value of the solution found
- `runtime`: The total runtime of the solver
- `gap`: The optimality gap of the solution
- `BNodes`: The number of branch-and-bound nodes explored
