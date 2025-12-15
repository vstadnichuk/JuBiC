# Mixed Integer Program (MIP) Solver

The Mixed Integer Program (MIP) solver provides an interface for solving bilevel problems by reformulating them as a single-level mixed-integer program. The user must perform this procedure explicitly. It provides a wrapper for passing a single-level MIP problem to JuBiC, which is especially useful if you want to benchmark other methods against an MIP solver. 


## Master Type

The `MIPMaster` type has the following fields:

- `mymip`: The MIP expression of the bilevel problem modeled using JuMP

When creating a MIP instance, set the `master` to the `MIPMaster` and set the `subproblems` to nothing.


## Solver Parameters

The MIP solver accepts the following parameters:

- `solver`: The optimization solver to be used
- `debbug_out`: Boolean flag to enable or disable the debug output
- `output_folder_path`: The path to the output directory
- `file_format_output`: The format of the output files (e.g., "lp", "mps")
- `stats`: The runtime statistics
- `runtime`: The maximum allowed runtime for the solver (in seconds)
- `threads_master`: Number of threads for the master problem

The parameters can be set using one of the following constructors:

```julia
MIPparam(solver, debbug_out, output_folder_path, file_format_output)
MIPparam(solver, debbug_out, output_folder_path, file_format_output, runtime)
MIPparam(solver, debbug_out, output_folder_path, file_format_output, stats, runtime, threads_master)
```


## Solver Output Statistics

- `Solver`: 'MIPSolver' the name of the solver used
- `time_limit`: The time limit set for the solver
- `Opt`: The objective value of the solution found
- `runtime`: The total runtime of the solver
- `gap`: The optimality gap of the solution
- `BNodes`: The number of branch-and-bound nodes explored
