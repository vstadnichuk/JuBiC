# Benders-like Cuts (BlC)

The Benders-like Cuts (BlC) solver first solves the High Point Relaxation of the bilevel problem, initially ignoring the second-level optimization functions. It then iteratively generates Benders-like cuts to ensure the optimality of the second-level problems using a user-provided big-M function explicitly provided by the user.


## Master Type

The `BlCMaster` type has the following fields:

- `hpr`: The High Point Relaxation JuMP model of the bilevel problem
- `A`: The set of shared resources
- `link_vars`: A dictionary mapping each shared resource to its corresponding first-level linking variable
- `big_m`: A function computing the big-M values for second-level problems given a shared resource $a \in A$, and sub-problem name
- `sub_names`: The names of the second-level problems
- `sub_objectives`: A dictionary mapping each second-level problem name to its objective function expression


## Solver Parameters

The BlC solver accepts the following parameters:

- `solver`: The optimization solver to be used
- `debbug_out`: Boolean flag to enable or disable the debug output
- `output_folder_path`: The path to the output directory
- `file_format_output`: The format of the output files (e.g., "lp", "mps")
- `stats`: The runtime statistics
- `runtime`: The maximum allowed runtime for the solver (in seconds)
- `threads_master`: Number of threads for the master problem
- `threads_sub_con`: Number of threads for the subproblems

The parameters can be set using one of the following constructors:

```julia
BLCparam(solver, debbug_out, output_folder_path, file_format_output)
BLCparam(solver, debbug_out, output_folder_path, file_format_output, runtime)
BLCparam(solver, debbug_out, output_folder_path, file_format_output, stats, runtime, threads_master, threads_sub_con)
```


## Solver Output Statistics

- `Solver`: 'BlCSolver' the name of the solver used
- `time_limit`: The time limit set for the solver
- `NSub`: The number of second-level problems
- `Blcuts`: The number of Benders-like cuts generated
- `SepaTime`: The time spent in the separator
- `BlCStatus`: The status of the BlC solver termination. Possible values are:
    - `Timeout_Submodel`: The solver terminated due to a timeout
    - `Terminate`: The solver terminated normally
- `Opt`: The objective value of the solution found
- `runtime`: The total runtime of the solver
- `gap`: The optimality gap of the solution
- `BNodes`: The number of branch-and-bound nodes explored
