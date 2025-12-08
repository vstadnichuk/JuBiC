# Benders-like Cuts using Lagrangian Duals (BlCLag)

The BlCLag solver extends the Benders-like Cuts (BlC) method by generating big-M values through solving corresponding Lagrangian dual problems.


## Master Type

The `BlCLagMaster` type has the following fields:

- `model`: The master problem without second level
- `A`: The set of shared resources
- `link_vars`: A dictionary mapping each shared resource to its corresponding first-level linking variable
- `sub_names`: The names of the second-level problems
- `second_level_obj`: A dictionary mapping each second-level problem name to its objective function expression


## Solver Parameters

The BlCLag solver accepts the following parameters:

- `solver`: The optimization solver to be used
- `debbug_out`: Boolean flag to enable or disable the debug output
- `output_folder_path`: The path to the output directory
- `file_format_output`: The format of the output files (e.g., "lp", "mps")
- `stats`: The runtime statistics
- `runtime`: The maximum allowed runtime for the solver (in seconds)
- `threads_master`: Number of threads for the master problem
- `threads_sub_con`: Number of threads for the subproblems
- `pareto`: The setting for the pareto cuts to be used. Possible values are:
    - `PARETO_NONE`: No pareto cuts are generated
    - `PARETO_OPTIMALITY_ONLY`: Only pareto optimality cuts are generated
    - `PARETO_OPTIMALITY_AND_FEASIBILITY`: Both pareto optimality cuts and pareto feasibility cuts are generated
- `warmstart`: Boolean flag to enable or disable warmstart in ConnectorLP (if false, ConnectorLP is reset at each iteration)
- `infinity_num`: The value used to represent infinity in the solver

The parameters can be set using one of the following constructors:

```julia
BlCLagparam(solver, debbug_out, output_folder_path, file_format_output, runtime)
BlCLagparam(solver, debbug_out, output_folder_path, file_format_output, pareto, runtime)
BlCLagparam(solver, debbug_out, output_folder_path, file_format_output, pareto, warmstart, runtime)
BlCLagparam(solver, debbug_out, output_folder_path, file_format_output, stats, runtime, threads_master, threads_sub_con, pareto, warmstart, infinity_num)
```


## Solver Output Statistics

- `Solver`: 'BlCLagSolver' the name of the solver used
- `time_limit`: The time limit set for the solver
- `NSub`: The number of second-level problems
- `BlCLagCuts`: The number of Benders-like cuts generated from Lagrangian duals
- `SepaTime`: The time spent in the separator
- `SepaTimeCut`: The time spent generating cuts in the separator
- `BlCLagStatus`: The status of the BlCLag solver termination. Possible values are:
    - `Timeout_Submodel`: The solver terminated due to a timeout
    - `Terminate`: The solver terminated normally
- `Opt`: The objective value of the solution found
- `runtime`: The total runtime of the solver
- `gap`: The optimality gap of the solution
- `BNodes`: The number of branch-and-bound nodes explored
