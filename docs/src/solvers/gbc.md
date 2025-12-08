# Generalized Benders Cut (GBC)

This is the hierarchical decomposition method for solving as described in the [preprint](https://optimization-online.org/?p=28877) and the major part of the implementation of the `JuBiC` package. It is a branch-and-cut method that iteratively generates second-level feasibility and optimality cuts, after solving the first-level problem.


## Master Type

The `Master` type the following fields:

- `model`: The JuMP model representing the master problem
- `A`: The set of shared resources between the first-level and second-level problems
- `link_vars`: A dictionary mapping each shared resource to its corresponding first-level linking variable
- `sub_name`: The names of the second-level problems
- `partial_decomposition`: Either `nothing` or a function
- `objL2`: A dictionary mapping each second-level problem name to its objective function

The `Master` type can be constructed using the following constructor:

```julia
Master(model, A, link_vars, sub_name)
Master(model, A, link_vars, sub_name, partial_decomposition)
Master(model, A, link_vars, sub_name, partial_decomposition, objL2)
```


## Solver Parameters

The GBC solver accepts the following parameters:

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
- `bigMwithLC`: Boolean flag to enable or disable BlC cuts
- `trim_coeff`: Boolean flag that determines whether to apply coefficient trimming optimization. When set to `true`, the solver uses bounds on cut coefficients to also trim k-coefficients and big M values generated from BlC coefficients
- `infinity_num`: The value used to represent infinity in the solver
- `g_round_digit`: The number of digits to round the coefficients in the solver

The parameters can be set using one of the following constructors:

```julia
GBCparam(solver, debbug_out, output_folder_path, file_format_output)
GBCparam(solver, debbug_out, output_folder_path, file_format_output, pareto)
GBCparam(solver, debbug_out, output_folder_path, file_format_output, pareto, runtime)
GBCparam(solver, debbug_out, output_folder_path, file_format_output, pareto, warmstart, bigMwithLC, trim_coeff, runtime)
GBCparam(solver, debbug_out, output_folder_path, file_format_output, stats, runtime, threads_master, threads_sub_con, pareto, warmstart, bigMwithLC, trim_coeff, infinity_num, g_round_digit)
```


## Solver Output Statistics

- `Solver`: 'GBCSolver' the name of the solver used
- `time_limit`: The time limit set for the solver
- `NSub`: The number of second-level problems
- `NFeasCuts`: The number of feasibility cuts generated
- `NOptCuts`: The number of optimality cuts generated
- `SepaTime`: The time spent in the separator
- `SepaTimeCut`: The time spent generating cuts in the separator
- `runtime_preprocessingGBC`: The time spent in preprocessing for GBC
- `GBCStatus`: The status of the GBC solver termination. Possible values are:
    - `Timeout_Submodel`: The solver terminated due to a timeout
    - `Terminate`: The solver terminated normally
- `Opt`: The objective value of the solution found
- `runtime`: The total runtime of the solver
- `gap`: The optimality gap of the solution
- `BNodes`: The number of branch-and-bound nodes explored
