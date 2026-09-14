# Numerics and Status Codes

JuBiC's decomposition solvers generate cuts from connector LPs and follower
subproblems. The master solver, connector LPs, follower solvers, and optional
Pareto refinements therefore exchange values that may be affected by solver
tolerances. Hence, numerical handling is essential and part of the normal solve pipeline.

## Current GBC pipeline

For each integer master solution, GBC performs the following operations:

1. Master linking-variable values are rounded to binary values before they are
   passed to follower separation.
2. Each follower subsolver is solved for the current master solution. Follower
   `y`-values are expected to be binary.
3. If a returned `y`-value is not within `10^-8` of zero or one, JuBiC rounds the value to a binary value and
   reevaluates both follower objective expressions using the rounded vector.
   The returned `y`-values and objective values are consequently consistent.
4. The connector LP is solved iteratively. A violated follower solution adds a
   connector constraint of the form

   ```math
   s-\sum_{a\in A'}k_a-\alpha g\le r,
   ```

   where `A'` is the resource pattern of the follower solution, `alpha` is its
   follower-objective value, and `r` is its first-level contribution.
5. The connector solution is used to construct an optimality or feasibility
   cut. Master and follower binary patterns are rounded for the cut. When
   `integer_obj=true`, coefficients are integerized as they are inserted into
   the cut.
6. If Pareto refinement is enabled, it is performed after the initial
   connector solve. The pre-Pareto connector snapshot is retained so that the
   standard cut can still be constructed if refinement fails.

GBC callback entry is serialized with a global callback lock because Gurobi may
invoke callbacks from multiple master threads. Distinct follower connector
separations may still run in parallel inside one callback. Each parallel
connector uses the configured connector/subsolver thread count; with parallel
separation enabled, connector solver calls use one thread.

## Coefficient and cut checks

Optimality-cut coefficients are expected to be nonnegative where required by
the formulation. Very small negative values caused by floating-point noise may
be treated as zero. A material sign violation raises
`NumericalIssueException` with status `Terminate_Numerics`.

JuBiC also checks whether a generated cut reproduces the incumbent reference
value implied by its construction. A discrepancy is recorded as a numerical
warning and contributes to the numerical status of the run.

The same resource can occur in more than one optimality-cut term, for example
in both a `k`-term and a `y`- or BlC-based term. Such overlap is retained as
part of the current cut-generation procedure but is reported as a warning for
inspection.

Note that all these situations can occur in a correct run, and are 
on their own no indicators for a false results. Hence, we only throw a warning.

## Repeated connector cuts

Generated connector cuts are retained when connector warm starts are enabled.
If separation returns a resource pattern and objective pair already present in
the connector, JuBiC checks whether the repetition can be explained by solver
tolerances.

For `alpha > 0`, the required value of `g` for the repeated constraint to be
nonviolated is

```math
g_{\mathrm{req}}=
\frac{s-\sum_{a\in A'}k_a-r}{\alpha}.
```

JuBiC estimates connector-side and follower-subsolver-side tolerances from the
scale of the connector row and the subsolver's absolute and relative MIP
tolerances. The larger estimate is used as the admissible `g`-movement. If the
current `g` plus this tolerance reaches `g_req`, the repetition is accepted as
numerically explainable, the duplicate constraint is not added again, and the
run receives a numerical status. If it cannot be explained, JuBiC raises a
numerical exception to prevent unresolved cycling.

## Pareto-refinement fallback

Pareto refinement constrains the connector's original objective to remain in a
tolerance band around the value obtained before refinement. If the refined
connector becomes infeasible or otherwise fails, JuBiC restores the retained
pre-Pareto snapshot and constructs the standard cut. The run is marked as
numerically affected even when the master solver subsequently reaches an
optimal solution.

## Status reporting

The final high-level result is stored in `Opt_status`. Solver-specific fields
such as `GBCStatus`, `BlCStatus`, `BlCLagStatus`, and `MibSStatus` provide the
corresponding method-level status.

The principal numerical statuses are:

- `Opt_Numerics`: an optimal master result was reached after numerical fallback
  or numerical warnings.
- `Timelimit_Numerics`: the runtime limit was reached after numerical fallback
  or numerical warnings.
- `Terminate_Numerics`: the solve terminated because a numerical condition was
  considered unsafe, such as an unexplained duplicate connector cut, an
  infeasible connector LP, or a material coefficient violation.

Ordinary statuses such as `Optimal`, `Timelimit`, and `Terminate` remain
available when no numerical annotation applies. A numerical annotation does
not by itself prove that the returned point is invalid; it indicates that the
result passed through a condition requiring numerical fallback or additional
verification.

## Diagnostics and logs

Primary numerical diagnostics are written to the Julia log associated with the
run. Warnings for binary rounding include the affected resource and the exact
conversion, for example `y[(i,j)]=2.0e-6 -> 1.0`, together with a statement
that the objective expressions were reevaluated.

When a numerical termination is caught, JuBiC may additionally write a
`numeric_termination.json` record containing the exception type, status,
message, solver context, and accumulated statistics. This file is intended for
cases where the log alone is insufficient.

## Exception types

The subsolver and connector layers use the following exception types:

- `TimeoutException`: an oracle, connector, subsolver, or external solver
  reached its available time limit.
- `NumericalIssueException`: cut generation or connector processing encountered
  a numerical condition represented by a numerical status.
- `MibSFailureException`: a MiBS-based solve failed or did not return the
  required solution.

Solver drivers catch these exceptions, preserve the underlying message and
context in the logs or diagnostic record, and propagate the corresponding
status to the result statistics.
