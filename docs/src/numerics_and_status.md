# Numerics and Status Codes

The GBC solver route proved to be highly numerical unsatble, and JuBiC contains a variaty of 
prevention steps that detect potential numerical instability and try to medigate the negative effects. 

## Current GBC pipeline

For each integer master solution, GBC performs the following operations:

1. Master linking-variable values are rounded to binary values before they are
   passed to follower separation.
2. Each follower subsolver is solved for the current master solution. Follower
   `y`-values are expected to be binary.
3. If a returned `y`-value is not within `10^-8` of zero or one, JuBiC default MIP subsolver rounds the value to a binary value and
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
   `integer_obj=true`, objective-derived values are rounded to the nearest
   integer only when they are within `10^-4` of that integer. Coefficients are
   integerized as they are inserted into the cut using the conservative rules
   described below.
6. If Pareto refinement is enabled, it is performed after the initial
   connector solve. The pre-Pareto connector snapshot is retained so that the
   standard cut can still be constructed if refinement fails.

For HNDP A* separation, connector-based transition costs are checked before
the labeling search. A negative value with magnitude at most `10^-4` is
treated as zero. This tolerance covers insignificant negative values caused by
the numerical solution of the nonnegative ConnectorLP `k` variables. The same
clamping is applied to the transition-cost evaluation and its shortest-path
heuristic. A negative value below `-10^-4` remains an error because A* requires
nonnegative transition costs.

## Parallel Callback Calls

GBC callback entry is serialized with a global callback lock because Gurobi may
invoke callbacks from multiple master threads. Distinct follower connector
separations may still run in parallel inside one callback. Each parallel
connector uses the configured connector/subsolver thread count; with parallel
separation enabled, connector solver calls use one thread.

## ConnectorLP bound

The GBC solver accepts an optional numeric `connector_s_bound` parameter for
the upper bound of the ConnectorLP variable `s`:

```json
"connector_s_bound": 1000000
```

If the parameter is omitted, the configured generic numerical bound is used.

## Integer-objective mode

`integer_obj=true` expresses a modeling assumption that objective values used in
the relevant cuts are integer-valued. Hence, fractional values close to integer 
can be safely rounded to closest integer.

The rule is applied to follower objective values (`optL2`) before they enter
connector optimality cuts, and to objective values and associated big-M terms
in the standalone and persistent BlC subsolver cuts. 

This objective rounding is separate from coefficient integerization. For
nonnegative cut coefficients such as `k`, `xi`, and BlC terms, the cut builder
uses upward rounding when `integer_obj=true`. Paired terms are assembled before the
master expression is finalized so that the cancellation in a term of the form
`q(1-x)` is preserved. 

The final integerized optimality-cut constant is rounded to the nearest
integer by `_adjust_optcut_constant` function inside `ConnectorLP`. If the resulting
cut cannot be reconciled with the current connector solution, the run receives
the corresponding numerical warning or numerical status:

## More Numerical Warnings 

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

For GBC, `Timeout` means that the global GBC runtime was exhausted, including
the safety buffer used to stop starting new connector separations. The status
`Timeout_Submodel` is reserved for a connector or follower subsolver that
actually reaches its own time limit before the global GBC deadline is exhausted.

## Diagnostics and logs

Primary numerical diagnostics are written to the Julia log associated with the
run. Warnings for binary rounding include the affected resource and the exact
conversion, for example `y[(i,j)]=2.0e-6 -> 1.0`, together with a statement
that the objective expressions were reevaluated.

When a numerical termination is caught, JuBiC may additionally write a
`numeric_termination.json` record containing the exception type, status,
message, solver context, and accumulated statistics.

## Exception types

The subsolver and connector layers use the following exception types:

- `TimeoutException`: an oracle, connector, subsolver, or external solver
  reached its available time limit.
- `NumericalIssueException`: cut generation or connector processing encountered
  a numerical condition represented by a numerical status.
- `MibSFailureException`: a MiBS-based solve failed or did not return the
  required solution (due to some MibS internal process, check MibS log for details then).
