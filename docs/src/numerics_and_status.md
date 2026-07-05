# Numerics and Status Codes

JuBiC's decomposition algorithms generate lazy cuts from auxiliary optimization
problems. This creates numerical challenges that are less visible in a single
compact MIP solve:

- master solutions that should be binary may arrive from the solver with small fractional noise,
- follower and connector objectives are compared across several models,
- repeated cuts can occur when a connector problem is numerically close to a previously generated solution,
- Pareto refinement can fail even when the non-Pareto connector solution is usable,
- and external solvers such as MiBS can fail or time out independently of the master MIP.

JuBiC addresses these issues with explicit rounding, coefficient validation,
fallbacks, duplicate-cut handling, and status propagation. The main mechanisms
are summarized below.

## Connector Row Seeding

In `GBC`, the connector LP is separated after the current first-level solution
has already been evaluated by the follower oracle. If
`GBCparam.connector_add_current_solution_cut=true`, JuBiC reuses this known
follower solution and adds the corresponding connector row before the connector
separation loop starts, unless that row is already present.

This can provide useful information earlier in the connector solve. It can also
increase the number of connector rows, because the current follower solution is
not necessarily a maximal or otherwise strongest row for the connector LP.

## Binary Rounding and Cut Validation

Some cut-generation routines need binary master or follower patterns even though
the solver may return values such as `0.999999999` or `1e-9`.

JuBiC rounds such values to the integer values they are intended to represent
before using them in cut construction. When `GBCparam.integer_obj=true`, JuBiC
also treats follower risk objectives as integer-valued and rounds generated cut
constants and coefficients consistently with that assumption.

If validation detects a material inconsistency, JuBiC throws
`NumericalIssueException(..., "Terminate_Numerics")`.

## Coefficient Safeguards

Generated cut coefficients are expected to satisfy sign conditions. JuBiC checks
these assumptions before accepting a cut. Very small sign violations caused by
numerical noise may be snapped to zero, while material violations terminate the
current solve path.

## Subsolver Numerical Preprocessing

Some connector LP solutions contain coefficients at the artificial numerical
upper bound `GBCparam.infinity_num`. Passing these coefficients directly into a
subsolver objective can dominate the intended objective structure and produce
unstable pricing behavior.

If `GBCparam.subsolver_numerical_preprocessing=true`, compatible subsolvers may
handle such cases explicitly. For `SubSolverJuMP`, this means:

- if the connector value `g` is at the numerical upper bound, solve with the original follower objective and evaluate the full connector objective afterward;
- otherwise, if some `k_a` values are at the numerical upper bound, temporarily forbid the corresponding follower-side linking variables;
- if the temporary forbidding step makes the subproblem infeasible, remove it and solve the original connector pricing problem.

The option is a numerical safeguard for difficult separation calls. It does not
change the mathematical instance.

## Duplicate-Cut Handling

JuBiC stores generated connector cuts. If a newly separated cut was already
generated before, JuBiC checks whether the repetition can be explained by
admissible numerical movement in the connector variables.

A repeated connector cut means that the separation oracle returns the same
resource pattern and objective coefficients although the corresponding connector
constraint was already added. Conceptually, the repeated constraint has the
form

```math
s - \sum_{a \in A'} k_a - \alpha g \le r,
```

where `A'` is the resource pattern of the repeated subproblem solution, `r` is
its first-level contribution, and ``\alpha`` is its follower-objective value.

For ``\alpha > 0``, JuBiC computes the value of ``g`` that would make this cut just
nonviolated:

```math
g_{\mathrm{req}} =
\frac{s - \sum_{a \in A'} k_a - r}{\alpha}.
```

It then estimates how much ``g`` may move because of connector and subsolver
tolerances. First, JuBiC defines a row scale

```math
R =
\max\left\{
1,\ |s|,\ \sum_{a \in A} |k_a|,\ |\alpha g|,\ |r|
\right\}.
```

The connector-side tolerance is

```math
\Delta g_{\mathrm{con}} =
\frac{10^{-6} R}{|\alpha|}.
```

The subsolver-side tolerance uses the absolute and relative optimality
tolerances of the subsolver. If `opt` is the objective value returned by the
subsolver, then

```math
\Delta g_{\mathrm{sub}} =
\frac{
\max\{\varepsilon_{\mathrm{abs}},
      \varepsilon_{\mathrm{rel}}\max(1,|\text{opt}|)\}
}{|\alpha|}.
```

JuBiC uses

```math
\Delta g =
\max\{\Delta g_{\mathrm{con}}, \Delta g_{\mathrm{sub}}\}.
```

If

```math
g + \Delta g \ge g_{\mathrm{req}},
```

then the repeated cut can be explained by numerical tolerance. JuBiC accepts the
connector solution with the stabilized value ``g + \Delta g`` and marks the run
with numerical status. Otherwise, JuBiC treats the repetition as unresolved
cycling and throws
`NumericalIssueException(..., "NumericalIssue_DuplicateCut")`.

## Pareto Fallbacks

After solving a connector LP, JuBiC may run Pareto refinement to strengthen a
generated cut. Pareto refinement keeps the original connector objective fixed
within a tolerance band while optimizing a secondary criterion:

```math
\text{current\_obj} - \epsilon
\le
\text{lp\_obj}
\le
\text{current\_obj} + \epsilon.
```

The relevant tolerances are `GBCparam.pareto_band_tolerance`,
`GBCparam.blc_pareto_band_tolerance`, and
`BlCLagparam.blc_pareto_band_tolerance`.

If Pareto refinement fails, JuBiC restores the pre-Pareto connector solution,
continues with the standard cut, and marks the run with numerical status
(`Opt_Numerics` when the underlying master solve is otherwise optimal).

## Cut-Overlap Warnings

In `GBC`, symmetric connector solutions can create unintuitive cuts where the
same resource appears in several cut terms. These cuts are not necessarily
wrong, but they can indicate that the algorithm is not behaving as intended on
that instance. JuBiC marks these runs with a warning for later inspection.

## MiBS Failure and Timeout Propagation

`SubSolverMiBS` and the direct MiBS wrapper distinguish between:

- timeout,
- MiBS execution failure,
- and successful completion.

Timeouts raise `TimeoutException`. MiBS stderr failures raise
`MibSFailureException`. The surrounding solver driver catches these exceptions
and records the corresponding status.

## Status Fields

JuBiC writes the final high-level solve status into:

- `Opt_status`

Solver families may also write more specific status keys:

- `GBCStatus`
- `BlCStatus`
- `BlCLagStatus`
- `MibSStatus`

The most relevant reported statuses are:

- `Optimal`: the solver reports an optimal solution.
- `Opt_Numerics`: the solver reports an optimal solution after using an explicit numerical fallback.
- `Timelimit`: the solver reached the runtime limit.
- `Timelimit_Numerics`: the solver reached the runtime limit after numerical fallback handling.

Other values indicate that the run terminated through a specific solver-side
condition or exception. Current values include:

- `Timeout_Submodel`: a submodel or follower evaluation timed out.
- `Timeout_Subsolver`: a bilevel-capable subsolver timed out.
- `Terminate`: the solver terminated after an unexpected internal error.
- `Terminate_MibS`: MiBS failed inside a subsolver call.
- `Terminate_Numerics`: numerical validation failed in a way JuBiC treats as unsafe.
- `OptimizeNotCalled`: the master optimization was not called or did not start normally.
- `NumericalIssue_DuplicateCut`: duplicate-cut analysis found an unexplained repeated connector cut.

## Exception Types

The subsolver layer defines custom exception types used for status propagation:

- `TimeoutException`: an oracle, connector, subsolver, or external solver reached the available time limit.
- `NumericalIssueException`: cut generation detected a numerical issue that is represented as a solver status.
- `MibSFailureException`: a MiBS-based solve failed or did not return the required solution.

Solver drivers catch these exceptions and translate them into the status fields
above.
