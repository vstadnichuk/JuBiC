# Generalized Benders Cuts (`GBC`)

`GBC` is JuBiC's native generalized Benders decomposition solver. Technical
details of the method are described in
[Stadnichuk and Koster](https://optimization-online.org/?p=28877).

## Mathematical Structure

For a first-level solution `x`, follower `k` solves its own lower-level problem

```math
\min\{ f_k(y) : y \in Y_k,\; y_a \le x_a \ \forall a \in A \}.
```

The contribution of follower `k` to the first-level objective is then evaluated
only over lower-level optimal solutions:

```math
\varphi_k(x) =
\min\{ r_k(y) :
    y \in \arg\min\{ f_k(z) : z \in Y_k(x) \}
\}.
```

The first-level problem solved by `GBC` can therefore be written as

```math
\begin{aligned}
\min_x\quad & F(x) + \sum_{k \in \mathcal{K}} \varphi_k(x) \\
\text{s.t.}\quad & x \in X \subseteq \{0,1\}^{|A|}.
\end{aligned}
```

## Lazy Cut Structure

The core idea of `GBC` is to approximate each unknown value function ``\varphi_k``
from below by linear constraints in the first-level variables. The master
contains one auxiliary variable `subObj[k]` per follower and lazy cuts enforce
increasingly tight lower bounds on these variables.

If the current first-level solution is feasible for follower `k`, `GBC`
generates an optimality cut. These cuts generalize
[Lagrangian cuts](https://doi.org/10.1007/s10107-018-1249-5) to the bilevel
setting and are added as linear constraints of the form

```math
\alpha_k + \sum_{a \in A} \beta_{ka} x_a \le \text{subObj}_k,
```

where the coefficients define a generated linear expression in the first-level
linking variables.

If the current first-level solution is infeasible for a follower, `GBC`
generates a feasibility cut to cut off this master solution. These are
automatically generated
[combinatorial Benders cuts](https://doi.org/10.1007/978-3-540-25960-2_14)
in the first-level variables.
They are added as linear constraints of the form

```math
\gamma + \sum_{a \in A} \delta_a x_a \ge 1.
```

## Certified Inexact Connector Pricing

For connector variables ``(s,g,k)``, pricing computes

```math
p^*(g,k)=\min_{p\in\mathcal P}
\left\{r_p+g c_p+\sum_{a\in A_p}k_a\right\}.
```

The connector point satisfies every omitted row exactly when ``s\le p^*``.
A pricing solver returns a feasible objective ``U`` and certified
lower bound ``L``, so ``L\le p^*\le U``. JuBiC uses three cases:

1. If ``U<s``, the feasible solution supplies a genuinely violated row.
2. If ``L\ge s``, the connector point is certified feasible.
3. Otherwise the maximum possible omitted-row violation is

```math
\delta=\max\{0,s-L\}.
```

Without a finite lower bound, ``\delta=+\infty``. Feasibility-cut pricing is
always exact; the policies below apply only to optimality cuts. Pareto
refinement and BlC coefficient strengthening are skipped whenever the returned
bounds leave a positive connector error, because those stages require exact
connector feasibility. No subsolver-level heuristic flag is consulted.

### Safe underestimation

`CONNECTOR_UNDERESTIMATION` subtracts ``\delta`` from the complete affine GBC
expression. This is equivalent to using ``s_{\mathrm{safe}}=s-\delta``:

```math
\ell^-(x)=\ell(x)-\delta\le\varphi(x)\qquad\forall x.
```

Every retained cut remains valid for the original value-function epigraph, so
the master objective bound is a direct lower bound on the bilevel optimum. If
no finite certificate exists, JuBiC omits the uncertified GBC expression and
falls back to the follower's global lower bound.

### Tight bounded overestimation

`CONNECTOR_OVERESTIMATION` retains the unshifted expression. Globally,

```math
\ell(x)\le\varphi(x)+\delta.
```

Such cuts are bounded-inexact cuts, not necessarily global majorants. Since
ordinary `solve_sub_for_x` may return any follower-optimal solution, the cached
fixed-``x`` response used during separation need not minimize the leader
contribution among follower optima. This does not invalidate the corrected
lower bound below, but an arbitrary response cannot provide a pricing-error
bound on the final interval width. JuBiC therefore treats optimistic
tie-breaking as a separate optional final-evaluation capability rather than a
requirement for cut generation.
Because the master enforces the maximum of all retained cuts, a newer exact cut does
not overwrite an older inaccurate cut. For follower ``k`` JuBiC tracks

```math
\bar\delta_k=\max_{j\in J_k}\delta_{kj},
\qquad E=\sum_k\bar\delta_k.
```

If ``B_M`` is the inexact master's MIP objective bound, ``B_M-E`` is a valid
lower bound on the original bilevel optimum.

### Objective interval

Let ``x^I`` be the final fully separated master incumbent and let ``V^+(x^I)``
denote the optimistic value obtained by minimizing the first-level contribution
over all follower-optimal responses. After the master terminates, GBC calls the
optional `solve_sub_for_x_optimistic` method for each follower. Binary master
values are normalized to exact zeros and ones before this evaluation and before
cache lookup. When every follower implements the method, JuBiC reports

```math
\begin{array}{ll}
\text{underestimation:} & B_M\le Z^*\le V^+(x^I),\\
\text{overestimation:}  & B_M-E\le Z^*\le V^+(x^I).
\end{array}
```

Let ``G_M`` be the master absolute MIP gap, and for underestimation let
``E_I=\sum_k\delta_k(x^I)`` be the sum of the final-incumbent connector errors.
Then the interval widths satisfy

```math
\begin{array}{ll}
\text{underestimation:} & W\le G_M+E_I,\\
\text{overestimation:}  & W\le G_M+E.
\end{array}
```

Thus, if there are ``m`` followers and every relevant absolute pricing error is
at most ``\delta``, then ``W\le G_M+m\delta``. These formulas also cover early
master termination.

The optimistic method is optional. For an unsupported follower, JuBiC reuses
the exact but arbitrarily tie-broken response cached during separation. This
still gives a valid upper endpoint. However, if ``\widehat V(x^I)`` is that
response and

```math
\tau(x^I)=\widehat V(x^I)-V^+(x^I)\ge 0,
```

then the width contains the additional, generally unknown term ``\tau(x^I)``.
JuBiC emits a warning and does not claim a pricing-error width bound. A unique
follower optimum, a leader contribution constant over the follower-optimal
set, or an independent bound on ``\tau`` would provide equivalent assurance.

Relevant `RunStats` keys are:

- `MasterObjective`, `MasterObjectiveBound`, and `MasterAbsoluteGap`;
- `ConnectorPricingMaxError` and `ConnectorPricingMaxErrorBySub`;
- `ConnectorApproximationErrorBound`, equal to ``E`` for overestimation;
- `FinalConnectorErrorSum` and `FinalConnectorUpperBound`;
- `IncumbentObjectiveUpperBound`;
- `OptimisticIncumbentObjective`, finite only when every follower completed
  optimistic evaluation;
- `FinalOptimisticEvaluationUsed`, `FinalOptimisticEvaluationComplete`,
  `FinalOptimisticEvaluationBySub`,
  `FinalOptimisticEvaluationFallbackFollowers`, and
  `FinalOptimisticEvaluationTime`;
- `ObjectiveIntervalLower`, `ObjectiveIntervalUpper`,
  `ObjectiveIntervalWidth`, and `ObjectiveIntervalCertified`.
- `ObjectiveIntervalWidthBoundedByPricingError` and
  `ObjectiveIntervalWidthBound`.
- `GBCSolutionType` (`Exact` or `Heuristic`), `GBCResultStatus` (for example
  `Optimal` or `HeuristicOptimal`), `MasterSolvedToOptimality`,
  `UsedInexactPricing`, and `NInexactPricingCalls`.

`Opt` remains the historical master-surrogate objective for benchmark
compatibility. The old `HeuristicMaster*` and `ExactIncumbentObjective` keys are
also retained as compatibility aliases. New benchmark code should use
`IncumbentObjectiveUpperBound` and the explicit GBC status markers.

## Master Representation Used by `GBC`

The first level is passed as a JuMP model through `Master`.

Important inputs are:

- `model`: the JuMP model containing the first-level variables, constraints, and the direct first-level objective term `F(x)`.
- `A`: the set over which the binary linking variables are defined.
- `link_vars`: a dictionary mapping each `a in A` to the corresponding first-level variable `x_a`.
- `sub_names`: names of the follower subproblems; these names are used to match master-side objects with subsolver objects.
- `objL2`: optional expressions for the lower-level objective contribution in the master representation, when this information is available and useful to the solver.
- `partial_decomposition`: optional callback for adding additional master-side variables or constraints when only part of the follower structure is decomposed.

The master model should contain the first-level structure only. The follower
value functions are represented through the generated lazy cuts.

## Compatible Subsolvers

`GBC` is decomposition-based, so one
[subsolver](../sub_solvers.md) is required for each follower. The currently
implemented supported subsolvers are:

- `SubSolverJuMP`
- `SubSolverMiBS`
- `AStarSolver`

## Solver Parameters

`GBC` is configured with `GBCparam`. The most relevant inputs are:

- `solver`: the MIP solver wrapper used for the master problem and auxiliary JuMP models.
- `debbug_out`: whether instance-level debug artifacts should be written.
- `output_folder_path`: output directory used when output logs are enabled.
- `file_format_output`: file format used for exported JuMP models, for example `"lp"` or `"mps"`.
- `pareto`: Pareto-cut mode, for example `PARETO_OPTIMALITY_ONLY`.
- `runtime`: runtime limit in seconds.
- `warmstart`: whether connector state is reused across callback solves.
- `bigMwithLC`: whether BlC-style subroutines are used to strengthen some GBC coefficients.
- `trim_coeff`: whether generated coefficients are trimmed by available bounds.
- `parallel_separation`: whether follower-side connector separation is parallelized across multiple Julia workers.
- `threads_master` and `threads_sub_con`: thread limits for the master MIP and follower-side solves.
- `connector_add_current_solution_cut`: if enabled, the follower solution computed for the current first-level point is inserted into the corresponding `ConnectorLP` before the connector separation loop starts, unless the same row already exists.
- `subsolver_numerical_preprocessing`: if enabled, compatible subsolvers may simplify numerically extreme connector objectives before solving the pricing problem.
- `connector_approximation`: use `CONNECTOR_UNDERESTIMATION` (default) for
  globally valid weakened cuts or `CONNECTOR_OVERESTIMATION` for point-tight,
  bounded-inexact cuts.

Additional constructor variants expose seed, thread, and numerical-tolerance
settings; see [Core API Reference](../solver_api.md).

For direct construction, the four-argument convenience constructor accepts the
policy as a keyword:

```julia
params = GBCparam(
    GurobiSolver(), false, output_directory, "lp";
    connector_approximation=CONNECTOR_OVERESTIMATION,
)
```

Batch and HNDP JSON configurations use
`"connector_approximation": "underestimation"` or `"overestimation"`.

If `parallel_separation = true`, JuBiC requires `threads_sub_con = 1`. In that
mode the parallelism comes from solving several follower-side models
concurrently, not from letting each worker model use multiple Gurobi threads.

The last two options are intended for numerically difficult `GBC` runs. They
can change the internal sequence of generated connector rows, but not the
mathematical problem being solved. The first option can give the connector LP
more information early. The second option is currently implemented for
[`SubSolverJuMP`](../subsolvers/subsolver_jump.md).

## Minimal Working Example

The following example is the simple unit-test bilevel problem used in
`test/gbc.jl`.

```math
\begin{aligned}
\min_{x_1,x_2,y_1,y_2}\quad & x_1 - x_2 + 10 y_2 \\
\text{s.t.}\quad & x_1, x_2 \in \{0,1\}, \\
& (y_1, y_2) \in \arg\min \Big\{-y_1 - y_2 : \\
& \qquad y_1 \le x_1, \\
& \qquad y_2 \le x_2, \\
& \qquad y_1 = 1, \\
& \qquad y_1, y_2 \in \{0,1\}\Big\}.
\end{aligned}
```

```julia
using JuBiC
using JuMP

function build_gbc_example()
    solver = GurobiSolver()
    optimizer = () -> get_next_optimizer(solver)

    A = [1, 2]
    sub_name = "Sub0"

    master_model = Model(optimizer)
    @variable(master_model, x[A], Bin)
    @objective(master_model, Min, x[1] - x[2])
    master = Master(master_model, A, Dict(a => x[a] for a in A), [sub_name])

    sub_model = Model(optimizer)
    set_silent(sub_model)
    @variable(sub_model, y[1:2], Bin)
    @constraint(sub_model, y[1] == 1)
    follower_obj = @expression(sub_model, -y[1] - y[2])
    @objective(sub_model, Min, follower_obj)
    master_obj_term = @expression(sub_model, 10 * y[2])

    subsolver = SubSolverJuMP(
        sub_name,
        sub_model,
        A,
        y,
        master_obj_term,
        follower_obj,
    )

    return Instance(master, [subsolver]), solver
end

instance, solver = build_gbc_example()
params = GBCparam(
    solver,
    false,
    mktempdir(),
    "lp",
    PARETO_OPTIMALITY_ONLY,
)

stats = solve_instance!(instance, params)
println(stats.data["Opt"])
```
