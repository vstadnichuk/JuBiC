# HNDP Solver Models

HNDP model builders translate an `HNDPwC` application instance into one of the
standard JuBiC `Instance` shapes. The solver is then called through the common
entry point:

```julia
stats = solve_instance!(instance, params)
```

The main bridge is `examples/HNDP/hndp_model_generation.jl`.

## JuBiC Decomposition Solvers

HNDP reuses JuBiC's native decomposition solvers without adding HNDP-specific
solver drivers.

- [`GBC`](../../solvers/gbc.md) is built by `build_hndp_gbc_instance(...)`.
- [`BlC`](../../solvers/blc.md) is built by `build_hndp_blc_instance(...)`.
- [`BlCLag`](../../solvers/blc_lag.md) is built by `build_hndp_blclag_instance(...)`.
- [`MiBS`](../../solvers/mibs.md) is built by `build_hndp_mibs_instance(...)`.

The HNDP code contributes the network-specific master and follower models, the
mapping from decision arcs to linking variables, big-M estimates where needed,
and optional custom follower oracles such as the A* subsolver.

## Strong-Duality Compact Model

`build_hndp_sd_instance(...)` builds a compact MIP reformulation for HNDP
instances without weight bounds.

The follower shortest-path problem can be written as a linear network-flow
problem. Because the network-flow relaxation is integral, the shortest-path
choice can be represented through primal flow variables and dual potentials. A
strong-duality equation enforces that the selected follower flow is optimal for
the user's cost objective.

For each user, the compact model contains:

- primal arc-flow variables for the user path,
- dual variables for the follower shortest-path problem,
- linking constraints that enforce `y^u_a <= x_a` on decision arcs,
- dual-feasibility constraints,
- and a strong-duality constraint equating primal and dual follower objective
  values.

The resulting model is a single-level MIP and is solved through JuBiC's
[`MIPMaster`](../../solvers/mip.md). JuBiC provides the generic MIP wrapper and
solver dispatch; the HNDP code implements the specific strong-duality
reformulation and the HNDP-specific big-M bounds.

Important options are:

- `big_m_mode`: chooses how HNDP derives big-M values. Current values include
  `HNDP_BIGM_FIXED_NETWORK_PATH` and `HNDP_BIGM_N_MINUS_ONE`.
- `indicator_constraints`: optionally replaces parts of the big-M logic by
  indicator constraints.
- `bound_duals`: bounds dual potentials using the derived big-M values.

Typical construction:

```julia
solver = GurobiSolver()
instance = build_hndp_sd_instance(
    hndp,
    solver;
    big_m_mode = HNDP_BIGM_FIXED_NETWORK_PATH,
    indicator_constraints = false,
    bound_duals = true,
)
params = MIPparam(solver, false, output_folder, "lp", 600.0)
stats = solve_instance!(instance, params)
```

## Path-Based Compact Model

`build_hndp_path_instance(...)` first enumerates feasible paths for each user
and then builds a compact path-selection MIP.

For user `u`, let `P_u` be the enumerated feasible path set. The model uses path
selection variables `z^u_p` and enforces

```math
\sum_{p \in P_u} z^u_p = 1 .
```

If path `p` uses decision arc `a`, the path can only be selected when `a` is
available:

```math
z^u_p \le x_a .
```

The objective coefficient of a path is the sum of construction-independent
operator-side contributions over the path. Construction cost remains attached
to `x`.

This route is also solved through JuBiC's [`MIPMaster`](../../solvers/mip.md).
JuBiC provides the generic compact-MIP solve path; HNDP implements path
enumeration, optional dominance pruning, and the conversion from paths to the
path-selection model.

Current builder options include:

- `enumeration_time_limit`: wall-clock budget for path-bound computation and
  path enumeration.
- `parallelize`: enumerates users in parallel during precomputation.
- `use_decision_arc_dominance`: enables the HNDP dominance rule that prunes path
  labels with no better cost/weight behavior and no better decision-arc usage.

Typical construction:

```julia
instance, enum_runtime, path_counts = build_hndp_path_instance(
    hndp,
    solver;
    enumeration_time_limit = 300.0,
    parallelize = true,
    use_decision_arc_dominance = true,
)
```

The remaining solve time is usually reduced by `enum_runtime` in benchmark
runs, so the configured total runtime still applies to the whole instance.

## Hybrid Path Models

The hybrid builders combine path enumeration with fallback formulations.

`build_hndp_hybrid_instance(...)` uses:

- path variables for users whose path sets were enumerated before the deadline,
- strong-duality fallback for users that timed out.

`build_hndp_hybrid_blc_instance(...)` uses:

- path variables inside the HPR for enumerated users,
- `BlC` subproblems for fallback users.

Both variants support:

- `enumeration_time_limit`,
- `use_decision_arc_dominance`,
- `fallback_if_paths_exceed_flow_vars`.

The last option triggers fallback even after successful enumeration if the
number of enumerated paths is larger than the number of arc-flow variables in
the fallback formulation. The benchmark metadata records these fallbacks
separately from enumeration timeouts.

## Subproblem Methods in Arc-Based Models

Arc-based HNDP builders can use different follower oracles:

- `HNDP_SUBPROBLEM_MIP`: generic JuMP follower model through
  [`SubSolverJuMP`](../../subsolvers/subsolver_jump.md).
- `HNDP_SUBPROBLEM_ASTAR`: custom HNDP A* / labeling oracle through
  [`AStarSolver`](../../subsolvers/astar.md).
- `HNDP_SUBPROBLEM_BLC_JUMP`: bilevel-capable JuMP follower through
  [`SubSolverBlCJuMP`](../../subsolvers/subsolver_blc_jump.md).
- `HNDP_SUBPROBLEM_MIBS`: bilevel follower solve through
  [`SubSolverMiBS`](../../subsolvers/subsolver_mibs.md).

The model configuration file selects these routes through fields such as
`model_type`, `subproblem_method`, `big_m_mode`, `parallelize`, and
`use_decision_arc_dominance`.
