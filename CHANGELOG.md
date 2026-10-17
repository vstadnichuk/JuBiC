# Changelog

## 0.6.0 (Unreleased)

Changes planned for version 0.6.0, relative to 0.5.0.

### New features and changes

- **Follower workers:** add `parallel_subsolvers` to configure follower worker
  counts independently of master and follower MIP thread limits, with thread
  capability checks and advisory warnings.
- **Cut controls:** add an explicit ConnectorLP `s` bound through
  `connector_s_bound`, integer-objective stabilization for BlC cuts, and big-M
  callbacks.
- **Branching policies:** support `default`, `linking_first`, and
  `fixed_linking_order` for GBC, BlC, and BlCLag.
- **HNDP instance generation:** add multimodal bicycle-station and transit-service
  expansion generators, grouped leader decisions and decision budgets; extend
  single-layer competition with decision-only and negative-risk options.
- **Reusable HNDP instances:** save and load complete networks as JSON, reuse
  saved instances across formulations, and resume eligible experiment runs.
- **Experiment configuration:** accept one combined instance/model/parameter JSON
  file or the existing three-file directory format. Place the small multimodal
  and general generation configurations beside the HNDP implementation.
- **Diagnostics and solver controls:** record structured numerical diagnostics
  and detailed connector timings.
- **Cleanup:** remove legacy StochasticMultipleKnapsack and BOBILib
  examples, paper-specific manifests, seed hunts, benchmark launchers, and
  archived experiment results.
- **Documentation:** update the documentation.

### Bug fixes

- **Parallel execution:** create and assign worker-specific Gurobi environments,
  serialize concurrent master callbacks, and protect mutable connector state to
  prevent races and unsafe environment sharing.
- **Time limits:** propagate remaining runtime through connector LP solves,
  pricing, callbacks, and Pareto refinement; stop separation at the deadline and
  preserve timeout status.
- **Numerical cuts:** keep rounded follower linking values consistent with their
  objective evaluations; reconcile connector fallback values and coefficients,
  and preserve paired linking terms when constructing cuts.
- **Integer objectives:** round objective-derived values only when they are
  already close to integers, instead of promoting genuinely fractional values.
  Preserve cost-aware big-M callbacks alongside coefficient stabilization and
  validate them on real follower solutions.
- **Persistent models:** remove temporary fixing and Pareto-band constraints on
  exceptional exits and flush Gurobi model updates before reusing models.
- **Statistics and termination:** avoid double-counting connector time, retain
  recoverable incumbents after interrupted solves, distinguish numerical
  termination from timeouts, and tolerate unavailable gap/node-count attributes.
- **MiBS input:** correct MPS/name handling, validate supported follower-instance
  formats, and improve time-limit handling.
- **A* search:** correct worker/model handling and reject unsupported negative
  transition costs before search.
- **Compatibility and branching:** preserve legacy positional `BLCparam` and
  `BlCLagparam` constructors. Apply branching priorities through JuMP variable
  attributes so cached/direct models and repeated application work correctly.
- **HNDP persistence and metadata:** preserve decision groups, budgets, and absent
  weight matrices in JSON. Write format version 2 while retaining version-1
  readers; store MiBS metadata through BilevelJuMP's object dictionary.
- **Reproducibility and example inputs:** propagate connector seeds consistently,
  preserve seeded origin/destination selection, include the missing
  Friedrichshain topology input, and keep example data visible to Git.
- **Regression checks:** avoid conflicting global temporary-directory helpers
  across test suites and handle module-qualified serialized-cut type names.
