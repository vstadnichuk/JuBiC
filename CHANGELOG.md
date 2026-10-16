# Changes

## 0.6.0 — branch integration (unreleased)

This release consolidates the four research branches below into `v0.6`, starting
from `68666a7` (the previous `main`/`v0.6` checkout). Changes are grouped by behavior,
so overlapping fixes are listed once rather than once per source commit.

### Bug fixes

- **Parallel Gurobi execution:** create worker environments before separation
  starts and bind persistent follower and connector models to their assigned
  environment. Serialize concurrent master callbacks and protect worker-owned
  mutable connector state to prevent races and unsafe environment sharing.
- **Time limits:** propagate remaining GBC time through connector LP solves,
  pricing, repeated callbacks, and Pareto refinement; cancel separation at the
  global deadline and preserve timeout status rather than treating it as an
  unrelated numerical failure.
- **Connector timing:** report disjoint LP, pricing, Pareto, and total connector
  times, avoiding double-counting work performed inside Pareto refinement.
- **Numerical cut construction:** keep fallback connector values and bound
  coefficients consistent, stabilize follower linking values and their objective
  evaluations, and retain diagnostics when a cut coefficient or connector state
  is numerically invalid.
- **Integer objectives:** round objective-derived values only when they are
  already close to an integer, and assemble paired linking terms consistently.
  A genuinely fractional value is no longer automatically promoted to an integer.
- **Persistent models:** remove temporary fixing and Pareto-band constraints on
  exceptional exits, and flush Gurobi's update bookkeeping after deleting
  temporary constraints before reusing a model.
- **Termination reporting:** retain recoverable incumbents after interrupted
  optimization, distinguish numerical termination from timeouts, and handle
  unavailable relative-gap/node-count attributes without losing other results.
- **MiBS input:** improve MPS/name handling and follower-instance validation,
  including additional parser checks and time-limit handling.
- **A* followers:** correct worker/model handling introduced by parallel
  separation and check unsupported negative transition costs before search.
- **Reproducibility:** pass connector seeds consistently and support reusing
  saved network instances across experiments and solver comparisons.
- **Cross-branch reconciliation:** preserve both cost-aware big-M callbacks and
  integer-coefficient stabilization; validate cost-aware callbacks on real
  follower solutions rather than calling them with `nothing` during setup.
- **Compatibility:** retain the former full positional `BLCparam` and
  `BlCLagparam` constructors while defaulting new worker controls to one worker.
- **Branching attributes:** set priorities through JuMP's variable-attribute
  interface so cached/direct models and repeated policy application do not
  attempt to reattach an already attached optimizer or use unmapped indices.
- **Test execution:** replace conflicting global `mktempdir` aliases with
  suite-specific helpers so GBC, BlCLag, batch, and HNDP tests can load together;
  tolerate module-qualified type names in the serialized-cut checks.
- **HNDP integration:** preserve seeded origin/destination selection and
  saved-instance support alongside the new negative-risk topology mode; correct
  the EMA reproduction launcher's relative include path. Limit the raw-data
  ignore rule to the repository root so new `examples/data/` inputs remain visible.
- **HNDP MiBS metadata:** store model-size metadata through BilevelJuMP's object
  dictionary instead of accessing the nonexistent `BilevelModel.ext` field.
- **Saved multimodal instances:** preserve decision groups, decision budgets,
  and absent weight matrices when writing and reading JSON. New files use
  format version 2 so older readers reject unsupported constraints instead of
  silently losing them. Legacy version-1 files remain readable.
- **Missing example input:** include the Friedrichshain network required by the
  named-topology generator and its regression tests, with upstream attribution.

### New features

- **Worker controls:** configure `parallel_subsolvers` separately from master and
  follower solver thread counts. Add advisory subsolver thread capabilities and
  warnings for unsupported or excessive requests.
- **Connector bounds:** configure the ConnectorLP `s` upper bound explicitly
  through `connector_s_bound`, without changing the existing default bound.
- **BlC cut controls:** enable integer-objective coefficient stabilization and
  supply big-M callbacks that depend on the current follower objective and
  leader solution, while keeping legacy callback signatures usable.
- **Branching policies:** select `default`, `linking_first`, or
  `fixed_linking_order` through solver configuration. Apply these policies
  consistently to GBC, BlC, and BlCLag and reject misspelled policy names.
- **Controlled solver experiments:** opt into the reduced Gurobi settings using
  `JUBIC_GUROBI_MINIMAL`; ordinary runs retain the existing settings.
- **Numerical diagnostics:** attach structured context to numerical exceptions,
  write termination diagnostics, and expose more detailed connector timing
  statistics.
- **HNDP models:** add cost-aware big-M rules, grouped leader decisions and
  decision budgets for multimodal bike-network and expansion models.
- **Multimodal generation:** generate car, bicycle, and transit layers from
  network attributes, including station selection, access costs, travel times,
  and expansion scenarios. Include additional named topology inputs and smoke
  and sweep runners.
- **Network persistence:** export/import generated networks as JSON, reuse saved
  graph inputs, and resume eligible experiment runs.
- **Experiment configuration:** add reusable HNDP manifests and runners,
  BOBILib compressed-input support, per-solver inputs, fresh-output sequential
  benchmark launchers, and EMA fixed-branching experiments.
- **Regression coverage:** add integration checks for constructor compatibility,
  typed cost-aware callbacks combined with integer cuts, branching priorities,
  and two concurrent follower workers.

### Example cleanup

- Remove the StochasticMultipleKnapsack example: its entry point, compact and
  decomposed model builders, SMPS reader, and example-specific tests are retired.
  No active solver or package test entry point depends on these files. No
  associated SMKP input data are present in the repository.
- The 421-test validation below predates this example cleanup.
- Remove legacy BOBILib example inputs and benchmark entry points, network-design
  benchmark configurations, and the EMA reproduction launcher.
- Remove the paper-specific HNDP seed-hunt handoff, paper manifests, and
  fixed-order reproduction batch. Retain generic experiment and solver helpers.
- Place `multimodal_example.json` and `network_generation_example.json` directly
  beside the HNDP code. Consolidate the multimodal instance/model/solver settings
  in one JSON file and let the general runner accept either that file or the
  existing three-file directory format. Update the usage documentation.
- Cross-check HNDP documentation against the builders and dispatch code: add
  formulation and generator comparison tables, BlCLag and dynamic big-M options,
  corrected topology/export lists, a local index, and runnable small examples.
- Targeted cleanup/documentation validation: 10 assertions passed, covering all
  five construction examples, save/reload, a direct strong-duality solve, both
  configuration formats, and the two successful multimodal experiment results.
  Local documentation links, builder names, and type coverage also passed.

### Integration scope and provenance

| Source branch | Integrated tip | Principal contribution |
| --- | --- | --- |
| `automated-benders-like-cut-generation-and-its-application-to-the-bilevel-network-design-problem` | `5784c8a` | Cost-aware big-M callbacks, branching controls, EMA experiments |
| `bilevel-feasible-paths-for-discrete-bilevel-network-design-problems` | `502b710` | Multimodal generation, resilient statistics (OSM/GTFS import excluded) |
| `catalog-formulations-multi-follower-discrete-bilevel-network-design` | `0a13584` | Shared HNDP paper manifests, runner, EMA topology |
| `solving-multi-follower-mixed-integer-bilevel-problems-with-binary-linking-variables` | `963e474` | Parallelism, time limits, numerical stability, parsers, saved networks |

The main README remains unchanged. Generated EMA result tables were omitted.
OSM/GTFS network construction is excluded from version 0.6: the importers,
pipeline guide, city-dependent benchmark entry points, and associated archived
result tables are removed, along with the unused `OpenStreetMapX` and `CodecXz`
direct dependencies. The historical
`examples/HNDP/aggregated_network_model_results/` folder is omitted entirely.
The active instance types, JSON persistence, grouped decisions, budgets, all
solver formulations, and TNTP-based multimodal generators remain available.
No solver implementation changes are required for this exclusion. The working
package version is `0.6.0`; this is not a published release.

### Validation

Post-removal dependency resolution and Julia syntax validation passed with
Julia 1.11.9. **421 assertions passed**: core/integration (95), HNDP formulation
and solver equivalence (252), and multimodal/JSON checks with solves (74).
Both retained TNTP-based multimodal smoke instances solved to optimality:
bike objective -18.0; expansion objective 0.0. Tests used two Julia threads,
Gurobi 13.0.3, and the repaired MiBS 1.2.2 installation. Paper-scale benchmarks
were not run. The active HNDP and core solver implementations are unchanged by
the OSM/GTFS exclusion.
