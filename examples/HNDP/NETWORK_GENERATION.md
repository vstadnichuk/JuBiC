# HNDP Network Generation Config

This document describes the JSON structure used by
[network_generation_example.json](run_settings/network_generation_example.json)
and
[hndp_network_generation.jl](hndp_network_generation.jl).

The generator creates HNDP network instances only. Solver-specific JuBiC models
are built later.

The intended workflow is JSON-driven: define instance families and parameter
grids in a network-generation JSON file, let the generator expand the list-valued
fields into concrete HNDP instances, and then pass these generated instances to
the model-generation and solver pipeline.

## Top-Level Structure

```json
{
  "instance_output": {
    "folder": "...",
    "graph_format": "gexf"
  },
  "parameter_seeds": [0, 1],
  "instances": [
    { "...": "..." }
  ]
}
```

### `instance_output`

- Type: object
- Optional: yes
- Meaning: controls whether generated instances are written to disk for
  debugging and inspection.

Supported fields:

- `folder`
  - Type: string
  - Meaning: base output folder. Each generated instance gets its own
    subfolder below it.
  - Reasonable values: writable paths such as
    `"examples/HNDP/generated_instances"` or `"/tmp/hndp_debug"`.
- `graph_format`
  - Type: string
  - Supported values: `"gexf"`
  - Meaning: graph export format for the instance graph.
  - Current behavior: each instance folder contains an `instance.gexf` file and
    a `users.json` companion file with the realized OD pairs and weight bounds.

### `parameter_seeds`

- Type: array of integers
- Optional: yes
- Default: `[0]`
- Meaning: random seeds used for generating arc parameters and user OD pairs.
- Reasonable values: small nonnegative integers such as `[0, 1, 2, 3, 4]`.

### `instances`

- Type: array of objects
- Optional: no
- Meaning: each object defines one batch of related HNDP network instances.

## Instance Spec Fields

Each entry in `instances` supports the following common fields.

### `name`

- Type: string
- Optional: yes
- Meaning: human-readable base name used in generated metadata and instance
  names.
- Reasonable values: short identifiers like `"csp_sioux"` or
  `"competition_small"`.

### `instance_type`

- Type: string
- Optional: no
- Supported values:
  - `"constrained_shortest_path"`
  - `"competition"`
- Meaning:
  - `constrained_shortest_path` builds a single-layer network.
  - `competition` builds either:
    - a two-layer competition network when the topology id starts with
      `"layered_"`, or
    - a single-layer network with a sampled set of decision arcs otherwise.

### `topologies`

- Type: array of strings
- Optional: yes
- Current supported values:
  - for `constrained_shortest_path`:
    - `"sioux_falls"`
  - for `competition`:
    - layered families:
      - `"layered_sioux_falls"`
      - `"layered_anaheim"`
      - `"layered_berlin_mitte_center"`
      - `"layered_ema"`
      - `"layered_friedrichshain_center"`
    - single-layer families:
      - `"sioux_falls"`
      - `"anaheim"`
      - `"berlin_mitte_center"`
      - `"ema"`
      - `"friedrichshain_center"`
- Meaning: selects the topology family used for this batch.

### `nusers`

- Type: array of integers
- Optional: yes
- Default: `[1]`
- Meaning: number of users generated in each network instance.
- Reasonable values: positive integers such as `[1]`, `[5, 10]`, `[25]`.

### `length_constrained`

- Type: array of booleans
- Optional: yes
- Default: `[true]`
- Meaning:
  - `true`: users receive a path-length bound and the lower level is a
    constrained shortest path problem.
  - `false`: no path-length bound is generated and the lower level is an
    unconstrained shortest path problem.

### `user_parameter_mode`

- Type: array of strings
- Optional: yes
- Default: `["shared"]`
- Supported values:
  - `"shared"`
  - `"per_user"`
- Meaning:
  - `shared`: all users share the same arc cost, risk, and weight matrices.
    Users differ only in origin, destination, and derived values such as the
    path-length bound.
  - `per_user`: each user gets its own randomly generated arc cost, risk, and
    weight matrices.
- Reasonable use:
  - use `shared` when users move on the same physical network with common arc
    parameters
  - use `per_user` when users represent heterogeneous commodities or
    user-specific valuations

### `max_cost`

- Type: integer
- Optional: yes
- Default: `100`
- Meaning: upper bound for randomly generated follower arc costs.
- Reasonable values: nonnegative integers such as `10`, `100`, `1000`.

### `max_risk`

- Type: integer
- Optional: yes
- Default: `100`
- Meaning: upper bound for randomly generated risk or profit coefficients
  before family-specific sign adjustments.
- Reasonable values: nonnegative integers such as `10`, `100`, `1000`.

### `max_weight`

- Type: integer
- Optional: yes
- Default: `100`
- Meaning: upper bound for randomly generated arc weights used in the
  path-length constraint.
- Reasonable values: nonnegative integers such as `10`, `100`, `1000`.

### `construction_cost`

- Type: integer
- Optional: yes
- Default: `0`
- Meaning: upper bound for random first-level construction cost on decision
  arcs. Values are drawn from `0:construction_cost`.
- Reasonable values:
  - `0` for no construction cost
  - small or medium nonnegative integers such as `10`, `50`, `100`

### `construction_cost_min`

- Type: integer
- Optional: yes
- Default: `0`
- Meaning: lower bound for random construction cost on decision arcs.
- Current behavior:
  - if omitted, the lower bound is `0`
  - if present together with `construction_cost_max`, costs are drawn from the
    inclusive range `[construction_cost_min, construction_cost_max]`

### `construction_cost_max`

- Type: integer
- Optional: yes
- Default: inherited from `construction_cost` when only that older field is used
- Meaning: upper bound for random construction cost on decision arcs.

### `availability_budget_fraction`

- Type: array of numbers
- Optional: yes
- Default: `[1.0]`
- Meaning: fraction of decision arcs that may be made available by the leader,
  measured relative to the number of decision arcs in the generated instance.
- Current behavior:
  - for a generated decision-arc set `A`, the integer budget is
    `floor(availability_budget_fraction * length(A))`
  - the generator stores both the requested fraction and the realized integer
    count in the instance metadata
  - if omitted, the effective budget is all decision arcs
  - specify `availability_budget_count` instead if the integer count should be
    fixed directly

### `availability_budget_count`

- Type: integer
- Optional: yes
- Meaning: explicit integer budget on the number of selectable decision arcs.
- Current behavior:
  - may not be specified together with `availability_budget_fraction`
  - must lie between `0` and the number of decision arcs
  - the corresponding fraction is recorded as
    `availability_budget_count / length(A)` in the metadata

### `parameter_seeds`

- Type: array of integers
- Optional: yes
- Meaning: local override of the top-level `parameter_seeds` for this instance
  spec only.

## Fields for `constrained_shortest_path`

### `alpha`

- Type: array of numbers
- Optional: yes
- Meaning: newer alias for `length_slack`.
- Current behavior:
  - if `alpha` is present, it is used instead of `length_slack`
  - values are written to metadata as both `alpha` and `length_slack` when
    `length_constrained = true`
  - if `length_constrained = false`, additional alpha values are ignored because
    no path-length bound is generated

### `length_slack`

- Type: array of numbers
- Optional: yes
- Default: `[1.0]`
- Meaning: controls how loose the path-length bound is relative to the shortest
  feasible path.
- Interpretation:
  - `0.0`: the bound is as tight as possible
  - `1.0`: the bound is maximally relaxed under the current formula
- Reasonable values: numbers in `[0, 1]`, commonly `0.0`, `0.1`, `0.25`,
  `0.5`, `1.0`
- Note: only relevant when `length_constrained = true`.

### `two_stage`

- Type: array of booleans
- Optional: yes
- Default: `[false]`
- Meaning:
  - `false`: risk and cost are generated independently
  - `true`: risk is set equal to cost, so the two levels are aligned
- Reasonable values: typically `[false]`; use `[true]` when testing a
  cooperative or aligned-objective variant.

## Fields for `competition`

### `alpha`

- Type: array of numbers
- Optional: yes
- Meaning: newer alias for `length_slack`, supported for layered competition,
  `single_layer_k_competition`, and `single_layer_k_decision_only`.
- Current behavior:
  - if `alpha` is present, it is used instead of `length_slack`
  - values are written to metadata as both `alpha` and `length_slack` when
    `length_constrained = true`
  - if `length_constrained = false`, additional alpha values are ignored because
    no path-length bound is generated

### `length_slack`

- Type: array of numbers
- Optional: yes
- Default: `[1.0]`
- Meaning: same meaning as for `constrained_shortest_path`, but applied to the
  generated competition network.
- Reasonable values: numbers in `[0, 1]`.

### `beta`

- Type: array of numbers
- Optional: yes
- Meaning: newer alias for `competitor_cost_factor`.
- Current behavior:
  - if `beta` is present, it is used instead of `competitor_cost_factor`
  - values are written to metadata as both `beta` and `competitor_cost_factor`

### `competitor_cost_factor`

- Type: array of numbers
- Optional: yes
- Default: `[0.8]`
- Meaning: multiplicative scaling applied to the user cost of decision arcs in
  the operator-controlled layer for layered competition instances, and to the
  sampled decision arcs in `single_layer_k_competition`.
- Interpretation:
  - values below `1.0` make those decision arcs cheaper for users
  - value `1.0` leaves their cost unchanged
  - values above `1.0` make them more expensive
- Reasonable values: positive numbers such as `0.5`, `0.8`, `1.0`, `1.2`
- Typical safe range: `(0, 2]`

### `single_layer_arc_mode`

- Type: array of strings
- Optional: yes
- Default: `["competition"]`
- Relevant only when:
  - `instance_type = "competition"`
  - and the topology id does **not** start with `"layered_"`
- Supported values:
  - `"competition"`
  - `"decision_only"`
- Meaning:
  - `"competition"`:
    - exactly `k` arcs are sampled as decision arcs,
    - the sampled arcs receive the `beta` cost scaling,
    - and only those sampled arcs contribute operator-side risk, with negative
      coefficients representing operator profit
  - `"decision_only"`:
    - exactly `k` arcs are sampled as decision arcs,
    - no `beta` scaling is applied,
    - sampled arcs are only the controllable arcs,
    - and operator risk remains on all arcs
- Important sign convention:
  - layered competition instances use negative operator-side coefficients on
    relevant operator arcs to represent profit earned by the operator
  - single-layer `"competition"` keeps nonzero operator-side coefficients only
    on sampled decision arcs, and these coefficients are negative
  - single-layer `"decision_only"` uses positive risk coefficients on all arcs

### `decision_arc_count`

- Type: array of integers
- Optional: yes
- Relevant only for single-layer competition instances
- Meaning: number of arcs sampled into the leader-controlled set `edgeA`
- Current behavior:
  - if the requested value exceeds the number of arcs in the base graph, the
    implementation clamps it to all arcs

## Competition Graph Styles

Generated competition metadata records:

- `competition_graph_style = "two_layer"`
  - for `layered_*` topologies
- `competition_graph_style = "single_layer_k_competition"`
  - for single-layer `k`-arc competition mode
- `competition_graph_style = "single_layer_k_decision_only"`
  - for single-layer `k`-arc decision-only mode

Additional metadata flags include:

- `layered_instance`
- `base_topology_family`
- `single_layer_arc_mode` for the single-layer competition family

## Batch Expansion

The generator creates one instance for every combination of the list-valued
parameters in a spec.

For example, if a spec contains:

```json
{
  "nusers": [5, 10],
  "length_constrained": [true, false],
  "user_parameter_mode": ["shared", "per_user"]
}
```

then this already gives `2 x 2 x 2 = 8` generated instances before considering
seeds or other fields.

## Example

```json
{
  "instance_output": {
    "folder": "examples/HNDP/generated_instances",
    "graph_format": "gexf"
  },
  "parameter_seeds": [0, 1],
  "instances": [
    {
      "name": "csp_sioux",
      "instance_type": "constrained_shortest_path",
      "topologies": ["sioux_falls"],
      "nusers": [5, 10],
      "length_constrained": [true, false],
      "length_slack": [0.2],
      "user_parameter_mode": ["shared"],
      "two_stage": [false],
      "max_cost": 100,
      "max_risk": 100,
      "max_weight": 100,
      "construction_cost": 0
    }
  ]
}
```
