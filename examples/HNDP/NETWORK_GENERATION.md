# HNDP Network Generation Config

This document describes the JSON structure used by
[network_generation_example.json](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/run_settings/network_generation_example.json)
and
[hndp_network_generation.jl](/home/stadnichuk/Documents/JuBiC/JuBiC/examples/HNDP/hndp_network_generation.jl).

The generator creates HNDP network instances only. Solver-specific JuBiC models
are built later.

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
  - `competition` builds a layered network where the competitor layer contains
    the decision arcs.

### `topologies`

- Type: array of strings
- Optional: yes
- Current supported values:
  - for `constrained_shortest_path`: `"sioux_falls"`
  - for `competition`: `"layered_sioux_falls"`
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

### `parameter_seeds`

- Type: array of integers
- Optional: yes
- Meaning: local override of the top-level `parameter_seeds` for this instance
  spec only.

## Fields for `constrained_shortest_path`

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

### `length_slack`

- Type: array of numbers
- Optional: yes
- Default: `[1.0]`
- Meaning: same meaning as for `constrained_shortest_path`, but applied to the
  layered competition network.
- Reasonable values: numbers in `[0, 1]`.

### `competitor_cost_factor`

- Type: array of numbers
- Optional: yes
- Default: `[0.8]`
- Meaning: multiplicative scaling applied to the cost of decision arcs in the
  competitor layer.
- Interpretation:
  - values below `1.0` make competitor-layer decision arcs cheaper
  - value `1.0` leaves their cost unchanged
  - values above `1.0` make them more expensive
- Reasonable values: positive numbers such as `0.5`, `0.8`, `1.0`, `1.2`
- Typical safe range: `(0, 2]`

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
