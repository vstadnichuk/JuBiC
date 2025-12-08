# Example 2: Solving an Instance in MibS Input File Format

This example illustrates how to solve a bilevel optimization problem using the MibS Input File Format, which is widely used for defining such problems. For more details on the format, please refer to the [MibS Input File Format](https://coin-or.github.io/MibS/input.html) and the [MPS Format](https://en.wikipedia.org/wiki/MPS_(format)).

---

## Problem Definition

We will consider the following bilevel problem with binary variables. The first-level problem is formulated as:

```math
\begin{aligned}
\min \quad  & x_1 + x_2 \\
s.t. \quad  & x_1, x_2 \in \{0,1\}
\end{aligned}
```

On the other hand, the second-level problem is defined as:

```math
\begin{aligned}
\min \quad  & y_1 \\
s.t. \quad  & y_1 \leq x_1 \\
            & y_2 \leq x_2 \\
            & x_1 + y_2 \geq 1 \\
            & y_1, y_2 \in \{0,1\}
\end{aligned}
```

An optimal solution for this bilevel problem is given by $x_1 = 0$ and $x_2 = 1$.

---

## Representation in MibS Input File Format

This bilevel problem can be represented in the MibS Input File Format. In the following is the corresponding AUX file `example2.aux`:

```aux
@NAME
Example2
@MPS
example2.mps
@NUMVARS
2
@NUMCONSTRS
3
@VARSBEGIN
Y1 1.
Y2 0.
@VARSEND
@CONSTRSBEGIN
C1
C2
C3
@CONSTRSEND
```

Additionally, here is the corresponding MPS file named `example2.mps`:

```mps
NAME          Example2
ROWS
 N  OBJ
 L  C1
 L  C2
 G  C3
COLUMNS
    X1                  OBJ                                       1
    X2                  OBJ                                       1
    Y1                  C1                                        1
    X1                  C1                                       -1
    Y2                  C2                                        1
    X2                  C2                                       -1
    X1                  C3                                        1
    Y2                  C3                                        1
RHS
    rhs                 C3                                        1
BOUNDS
 BV bnd                 X1
 BV bnd                 X2
 BV bnd                 Y1
 BV bnd                 Y2
ENDATA
```

---

## Step-by-Step Implementation

### Step 1: Load Required Packages

Begin by loading the necessary packages:

```julia
using JuBiC, Gurobi
```

### Step 2: Read the Instance

Next, read the bilevel instance from the MibS input files. Make sure to specify the correct paths to your MPS and AUX files:

```julia
# Define file paths for MPS and AUX input files
mps_file = "docs/src/examples/example2.mps"
aux_file = "docs/src/examples/example2.aux"

# Read the bilevel instance
instance = get_GBC_instance(mps_file, aux_file, Gurobi.Optimizer)
```

### Step 3: Define the Parameters

Define the parameters for solving the bilevel instance using the GBC method:

```julia
# Set up the solver parameters
params = GBCparam(
    GurobiSolver(Gurobi.Env()),     # Solver used
    false,                          # Disable debug output
    "./output",                     # Output directory
    "lp",                           # Output file format
    PARETO_OPTIMALITY_ONLY          # Cut type
)

# Create the output directory
mkpath("./output")  
```

### Step 4: Solve the Instance

Now it's time to solve the instance:

```julia
# Solve the bilevel instance
result = solve_instance!(instance, params)
```

---

## Background: Transforming MibS Input File Format

It should be noted that JuBiC internally transforms the original bilevel problem formulation into an equivalent formulation with linking variables, as outlined in the [Problem Formulation](../index.md#problem-formulation) section. A key observation is that the constraint $x_1 + y_2 \geq 1$ in the second-level problem does not conform to this formulation as it directly involves the first-level variable $x_1$. To address this issue, JuBiC introduces auxiliary variables. Specifically, it creates a copy of the variable $x_1$ in the second-level problem and replaces it in the constraint accordingly. To facilitate this, we introduce the binary variable $\bar{x}_1$ in the first-level problem to represent the complement of $x_1$, which can be expressed by the equation $x_1 + \bar{x}_1 = 1$. In the second-level problem, the corresponding second-level binary variables $x_1'$ for $x_1$ and $\bar{x}_1'$ for $\bar{x}_1$ are defined. By enforcing that $\bar{x}_1'$ serves as the complement of $x_1'$, along with the binary linking constraints ($x_1' \leq x_1$ and $\bar{x}_1' \leq \bar{x}_1$), we establish that $x_1 = x_1'$.

The formalized structure of both problems is presented below. The first-level problem is defined as:

```math
\begin{aligned}
\min \quad  & x_1 + x_2 \\
s.t. \quad  & x_1 + \bar{x}_1 = 1 \\
     \quad  & x_1, \bar{x}_1, x_2 \in \{0,1\}
\end{aligned}
```

The second-level problem is reformulated as:

```math
\begin{aligned}
\min \quad  & y_1 \\
s.t. \quad  & y_1 \leq x_1 \\
            & y_2 \leq x_2 \\
            & x_1' \leq x_1 \\
            & \bar{x}_1' \leq \bar{x}_1 \\
            & x_1' + \bar{x}_1' = 1 \\
            & x_1' + y_2 \geq 1 \\
            & y_1, y_2, x_1', \bar{x}_1' \in \{0,1\}
\end{aligned}
```