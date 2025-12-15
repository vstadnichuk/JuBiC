# Example 1: Simple Bilevel Problem

This example demonstrates how to define and solve a simple bilevel optimization problem using binary variables at both the first- and second-level.

---

## Problem Definition

Consider the following simple bilevel problem. The first-level problem involves two binary variables $x_1$ and $x_2$. It is formulated as follows:

```math
\begin{aligned}
\min \quad  & x_1 - x_2 + 10 y_2 \\
s.t. \quad  & x_1, x_2 \in \{0,1\}
\end{aligned}
```

Here, $y_2$ is a variable from the second-level problem that influences the first-level objective. The second-level problem involves the binary variables $y_1$ and $y_2$, and is formulated as:

```math
\begin{aligned}
\min \quad  & - y_1 - y_2 \\
s.t. \quad  & y_1 \leq x_1 \\
            & y_2 \leq x_2 \\
            & y_1 = 1 \\
            & y_1, y_2 \in \{0,1\}
\end{aligned}
```

The optimal solution for this bilevel problem is given by $x_1 = 1$ and $x_2 = 0$.

---

## Step-by-Step Implementation

### Step 1: Load Required Packages

Start by loading the necessary packages for modeling and solving the optimization problems:

```julia
using JuBiC, JuMP
using Gurobi        # Or any other supported solver
```

### Step 2: Define the First-Level Problem

Define the first-level (master) problem as follows:

```julia
# Define the set of shared resources and the subproblem name
A = [1, 2]
sub_name = "Sub0"

# Create the first-level model
master_model = Model()
set_optimizer(master_model, Gurobi.Optimizer)

# Define binary variables of first level
@variable(master_model, x[A], Bin)

# Define the first-level objective function (without the second-level variables)
@objective(master_model, Min, x[1] - x[2])

# Create a dictionary for the linking variables
xdict = Dict(a => x[a] for a in A)

# Create the Master object
master = Master(master_model, A, xdict, [sub_name])
```

### Step 3: Define the Second-Level Problem

Define the second-level (sub) problem as follows:
```julia
# Create the second-level model
sub_model = Model()
set_optimizer(sub_model, Gurobi.Optimizer)

# Define binary variables for the second level
@variable(sub_model, y[1, 2], Bin)

# Define copies of first-level (linking) variables
@variable(sub_model, xc[A], Bin)

# Define second level constraints (including linking constraints)
@constraint(sub_model, c1, y[1] <= xc[1])
@constraint(sub_model, c2, y[2] <= xc[2])
@constraint(sub_model, cB, y[1] == 1)

# Define the second-level objective function
sub_obj = @expression(sub_model, -y[1] - y[2])
@objective(sub_model, Min, sub_obj)

# Define the part of the first level objective function that depends on the second-level variables
master_sub_obj = 10 * y[2]

# Set capacities of linking constraints (here all are 1)
link_constraints_capacities = Dict(a => 1 for a in A)
```

### Step 4: Create the SubSolver

We now have to define a `SubSolver` to handle the second-level problem. Custom pricing solvers can also be implemented (see [Example 3](../examples/example3.md)), but here we will use the built-in `SubSolverJuMP` to handle the second-level problem:

```julia
sub_solver = SubSolverJuMP(
    sub_name,                           # Name of the subproblem
    sub_model,                          # The second-level model
    A,                                  # Set of shared resources
    xc,                                 # Copies of first-level (linking) variables
    y,                                  # second-level variables
    link_constraints_capacities,        # Capacities of linking constraints
    master_sub_obj,                     # Part of the first-level objective function that depends on second-level variables
    sub_obj,                            # second-level objective function
    timelimit -> (false, 0)             # Custom function which adds additional cuts (not used in this example)
)
```

### Step 5: Define the Parameters

Define the parameters for solving the bilevel instance using the GBC method and creating an output directory:

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

### Step 6: Create and Solve the Instance

Finally, create the bilevel instance and solve it:

```julia
# Create the bilevel instance
instance = Instance(master, [sub_solver])

# Solve the instance
result = solve_instance!(instance, params)  
```
