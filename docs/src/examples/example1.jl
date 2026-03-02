using JuBiC, JuMP
using Gurobi        # Or any other supported solver


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


# Create the second-level model
sub_model = Model()
set_optimizer(sub_model, Gurobi.Optimizer)

# Define binary variables for the second level
@variable(sub_model, y[A], Bin)

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


# Create the bilevel instance
instance = Instance(master, [sub_solver])

# Solve the instance
result = solve_instance!(instance, params)  
