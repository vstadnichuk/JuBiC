using JuBiC, Gurobi


# Define file paths for MPS and AUX input files
mps_file = "docs/src/examples/example2.mps"
aux_file = "docs/src/examples/example2.aux"

# Read the bilevel instance
instance = get_GBC_instance(mps_file, aux_file, Gurobi.Optimizer)


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


# Solve the bilevel instance
result = solve_instance!(instance, params)
