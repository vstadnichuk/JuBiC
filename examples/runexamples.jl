# Add the parent directory to the LOAD_PATH
# TODO: Seems like quite the old file. Maybe delete at some point or integrate into a "Getting Started" example file 
push!(LOAD_PATH, joinpath(@__DIR__, ".."))

using JuBiC, JuMP, Gurobi
optimizer = Gurobi.Optimizer


include("logging.jl")
include("simple_bilevel_example.jl")
include("HNDP/hndp_tests.jl")

# Create a folder for the run
logging_folder = init_logging_folder()


hndpt = build_toy_HNDPwC()
myfolderBlC = logging_folder * "/BlCSolver"
mkpath(myfolderBlC)

# create instance
inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()))
blc_param = BLCparam(GurobiSolver(Gurobi.Env()), false, myfolderBlC, "lp", 10000)

# set parameter of instance
JuBiC.new_stat!(JuBiC.get_stats(blc_param), "seed", 42)

statsistics = solve_instance!(inst, blc_param)
println(statsistics)
print_stats_to_csv([statsistics], myfolderBlC * "/stats.csv")





# Define instances and parameters to be solved
instances = []
parameters = []

# Append the current examples to the instances and parameters
simple_bilevel_instance = simple_bilevel_example()
simple_bilevel_parameter = GBCparam(
    GurobiSolver(Gurobi.Env()),
    false,
    logging_folder * "/simple_bilevel_example",
    "lp",
    PARETO_OPTIMALITY_ONLY,
)
mkpath(logging_folder * "/simple_bilevel_example")
push!(instances, simple_bilevel_instance)
push!(parameters, simple_bilevel_parameter)

# Solve the instances
statsistics = []
@assert length(instances) == length(parameters)

for i in eachindex(instances)
    instance = instances[i]
    parameter = parameters[i]
    statsistic = solve_instance!(instance, parameter)
    push!(statsistics, statsistic)
end

# Print statistics to the logging folder
print_stats_to_csv(statsistics, logging_folder * "/stats.csv")
