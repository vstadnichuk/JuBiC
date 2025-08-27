### Main file to execute tests for stochastic multiple knapsack instances
using Gurobi 
using JuBiC

include("read_smps.jl")
include("knapsack_models.jl")
include("../HNDP/hndp_tests.jl")
include("../logging.jl")
include("knapsack_tests.jl")




compact_two_stage_model = compact_two_stage_from_file("examples/data/smkp/smkp_1.cor", "examples/data/smkp/smkp_1.sto")
write_to_file(compact_two_stage_model, "test_read_with_L2.lp")

#compact_two_stage_neg = compact_two_stage_from_file_neg("examples/data/smkp/smkp_1.cor", "examples/data/smkp/smkp_1.sto")
#write_to_file(compact_two_stage_neg, "test_read_with_L2_neg.lp")

# test construction of GBCSolver instances
#=instance1 = two_stage_JuBiC("examples/data/smkp/smkp_1.cor", "examples/data/smkp/smkp_1.sto", GurobiSolver(Gurobi.Env()))
write_to_file(instance1.master.model, "test_read_jmaster.lp")
for sub in instance1.subproblems 
    write_to_file(sub.mip_model, "test_read_S$(sub.name).lp")
end =#


# solve both models with negative and without negative variables
#test_neg_compact(50)
#test_twosatge_JuBiC(50)






