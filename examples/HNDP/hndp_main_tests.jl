# Execute this file to reproduce the tests for the Paper on different path-based and arc-based models

include("hndp_tests.jl")
include("../logging.jl")

# Tests for simple instances checking basic functionalities
#test_toy_HNDPwC()
#test_negative_HNDP()

# Tests with range restrictions for users
unumber = 350
test_HNDPwC(unumber, [0.4, 0.5, 0.6, 0.7, 0.8], 10; hsolver=["GBC", "BlC"], time_limit=600, two_stage=false, partial_decomposition=false, constrac_cost=unumber*30)

# Tests for layer graph structure 
#test_HNDPfix("examples/HNDP/run_settings/mytests.json")






