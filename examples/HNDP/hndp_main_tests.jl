# Execute this file to reproduce the tests for the Paper on different path-based and arc-based models

include("hndp_tests.jl")
include("../logging.jl")

# Tests for simple instances checking basic functionalities
#test_toy_HNDPwC(["GBC", "GBCLag", "BlC", "BlCLag"])
test_negative_HNDP(["GBCLag"])

# Test for GBC and blc solvers
#test_HNDPwC("examples/HNDP/run_settings/GBCparamsTest.json")

# Tests for layer graph structure 
#test_HNDPfix("examples/HNDP/run_settings/mytests.json")






