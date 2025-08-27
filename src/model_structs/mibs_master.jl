# The master problem for the MibS wrapper
using BilevelJuMP

struct MibSMaster
    model::BilevelJuMP.BilevelModel  # The Bilevel model you want to solve
end