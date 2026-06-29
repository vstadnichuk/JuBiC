# The master problem for the MibS wrapper
using BilevelJuMP

"""
Master wrapper for JuBiC's direct MiBS route.

`MibSMaster` stores a `BilevelJuMP.BilevelModel` that is passed to the MiBS
interface through `MibSparam`.
"""
struct MibSMaster
    model::BilevelJuMP.BilevelModel  # The Bilevel model you want to solve
end
