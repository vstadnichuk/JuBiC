# The master problem for the MIP wrapper. Just put your MIP you want to solve here

"""
Master wrapper for JuBiC's compact MIP route.

`MIPMaster` contains a JuMP model that should be solved directly by the MIP
solver configured through `MIPparam`.
"""
struct MIPMaster
    mymip::JuMP.Model  # The MIP you want to solve
end

# no automated check needed (at least I would not know what to check)
