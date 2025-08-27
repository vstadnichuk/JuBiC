# The master problem for the MIP wrapper. Just put your MIP you want to solve here

struct MIPMaster
    mymip::JuMP.Model  # The MIP you want to solve
end

# no automated check needed (at least I would not know what to check)