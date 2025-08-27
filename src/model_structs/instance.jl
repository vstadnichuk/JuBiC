"""
A single problem instance consisting of a master and (multiple independent) subproblems. 
Note that an instance is independent of a solver.
"""
struct Instance
    master::Union{Master,BlCMaster,MIPMaster,MibSMaster}  # Master, BlCMaster, MIPmaster or MibSMaster instance
    subproblems  # List of subproblems (typically Vector{SubSolver}). Set to nothing for MIPMaster
end



# TODO: Just realized that we do have checks for Master and Subproblem but no checks if within an Instance the master and subproblems fit 
## (e.g., master and subproblem have the same number of variables). Maybe add such checks when someone finds the time for this.


