# A wrapper for the MibS solver
using BilevelJuMP, MibS_jll

function solve_with_MibS!(inst::Instance, param::MibSparam)
    model = inst.master.model

    # register solver
    @debug "starting setup of the MibS solver."
    new_stat!(param.stats, "Solver", "MibSSolver")

    @debug "Finished model construction. Now proceeding to optimization process with MibS solver."
    solution = BilevelJuMP.solve_with_MibS(model, MibS_jll.mibs)
    status = solution.status
    objective = solution.objective

    if status == true
        new_stat!(param.stats, "Opt", objective)
    else
        @debug "The MibS solver did not find an optimal solution."
    end
end