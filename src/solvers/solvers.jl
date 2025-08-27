# You can find a list of supported solvers for hierarchical problems with two levels
# You can also find functions to directly execute them here

####### Types of Solvers available #######

# GBCSolver: Use hierarchical decomposition and generate generalized Benders cuts to inlcude lower level decisions.


####### Solve your problem by running one of these solvers #######

"""
    solve_instance(inst::Instance, params::GBCparam, hsolver::GBCSol)

Solve the passed insatnce with our GBCSolver employing hierarchical decomposition.

# Arguments
- 'inst::Instance': The problem instance.
- 'params::GBCparam': The parameters for the solver

# Returns
- 'stats': The statitics collected while solving.
"""
function solve_instance!(inst::Instance, params::GBCparam)
    logger = new_file_logger(params.output_folder_path * "/debuglogGBC.txt", params.debbug_out)
    with_logger(logger) do
        solve_with_GBC!(inst, params)
    end
    return params.stats
end

"""
    solve_instance(inst::Instance, params::BLCparam)

Solve the passed insatnce with our BlCSolver employing benders-like decomposition.

# Arguments
- 'inst::Instance': The problem instance.
- 'params::GBCparam': The parameters for the solver

# Returns
- 'stats': The statitics collected while solving.
"""
function solve_instance!(inst::Instance, params::BLCparam)
    logger = new_file_logger(params.output_folder_path * "/debuglogBLC.txt", params.debbug_out)
    with_logger(logger) do
        solve_with_BLC!(inst, params)
    end
    return params.stats
end

"""
    solve_instance(inst::Instance, params::MIPparam)

Solve the passed insatnce with our MIPsolver employing that just executes the underlying MIP solver.

# Arguments
- 'inst::Instance': The problem instance.
- 'params::MIPparam': The parameters for the solver

# Returns
- 'stats': The statitics collected while solving.
"""
function solve_instance!(inst::Instance, params::MIPparam)
    logger = new_file_logger(params.output_folder_path * "/debuglogMIP.txt", params.debbug_out)
    with_logger(logger) do
        solve_with_MIP!(inst, params)
    end
    return params.stats
end

"""
    solve_instance(inst::Instance, params::MibSparam)

Solve the passed instance with the MibS solver.

# Arguments
- 'inst::Instance': The problem instance.
- 'params::MibSparam': The parameters for the solver

# Returns
- 'stats': The statitics collected while solving.
"""
function solve_instance!(inst::Instance, params::MibSparam)
    logger = new_file_logger(params.output_folder_path * "/debuglogMibS.txt", params.debbug_out)
    with_logger(logger) do
        solve_with_MibS!(inst, params)
    end
    return params.stats
end

