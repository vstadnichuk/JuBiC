### Contains fixed testsets for stochastic multiple knapsack

paths_smkp = [("examples/data/smkp/smkp_$i.cor", "examples/data/smkp/smkp_$i.sto") for i in 1:30]

"""
    test_neg_compact(time_limit::Number)

Test all 30 instances of stochastic multiple knapsack from SIPLIB (https://www2.isye.gatech.edu/~sahmed/siplib/smkp/smkp.html) using the compact formulation, i.e., one model including all scenarios.
Uses negated original variables for second-stage constraints. 
    - 'time_limit': The time limit set for each instance.
    - 'positive': If true, use original formulation without negative variables
"""
function test_neg_compact(time_limit::Number; positive=true)
    instances = []
    params = []

    # solve both models with negative and without negative variables
    myfolder = init_logging_folder()
    for i in 1:30
        paths = paths_smkp[i]
        if positive
            compact_two_stage_neg = compact_two_stage_from_file(paths[1], paths[2])
        else
            compact_two_stage_neg = compact_two_stage_from_file_neg(paths[1], paths[2])
        end
        subfolder = myfolder * "/k_$i"
        create_folder_if_not_exists(subfolder)
        inst_neg, param_neg = JuMP_to_instance(compact_two_stage_neg, "Neg_$i", subfolder; timelimit=time_limit)

        push!(instances, inst_neg)
        push!(params, param_neg)
    end
    testrun!(instances, params, myfolder)
end

"""
    test_twosatge_JuBiC(time_limit::Number)

Test all 30 instances of stochastic multiple knapsack from SIPLIB (https://www2.isye.gatech.edu/~sahmed/siplib/smkp/smkp.html) using the GBCSolver of JuBiC.
Uses negated original variables for second-stage constraints to link both stages.
    - 'time_limit': The time limit set for each instance.
"""
function test_twosatge_JuBiC(time_limit::Number)
    instances = []
    params = []

    # solve both models with negative and without negative variables
    myfolder = init_logging_folder()
    for i in 1:30
        # create Instance
        paths = paths_smkp[i]
        two_stage_JuBiC_instance = two_stage_JuBiC(paths[1], paths[2], GurobiSolver(Gurobi.Env()))

        # create Params 
        subfolder = myfolder * "/k_$i"
        create_folder_if_not_exists(subfolder)

        gbc_param = GBCparam(
            GurobiSolver(Gurobi.Env()),
            true,
            subfolder,
            "lp",
            PARETO_OPTIMALITY_ONLY,
            time_limit,
        )

        # save instance and params
        push!(instances, two_stage_JuBiC_instance)
        push!(params, gbc_param)

        # set parameter of instance
        new_stat!(get_stats(gbc_param), "Instance_name", "SMK_GBC_$i")
    end
    testrun!(instances, params, myfolder)
end




# -------------------------
# Auxiliary functions
# -------------------------
function JuMP_to_instance(model::JuMP.Model, name::String, myfolder; timelimit=300)
    create_folder_if_not_exists(myfolder)

    # create instance
    master = MIPMaster(model)
    inst = Instance(master, nothing)
    
    ca_param = MIPparam(
        GurobiSolver(Gurobi.Env()),
        true,
        myfolder,
        "lp",
        timelimit,
    )

    # set parameter of instance
    new_stat!(get_stats(ca_param), "Instance_name", name)

    # save generated and continue
    return inst, ca_param
end



