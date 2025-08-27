#Functions providing wrappers for tests of implementation and algorithms for various HNDP instances and solvers
using JuBiC
using Gurobi
using Logging
using CSV
using DataFrames
using JSON
include("hndp_solvers.jl")



function testrun!(instances, parameters, logfolder)
    # TODO: Right place for function? Is not limited towards HNDP?
    @assert length(instances) == length(parameters)
    logfile = logfolder * "/stats.csv"

    # start solving 
    my_stats = []
    last_save_time = time()
    for (ind, inst) in enumerate(instances)
        try
            para = parameters[ind]
            statsistics = solve_instance!(inst, para)
            push!(my_stats, statsistics)

            # make a backup every 5 minutes (or longer if single run takes longer)
            if last_save_time + 300 < time()
                @info "Doing an intermediate save of results till now. Note that we overwrite the file $logfile"
                print_stats_to_csv(my_stats, logfile)
                last_save_time = time()
            end
        catch err
            @error "When solving instance number $ind we incurred the error $err. Continue solving other instances. "
            @error stacktrace()
            #rethrow(err)
        end
    end

    print_stats_to_csv(my_stats, logfile)
end


"""
    test_toy_HNDPwC(hsolver =["GBC", "BlC"], time_limit=3600, partial_decomposition = false)

Create and solve a simple HNDPwC instance. Mainly intended as test instance which can be verified by hand. 
"""
function test_toy_HNDPwC(hsolver=["GBC", "BlC"], time_limit=3600, partial_decomposition=false)
    instances = []
    parameters = []

    myfolder = init_logging_folder()
    logger = JuBiC.new_file_logger(myfolder * "/debuglogInstances.txt", true)
    @info "This debuger only contains information on the generation of the instances and models used in the function 'test_toy_HNDPwC'. For loggers of the different model runs, look in the respective folders. "
    
    with_logger(logger) do
        hndpt = build_toy_HNDPwC()

        if "GBC" in hsolver
            myfolderGBC = myfolder * "/GBCSolver"
            create_folder_if_not_exists(myfolderGBC)

            # create instance
            inst = to_GBCInstance(
                hndpt,
                GurobiSolver(Gurobi.Env());
                partial_dec=partial_decomposition,
            )
            gbc_param = GBCparam(
                GurobiSolver(Gurobi.Env()),
                true,
                myfolderGBC,
                "lp",
                PARETO_OPTIMALITY_ONLY,
                time_limit,
            )

            # set parameter of instance
            new_stat!(get_stats(gbc_param), "seed", 42)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, gbc_param)
        end

        # build BlC
        if "BlC" in hsolver
            myfolderBlC = myfolder * "/BlCSolver"
            create_folder_if_not_exists(myfolderBlC)

            # create instance
            inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()))
            blc_param =
                BLCparam(GurobiSolver(Gurobi.Env()), true, myfolderBlC, "lp", time_limit)

            # set parameter of instance
            new_stat!(get_stats(blc_param), "seed", 42)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, blc_param)
        end
    end # end logger

    # test Sioux toy instance
    testrun!(instances, parameters, myfolder)
end


"""
    test_negative_HNDP(hsolver=["GBC", "BlC", "CA", "CP", "CH"], time_limit=3600, partial_decomposition=false, cycle_free_GBC=false)

Test the algorithm for a toy example of the HNDP with negative cycle according to first-level objective.
# Arguments
    - 'hsolver': The type of solver used for the bilevel problem. Currently supported are "GBC" for GBCSolver, "BlC" for the BlCSolver, 
        "CA" for compact arc-based model, "CP" for compact path-based model, and "CH" for a hybrid version of the previous two.
        If multiple solvers are passed, solve instance with each. 
    - 'time_limit': Time limit for the solvers.
    - 'partial_decomposition': If true, employ partial decomposition for GBCSolver.
    - 'cycle_free_GBC': If true, additionaly generate Benders-like cuts within the MIP subproblem for GBCSolver. Has no effect on other solvers.
"""
function test_negative_HNDP(hsolver=["GBC", "BlC", "CA", "CP", "CH"], time_limit=3600, partial_decomposition=false, cycle_free_GBC=true)
    instances = []
    parameters = []

    myfolder = init_logging_folder()
    logger = JuBiC.new_file_logger(myfolder * "/debuglogInstances.txt", true)
    @info "This debuger only contains information on the generation of the instances and models used in the function 'test_negative_HNDP'. For loggers of the different model runs, look in the respective folders. "
    
    with_logger(logger) do
        hndpt = build_toy_HNDPneg()

        # build Hierarchical Decomposition model
        if "GBC" in hsolver
            myfolderGBC = myfolder * "/GBCSolver"
            create_folder_if_not_exists(myfolderGBC)

            # create instance
            inst = to_GBCInstance(
                hndpt,
                GurobiSolver(Gurobi.Env());
                partial_dec=partial_decomposition, 
                MIPsubsolver = true, # we need MIP subsolver for negative cycles right now
                cycle_free_sub = cycle_free_GBC
            )
            gbc_param = GBCparam(
                GurobiSolver(Gurobi.Env()),
                true,
                myfolderGBC,
                "lp",
                PARETO_OPTIMALITY_ONLY,
                time_limit,
            )

            # set parameter of instance
            new_stat!(get_stats(gbc_param), "seed", 42)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, gbc_param)
        end

        # build BlC
        if "BlC" in hsolver
            myfolderBlC = myfolder * "/BlCSolver"
            create_folder_if_not_exists(myfolderBlC)

            # create instance
            inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()); MIPsubsolver = true) # we need MIP subsolver for negative cycles right now
            blc_param =
                BLCparam(GurobiSolver(Gurobi.Env()), true, myfolderBlC, "lp", time_limit)

            # set parameter of instance
            new_stat!(get_stats(blc_param), "seed", 42)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, blc_param)
        end

        # build arc compact solver instance
        if "CA" in hsolver
            myfolderCA = myfolder * "/CASolver"
            create_folder_if_not_exists(myfolderCA)

            # create instance
            inst = to_MIPInstance_arc(
                hndpt,
                GurobiSolver(Gurobi.Env())
            )
            ca_param = MIPparam(
                GurobiSolver(Gurobi.Env()),
                true,
                myfolderCA,
                "lp",
                time_limit,
            )

            # set parameter of instance
            new_stat!(get_stats(ca_param), "seed", 42)
            new_stat!(get_stats(ca_param), "CompactType", "Arc")

            # save generated and continue
            push!(instances, inst)
            push!(parameters, ca_param)
        end 
    
        # build arc compact solver instance
        if "CP" in hsolver
            myfolderCP = myfolder * "/CPSolver"
            create_folder_if_not_exists(myfolderCP)

            # create instance
            inst, enumtime, _, _, _ = to_MIPInstance_path(
                hndpt,
                GurobiSolver(Gurobi.Env()), 
                time_limit
            )
            cp_param = MIPparam(
                GurobiSolver(Gurobi.Env()),
                true,
                myfolderCP,
                "lp",
                time_limit-enumtime,
            )

            # set parameter of instance
            new_stat!(get_stats(cp_param), "seed", 42)
            new_stat!(get_stats(cp_param), "CompactType", "Path")
            new_stat!(get_stats(cp_param), "Enumtime", enumtime)


            # save generated and continue
            push!(instances, inst)
            push!(parameters, cp_param)
        end 

        # build arc compact solver instance
        if "CH" in hsolver
            myfolderCP = myfolder * "/CHSolver"
            create_folder_if_not_exists(myfolderCP)

            # create instance
            inst, enumtime, nbad_users = to_MIPInstance_hybrid(
                hndpt,
                GurobiSolver(Gurobi.Env()), 
                0.1 # time for enumeration of paths. You can play around with it to obtain paths or arcs for the one user we have in the toy example
            )
            ch_param = MIPparam(
                GurobiSolver(Gurobi.Env()),
                true,
                myfolderCP,
                "lp",
                time_limit-enumtime,
            )

            # set parameter of instance
            new_stat!(get_stats(ch_param), "seed", 42)
            new_stat!(get_stats(ch_param), "CompactType", "Hybrid")
            new_stat!(get_stats(ch_param), "nbadusers", nbad_users)
            new_stat!(get_stats(ch_param), "Enumtime", enumtime)

            

            # save generated and continue
            push!(instances, inst)
            push!(parameters, ch_param)
        end 
    end # end logger

    # test Sioux toy instance
    testrun!(instances, parameters, myfolder)
end


"""
    test_HNDPwC(users, alphas, nruns::Number; hsolver=["GBC", "BlC"], time_limit=3600, partial_decomposition=false)

Basic test for the Sioux Falls network Hazmat Network Design Problem with Capacity Constraints.
Test all combinations of the following parameters.
- 'users': A list with the number of users for the individual runs.
- 'alphas': The different alpha values that may appear. The alpha should be a parameter between 0 and 1, and sets the capacity limit for the test instance. 0 indicates very tight bounds, and 1 that there is no capacity limit. 
- 'nruns': Number of runs with different seed for each instance.
- 'hsolver': The type of solver used for the bilevel problem. Currently supported are "GBC" for GBCSolver and "BlC" for the BlCSolver. 
If multiple solvers are passed, solve instance with each. 
- 'time_limit': The time limit set on (each) solver (in sec).
- 'partial_decomposition': If true, apply partial decomposition where applicable, i.e., when using our GBCSolver
- 'two_stage': If true, use same objective function for risk and cost, i.e., both levels cooperate to find cost-minimal paths 
- 'constrac_cost': The cost for including an arc are generated randomly within 0:constructioncost
"""
function test_HNDPwC(users, alphas, nruns::Number; hsolver=["GBC", "BlC"], time_limit=3600, partial_decomposition=false, two_stage=false, constrac_cost=0)
    myfolder = init_logging_folder()
    instances = []
    parameters = []

    # generate the underlying HNDPwC instances
    hndps = Dict(
        (u, al, nr) => build_random_SiouxFalls(u, al; seed=nr, two_stage=two_stage, constructioncost=constrac_cost) for
        (u, al, nr) in Base.product(users, alphas, 1:nruns)
    )

    # build GBC instances
    if "GBC" in hsolver
        myfolderGBC = myfolder * "/GBCSolver"
        create_folder_if_not_exists(myfolderGBC)
        for u in users
            for al in alphas
                for nr = 1:nruns
                    # create output folder
                    myfolderrun = myfolderGBC * "/S$(u)_$(al)_$(nr)"
                    create_folder_if_not_exists(myfolderrun)

                    # create instance
                    hndpt = hndps[u, al, nr]
                    inst = to_GBCInstance(
                        hndpt,
                        GurobiSolver(Gurobi.Env());
                        partial_dec=partial_decomposition,
                    )
                    gbc_param = GBCparam(
                        GurobiSolver(Gurobi.Env()),
                        false,
                        myfolderrun,
                        "lp",
                        PARETO_OPTIMALITY_ONLY,
                        time_limit,
                    )

                    # set parameter of instance
                    new_stat!(get_stats(gbc_param), "U", u)
                    new_stat!(get_stats(gbc_param), "alpha", al)
                    new_stat!(get_stats(gbc_param), "constructioncost", constrac_cost)
                    new_stat!(get_stats(gbc_param), "partial_decomposition", partial_decomposition)
                    new_stat!(get_stats(gbc_param), "seed", nr)

                    # save generated and continue
                    push!(instances, inst)
                    push!(parameters, gbc_param)
                end
            end
        end
    end

    # build BlC
    if "BlC" in hsolver
        myfolderBlC = myfolder * "/BlCSolver"
        create_folder_if_not_exists(myfolderBlC)
        for u in users
            for al in alphas
                for nr = 1:nruns
                    # create output folder
                    myfolderrun = myfolderBlC * "/S$(u)_$(al)_$(nr)"
                    create_folder_if_not_exists(myfolderrun)

                    # create instance
                    hndpt = hndps[u, al, nr]
                    inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()))
                    blc_param = BLCparam(GurobiSolver(Gurobi.Env()), false, myfolderrun, "lp", time_limit)

                    # set parameter of instance
                    new_stat!(get_stats(blc_param), "U", u)
                    new_stat!(get_stats(blc_param), "alpha", al)
                    new_stat!(get_stats(blc_param), "constructioncost", constrac_cost)
                    new_stat!(get_stats(blc_param), "seed", nr)

                    # save generated and continue
                    push!(instances, inst)
                    push!(parameters, blc_param)
                end
            end
        end
    end

    # test Sioux
    testrun!(instances, parameters, myfolder)
end



"""
    test_HNDPfix(json_file_path)
Basic test for the Sioux Falls network Hazmat Network Design Problem with fixed network.
Test all combinations of the following parameters, which should be passed over the json file given by path 'json_path'.
# Parameters
    - 'users': A list with the number of users for the individual runs.
    - 'nruns': Number of runs with different seed for each instance.
    - 'betas': A list of beta vaölues employed when generating instances 

    - 'hsolver': The type of solver used for the bilevel problem. Currently supported are "GBC" for GBCSolver, "BlC" for the BlCSolver, 
    "CA" for compact arc-based model, "CP" for compact path-based model, and "CH" for a hybrid version of the previous two.
    If multiple solvers are passed, solve instance with each. 
    - 'bigMsetting': A list of setting to try for generating different kinds of big M. Only considered if model uses big M's. Options are ["F", "FB", "C", "CB", "I", "IFB", "ICB"]:
        -- "F": Use fixed network for big M if available (BLC and CA solver)
        -- "FB": Additionaly provide bounds on dual variables
        -- "C": Basic big M's. 
        -- "CB": Basic big M's with dual variable bounds
        -- "I": Use indicator constraints for dual.
        -- "IFB": Combine indicator constraints with "F" dual vairable bounds.
        -- "ICB": Combine indicator constraints with "C" dual vairable bounds
    - 'time_limit': The time limit set on (each) solver (in sec).
    - 'time_enum_hybrid': Time which is at most invested in hybrid model for path enumeration.
    - 'partial_decomposition': If true, apply partial decomposition where applicable, i.e., when using our GBCSolver
    - 'cycle_free_GBC': If true, additionaly generate Benders-like cuts within the MIP subsolver of the GBC solver. Has no effect on other solvers.
    - 'debug_mode': If true, debug output is enabled for all solvers. 
"""
function test_HNDPfix(json_file_path)
    
    # Parse the JSON file into a dictionary.
    params = JSON.parsefile(json_file_path)
    
    # Extract required parameters.
    # If any key is missing these will throw an error; you may wish to handle that accordingly.
    users  = params["users"]
    nruns  = params["nruns"]
    betas  = params["betas"]

    # Extract optional parameters with default fallbacks.
    hsolver             = get(params, "hsolver", ["GBC", "BlC", "CA", "CP", "CH"])
    bigMsetting         = get(params, "bigMsetting", ["F", "C"])
    time_limit          = get(params, "time_limit", 3600)
    time_enum_hybrid    = get(params, "time_enum_hybrid", 300)
    partial_decomposition = get(params, "partial_decomposition", false)
    cycle_free_GBC      = get(params, "cycle_free_GBC", false)
    debug_mode          = get(params, "debug_mode", false)

     # Optionally, print the loaded parameters for confirmation.
    println("Executing test_HNDPfix with the following parameters:")
    println("  users: ", users)
    println("  nruns: ", nruns)
    println("  betas: ", betas)
    println("  hsolver: ", hsolver)
    println("  bigMsetting: ", bigMsetting)
    println("  time_limit: ", time_limit)
    println("  time_enum_hybrid: ", time_enum_hybrid)
    println("  partial_decomposition: ", partial_decomposition)
    println("  cycle_free_GBC: ", cycle_free_GBC)
    println("  debug_mode: ", debug_mode)
    
    # start main part of function
    instances = []
    parameters = []
    
    # Generate HNDP multi-layer graph instances
    myfolder = init_logging_folder()
    logger = JuBiC.new_file_logger(myfolder * "/graph_generation_log.txt", true)
    @info "This debuger only contains information on the generation of the HNDP graphs used in the function 'test_HNDPfix'. For loggers of the different model runs, look in the respective folders."
    hndps = with_logger(logger) do
        hndps = Dict(
            (u, nr, be) => build_random_layer_SiouxFalls(u, 0; seed=nr, beta=be) for
            (u, nr, be) in Base.product(users, 1:nruns, betas)
        )
    end

    # build GBC instances
    if "GBC" in hsolver
        myfolderGBC = myfolder * "/GBCSolver"
        create_folder_if_not_exists(myfolderGBC)
        for u in users
            for nr = 1:nruns
                for be in betas
                    # create output folder
                    myfolderrun = myfolderGBC * "/S$(u)_$(nr)_$be"
                    create_folder_if_not_exists(myfolderrun)

                    loggerGBC = JuBiC.new_file_logger(myfolderrun * "/setupGBC$(u)_$(nr)_$be.txt", true)
                    @info "This debuger only contains information on the generation of the GBCSolver instances used in the function 'test_HNDPfix'."
                    with_logger(loggerGBC) do
                        # create instance
                        hndpt = hndps[u, nr, be]
                        inst = to_GBCInstance(
                            hndpt,
                            GurobiSolver(Gurobi.Env());
                            partial_dec=partial_decomposition,
                            MIPsubsolver = true,
                            cycle_free_sub = cycle_free_GBC
                        )
                        gbc_param = GBCparam(
                            GurobiSolver(Gurobi.Env()),
                            debug_mode,
                            myfolderrun,
                            "lp",
                            PARETO_OPTIMALITY_ONLY,
                            time_limit,
                        )

                        # set parameter of instance
                        new_stat!(get_stats(gbc_param), "U", u)
                        new_stat!(get_stats(gbc_param), "partial dec", partial_decomposition)
                        new_stat!(get_stats(gbc_param), "seed", nr)
                        new_stat!(get_stats(gbc_param), "beta", be)

                        # save generated and continue
                        push!(instances, inst)
                        push!(parameters, gbc_param)
                    end
                end
            end
        end
    end #end GBC

    # build BlC
    if "BlC" in hsolver
        myfolderBlC = myfolder * "/BlCSolver"
        create_folder_if_not_exists(myfolderBlC)
        for u in users
            for nr = 1:nruns
                for be in betas
                    for bm in bigMsetting
                        # scip settings taht rely on dual variables
                        if bm != "F" && bm != "C"
                            @info "Sciping bigM setting $bm for BlC solver as it is not directly applicatble for it"
                            continue
                        end
                        # transform big M setting
                        fixBigM = false
                        if bm == "F"
                            fixBigM = true
                        end

                        # create output folder
                        myfolderrun = myfolderBlC * "/S$(u)_$(nr)_$(be)_$bm"
                        create_folder_if_not_exists(myfolderrun)

                        loggerBlC = JuBiC.new_file_logger(myfolderrun * "/setupBlC$(u)_$(nr)_$be.txt", true)
                        @info "This debuger only contains information on the generation of the BlCSolver instances used in the function 'test_HNDPfix'."
                        with_logger(loggerBlC) do
                            # create instance
                            hndpt = hndps[u, nr, be]
                            inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()); MIPsubsolver = true, fixedBigM=fixBigM)
                            @debug "Time limit is set to $time_limit"
                            blc_param = BLCparam(GurobiSolver(Gurobi.Env()), debug_mode, myfolderrun, "lp", time_limit)

                            # set parameter of instance
                            new_stat!(get_stats(blc_param), "U", u)
                            new_stat!(get_stats(blc_param), "seed", nr)
                            new_stat!(get_stats(blc_param), "TypeBigM", bm)
                            new_stat!(get_stats(blc_param), "beta", be)

                            # save generated and continue
                            push!(instances, inst)
                            push!(parameters, blc_param)
                        end
                    end
                end
            end
        end
    end # end BlC

    # build compact arc-based MIP model
    if "CA" in hsolver
        myfolderBlC = myfolder * "/CASolver"
        create_folder_if_not_exists(myfolderBlC)
        for u in users
            for nr = 1:nruns
                for be in betas
                    for bm in bigMsetting
                        # transform big M setting to boolean parameters
                        fixBigM = false
                        if occursin("F", bm)
                            fixBigM = true
                        end
                        
                        indicator = false
                        if occursin("I", bm)
                            indicator = true
                        end

                        boundsL2vars = false
                        if occursin("B", bm)
                            boundsL2vars = true
                        end


                        # create output folder
                        myfolderrun = myfolderBlC * "/S$(u)_$(nr)_$(be)_$bm"
                        create_folder_if_not_exists(myfolderrun)

                        loggerCA = JuBiC.new_file_logger(myfolderrun * "/setupCA$(u)_$(nr)_$(be)_$bm.txt", true)
                        @info "This debuger only contains information on the generation of the arc-based MIP instances used in the function 'test_HNDPfix'."
                        with_logger(loggerCA) do
                            # create instance
                            hndpt = hndps[u, nr, be]
                            inst = to_MIPInstance_arc(hndpt, GurobiSolver(Gurobi.Env()); fixedBigM=fixBigM, indicator=indicator, boundsL2vars=boundsL2vars)
                            ca_param = MIPparam(
                                GurobiSolver(Gurobi.Env()),
                                debug_mode,
                                myfolderrun,
                                "lp",
                                time_limit,
                            )

                            # set parameter of instance
                            new_stat!(get_stats(ca_param), "U", u)
                            new_stat!(get_stats(ca_param), "seed", nr)
                            new_stat!(get_stats(ca_param), "beta", be)
                            new_stat!(get_stats(ca_param), "CompactType", "Arc")
                            new_stat!(get_stats(ca_param), "TypeBigM", bm)


                            # save generated and continue
                            push!(instances, inst)
                            push!(parameters, ca_param)
                        end
                    end
                end
                
            end
        end
    end # end compact arc-based MIP

    # build compact arc-based MIP model
    if "CP" in hsolver
        myfolderCP = myfolder * "/CPSolver"
        create_folder_if_not_exists(myfolderCP)
        enumtime_collection = Dict((u, be) => [] for (u, be) in Base.product(users, betas))
        npaths_collection = Dict((u, be) => [] for (u, be) in Base.product(users, betas))
        for u in users
            for nr = 1:nruns
                for be in betas
                    # create output folder
                    myfolderrun = myfolderCP * "/S$(u)_$(nr)_$be"
                    create_folder_if_not_exists(myfolderrun)

                    loggerCP = JuBiC.new_file_logger(myfolderrun * "/setupCP$(u)_$(nr)_$be.txt", true)
                    @info "This debuger only contains information on the generation of the path-based MIP instances used in the function 'test_HNDPfix'."
                    with_logger(loggerCP) do
                        # create instance
                        hndpt = hndps[u, nr, be]
                        inst, enumtime, sp_enumtimes, enum_enumtimes, npaths = to_MIPInstance_path(hndpt, GurobiSolver(Gurobi.Env()), time_limit)
                        cp_param = MIPparam(
                            GurobiSolver(Gurobi.Env()),
                            debug_mode,
                            myfolderrun,
                            "lp",
                            time_limit-enumtime,
                        )

                        # set parameter of instance
                        new_stat!(get_stats(cp_param), "U", u)
                        new_stat!(get_stats(cp_param), "seed", nr)
                        new_stat!(get_stats(cp_param), "beta", be)
                        new_stat!(get_stats(cp_param), "CompactType", "Path")
                        new_stat!(get_stats(cp_param), "Enumtime", enumtime)
                        new_stat!(get_stats(cp_param), "SumPaths", sum(npaths))
                        new_stat!(get_stats(cp_param), "npath_users", u)
                        
                        for (i, val) in pairs(sp_enumtimes)
                            push!(enumtime_collection[u, be], val + enum_enumtimes[i])
                        end
                        append!(npaths_collection[u, be], npaths)

                        # save generated and continue
                        push!(instances, inst)
                        push!(parameters, cp_param)
                    end
                end

                
            end
        end

        # print collected times for enumeration to file
        enum_to_file(enumtime_collection, myfolderCP * "/enumtimesAggregated.csv")
        enum_to_file(npaths_collection, myfolderCP * "/npathsAggregated.csv")
    end # end compact arc-based MIP

    # build compact arc-based MIP model
    if "CH" in hsolver
        myfolderCH = myfolder * "/CHSolver"
        create_folder_if_not_exists(myfolderCH)
        for u in users
            for nr = 1:nruns
                for be in betas
                    # create output folder
                    myfolderrun = myfolderCH * "/S$(u)_$(nr)_$(be)"
                    create_folder_if_not_exists(myfolderrun)

                    loggerCH = JuBiC.new_file_logger(myfolderrun * "/setupCH$(u)_$(nr).txt", true)
                    @info "This debuger only contains information on the generation of the hybrid MIP instances used in the function 'test_HNDPfix'."
                    with_logger(loggerCH) do
                        # create instance
                        hndpt = hndps[u, nr, be]
                        inst, enumtime, nbad_users = to_MIPInstance_hybrid(hndpt, GurobiSolver(Gurobi.Env()), time_enum_hybrid; fixedBigM=true)
                        ch_param = MIPparam(
                            GurobiSolver(Gurobi.Env()),
                            debug_mode,
                            myfolderrun,
                            "lp",
                            time_limit-enumtime,
                        )

                        # set parameter of instance
                        new_stat!(get_stats(ch_param), "U", u)
                        new_stat!(get_stats(ch_param), "seed", nr)
                        new_stat!(get_stats(ch_param), "beta", be)
                        new_stat!(get_stats(ch_param), "CompactType", "Hybrid")
                        new_stat!(get_stats(ch_param), "TypeBigM", "F")
                        new_stat!(get_stats(ch_param), "Enumtime", enumtime)
                        new_stat!(get_stats(ch_param), "npath_users", nbad_users)

                        # save generated and continue
                        push!(instances, inst)
                        push!(parameters, ch_param)
                    end
                end
            end
        end
    end # end compact arc-based MIP

    # test Sioux
    testrun!(instances, parameters, myfolder)
end




###### Auxiliary functions ########
function enum_to_file(enumtime_collection, filename)
    # Determine the maximum number of values in any of the vectors.
    max_rows = maximum(length.(values(enumtime_collection)))

    # (Optional) sort the keys if you want a defined column order.
    keys_sorted = sort(collect(keys(enumtime_collection)))

    # Create an empty DataFrame.
    df = DataFrame()

    # For each key, create a column.
    # The first row (as header) will be the key and following rows the vector elements.
    for k in keys_sorted
        # For each key, build a column of length `max_rows`.
        # If the current vector is too short, we fill with `missing`.
        column = [ i <= length(enumtime_collection[k]) ? enumtime_collection[k][i] : missing for i in 1:max_rows ]
        # Use the key (converted to string) as the column header.
        df[!, string(k)] = column
    end

    # Write the DataFrame to a CSV file.
    CSV.write(filename, df; delim=";")
end