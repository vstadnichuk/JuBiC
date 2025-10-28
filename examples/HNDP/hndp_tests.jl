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
    test_toy_HNDPwC(hsolver=["GBC", "BlC", "GBCLag", "BlCLag", "BlCLagMiBS"], time_limit=3600, partial_decomposition=false)

Create and solve a simple HNDPwC instance. Mainly intended as test instance which can be verified by hand. 
"""
function test_toy_HNDPwC(hsolver=["GBC", "BlC", "GBCLag", "BlCLag", "BlCLagMiBS"], time_limit=3600, partial_decomposition=false)
    instances = []
    parameters = []

    myfolder = init_logging_folder()
    logger, io = JuBiC.new_file_logger(myfolder * "/debuglogInstances.txt", true)
    @info "This debuger only contains information on the generation of the instances and models used in the function 'test_toy_HNDPwC'. For loggers of the different model runs, look in the respective folders. "
    
    try 
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
                new_stat!(get_stats(gbc_param), "subsolver", SGBC_LABEL)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, gbc_param)
            end

            # build Hierarchical Decomposition model with Lagrangian cuts 
            if "GBCLag" in hsolver
                myfolderGBCLag = myfolder * "/GBCLagSolver"
                create_folder_if_not_exists(myfolderGBCLag)

                # create instance
                subsolvertype = SGBC_MIP_CYCLEFREE 
                inst = to_GBCInstance(
                    hndpt,
                    GurobiSolver(Gurobi.Env());
                    partial_dec=partial_decomposition, 
                    partial_objL2=partial_decomposition,
                    subtype = subsolvertype
                )
                gbclag_param = GBCparam(
                    GurobiSolver(Gurobi.Env()),
                    true,
                    myfolderGBCLag,
                    "lp",
                    PARETO_OPTIMALITY_ONLY,
                    true, # warm start
                    true,  # use Lagrangian cuts for big M
                    true, # trim big M coef
                    time_limit
                )

                # set parameter of instance
                new_stat!(get_stats(gbclag_param), "seed", 42)
                new_stat!(get_stats(gbclag_param), "subsolver", subsolvertype)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, gbclag_param)
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
                new_stat!(get_stats(blc_param), "subsolver", SBlC_MIP)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, blc_param)
            end

            if "BlCLag" in hsolver
                myfolderBlCLag = myfolder * "/BlCLagSolver"
                create_folder_if_not_exists(myfolderBlCLag)

                # create instance
                inst = to_BlCInstance(
                    hndpt,
                    GurobiSolver(Gurobi.Env());
                    subsolver=SBlCLAG_MIP_CYCLEFREE
                )
                blclag_param = BlCLagparam(
                    GurobiSolver(Gurobi.Env()),
                    true,
                    myfolderBlCLag,
                    "lp",
                    PARETO_OPTIMALITY_ONLY,
                    true,
                    time_limit
                )

                # set parameter of instance
                new_stat!(get_stats(blclag_param), "seed", 42)
                new_stat!(get_stats(blclag_param), "subsolver", SBlCLAG_MIP_CYCLEFREE)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, blclag_param)
            end

            if "BlCLagMiBS" in hsolver
                myfolderBlCLagMiBS = myfolder * "/BlCLagMiBSSolver"
                create_folder_if_not_exists(myfolderBlCLagMiBS)

                # create instance
                inst = to_BlCInstance(
                    hndpt,
                    GurobiSolver(Gurobi.Env());
                    subsolver=SBlCLAG_MiBS
                )
                blclagmibs_param = BlCLagparam(
                    GurobiSolver(Gurobi.Env()),
                    true,
                    myfolderBlCLagMiBS,
                    "lp",
                    PARETO_OPTIMALITY_ONLY,
                    true,
                    time_limit
                )

                # set parameter of instance
                new_stat!(get_stats(blclagmibs_param), "seed", 42)
                new_stat!(get_stats(blclagmibs_param), "MIPsubsolver", SBlCLAG_MiBS)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, blclagmibs_param)
            end
        end # end logger
    finally
        close(io)
    end

    # test Sioux toy instance
    testrun!(instances, parameters, myfolder)
end


"""
    test_negative_HNDP(hsolver=["GBC", "GBCLag", "BlC", "BlCLag", "BlCLagMiBS", "CA", "CP", "CH"], time_limit=50, partial_decomposition=false, cycle_free_GBC=true, withweights=false)

Test the algorithm for a toy example of the HNDP with negative cycle according to first-level objective.
# Arguments
    - 'hsolver': The type of solver used for the bilevel problem. Currently supported are "GBC" for GBCSolver, "GBCLag" is the GBCSolver but Lagrian cuts are used for generating big M constants in coef.,
        "BlC" for the BlCSolver, "BlCLag" for BlCLagSolver where Lagrangian cuts are used for BlCuts, "BlCLagMiBS" if the Lagragian subproblem should be solved with MiBS,
        "CA" for compact arc-based model, "CP" for compact path-based model, and "CH" for a hybrid version of the previous two.
        If multiple solvers are passed, solve instance with each. 
    - 'time_limit': Time limit for the solvers.
    - 'partial_decomposition': If true, employ partial decomposition for GBCSolver.
    - 'cycle_free_GBC': If true, additionaly generate Benders-like cuts within the MIP subproblem for GBCSolver. Has no effect on other solvers.
    - 'withweights': If true, the HNDP instance will contain an additional weight knapsack constraint. It will throw an exception if used with a solver who does not support this setting.
"""
function test_negative_HNDP(hsolver=["GBC", "GBCLag", "BlC", "BlCLag", "BlCLagMiBS", "CA", "CP", "CH"], time_limit=30, partial_decomposition=true, cycle_free_GBC=true, withweights=false)
    instances = []
    parameters = []

    myfolder = init_logging_folder()
    logger, io = JuBiC.new_file_logger(myfolder * "/debuglogInstances.txt", true)
    @info "This debuger only contains information on the generation of the instances and models used in the function 'test_negative_HNDP'. For loggers of the different model runs, look in the respective folders. "
    
    try

        with_logger(logger) do
            hndpt = build_toy_HNDPneg()

            # build Hierarchical Decomposition model
            if "GBC" in hsolver
                myfolderGBC = myfolder * "/GBCSolver"
                create_folder_if_not_exists(myfolderGBC)

                # create instance
                subsolvertype = if cycle_free_GBC SGBC_MIP else SGBC_MIP_CYCLEFREE end
                inst = to_GBCInstance(
                    hndpt,
                    GurobiSolver(Gurobi.Env());
                    partial_dec=partial_decomposition, 
                    subtype=subsolvertype
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
                new_stat!(get_stats(gbc_param), "subsolver", subsolvertype)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, gbc_param)
            end

            # build Hierarchical Decomposition model with Lagrangian cuts 
            if "GBCLag" in hsolver
                myfolderGBCLag = myfolder * "/GBCLagSolver"
                create_folder_if_not_exists(myfolderGBCLag)

                # we need to generate bilevel feasible solutions for the big M generation of BlC cuts
                if cycle_free_GBC
                    # create instance
                    subsolvertype = SGBC_MIP_CYCLEFREE
                    inst = to_GBCInstance(
                        hndpt,
                        GurobiSolver(Gurobi.Env());
                        partial_dec = partial_decomposition, 
                        partial_objL2 = partial_decomposition,
                        subtype = subsolvertype
                    )
                    gbclag_param = GBCparam(
                        GurobiSolver(Gurobi.Env()),
                        true,
                        myfolderGBCLag,
                        "lp",
                        PARETO_OPTIMALITY_ONLY,
                        true, # warm start
                        true,  # use Lagrangian cuts for big M
                        true, # trim big M coef
                        time_limit
                    )

                    # set parameter of instance
                    new_stat!(get_stats(gbclag_param), "seed", 42)
                    new_stat!(get_stats(gbclag_param), "subsolver", subsolvertype)

                    # save generated and continue
                    push!(instances, inst)
                    push!(parameters, gbclag_param)
                else
                    @info "We scip GBCLag solver because we need cycle breaking constraints in the subsolver which is not given by current settings."
                end

                
            end

            # build BlC
            if "BlC" in hsolver
                myfolderBlC = myfolder * "/BlCSolver"
                create_folder_if_not_exists(myfolderBlC)

                # create instance
                inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()); subsolver = SBlC_MIP) # we need MIP subsolver for negative cycles right now
                blc_param =
                    BLCparam(GurobiSolver(Gurobi.Env()), true, myfolderBlC, "lp", time_limit)

                # set parameter of instance
                new_stat!(get_stats(blc_param), "seed", 42)
                new_stat!(get_stats(blc_param), "subsolver", SBlC_MIP)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, blc_param)
            end

            # build BlCLag
            if "BlCLag" in hsolver
                myfolderBlCLag = myfolder * "/BlCLagSolver"
                create_folder_if_not_exists(myfolderBlCLag)

                # create instance
                subsolvertype = if cycle_free_GBC SBlCLAG_MIP_CYCLEFREE else SGBC_MiBS end
                inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()); subsolver=subsolvertype) 
                blclag_param =
                    BlCLagparam(GurobiSolver(Gurobi.Env()), true, myfolderBlCLag, "lp", time_limit)

                # set parameter of instance
                new_stat!(get_stats(blclag_param), "seed", 42)
                new_stat!(get_stats(blclag_param), "subsolver", subsolvertype)

                # save generated and continue
                push!(instances, inst)
                push!(parameters, blclag_param)
            end

            # build arc compact solver instance
            if "CA" in hsolver
                if withweights
                    throw(ArgumentError("Solver CA does not support weights on the second level as then the second level is non-convex."))
                end

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
                if withweights
                    throw(ArgumentError("Solver CP does not support weights on the second level as the path enumeration algorithm was not tested for thissetting yet."))
                end

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
                if withweights
                    throw(ArgumentError("Solver CH does not support weights on the second level as then the second level is non-convex."))
                end

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
    finally
        close(io)
    end

    # test Sioux toy instance
    testrun!(instances, parameters, myfolder)
end



"""
    test_HNDPwC(users, alphas, nruns::Number; hsolver=["GBC", "BlC"], time_limit=3600, partial_decomposition=false)

Basic test for the Sioux Falls network Hazmat Network Design Problem with Capacity Constraints. 
You should path a .json file with the following parameter lists. All combinations of the parameter contained in the lsits will be tested. 
- 'fixnetwork': A Bool. If true, use instances with a fixed network and negative arc cost
- 'users': A list with the number of users for the individual runs.
- 'alphas': The different alpha values that may appear. The alpha should be a parameter between 0 and 1, and sets the capacity limit for the test instance. 0 indicates very tight bounds, and 1 that there is no capacity limit. 
- 'nruns': Number of runs with different seed for each instance.
- 'hsolver': The type of solver used for the bilevel problem. Currently supported are "GBC" for GBCSolver, "BlC" for the BlCSolver, and "BlCLag" for the BlCLagSolver.  
If multiple solvers are passed, solve instance with each. 
- 'time_limit': The time limit set on (each) solver (in sec).
- 'partial_decomposition': A list containing true or false (or both). If true, apply partial decomposition where applicable, i.e., when using our GBCSolver
- 'two_stage': If true, use same objective function for risk and cost, i.e., both levels cooperate to find cost-minimal paths 
- 'constrac_cost': The cost for including an arc are generated randomly within 0:constructioncost
- 'warmstartGBC': A list containing true or false (or both). If false, ConnectorLP is reset after each iteration (what is very stupid, so only use this for test purpose)
- 'paretooptcuts': A list containing true or false (or both). If true, use Pareto optimal cuts are used. Note that this applies to all cuts that are generated by solving Lagrangian dual.
- 'subsolvertype': This is a list of solver types. The currently supported types are "MIP" for a simple MIP solver, "MIP_CYCLE" for feasible solutions only, and "MiBS" for MiBS. Note that this setting only affects GBC and BlC solvers.
- 'LCinGBC': A list containing true or false (or both). If true, when solving ConnectorLP for GBCSolver, use a second ConnectorLP to compute the big M coefficients needed for Bilevel Lagrangian cuts.  
- 'blCandGBCcuts': A list containing true or false (or both). If true, when generating GBC, also add BlC to master if coef. were computed in GBC procedure (i.e., 'LCinGBC' setting was true).
"""
function test_HNDPwC(json_path::String)
    # parse passed json file
    cfg = JSON.parsefile(json_path)

    fixnetwork            = get(cfg, "fixnetwork", true)
    users                 = cfg["users"]
    alphas                = cfg["alphas"]
    beta                  = get(cfg, "beta", 1)
    nruns                 = cfg["nruns"]
    hsolver               = get(cfg, "hsolver", ["GBC", "BlC", "BlCLag"])
    time_limit            = get(cfg, "time_limit", 3600)
    partial_decomposition = get(cfg, "partial_decomposition", [true])
    two_stage             = get(cfg, "two_stage", false)
    constrac_cost         = get(cfg, "constrac_cost", 0)
    warmstartGBC          = get(cfg, "warmstartGBC", [true])
    paretooptcuts         = get(cfg, "paretooptcuts", [true])
    subsolvertype         = get(cfg, "subsolvertype", ["MIP", "MIP_CYCLE", "MiBS"])
    LCinGBC               = get(cfg, "LCinGBC", [false])
    blCandGBCcuts         = get(cfg, "blCandGBCcuts", [false])
    trim_coeff_opt        = get(cfg, "trim_coeff", [true])
    debug_mode            = get(cfg, "debug_mode", false)


    # init folder and lists
    myfolder = init_logging_folder()
    instances = []
    parameters = []
    pareto_options = vcat(
    (true  in paretooptcuts ? [PARETO_OPTIMALITY_AND_FEASIBILITY] : []),
    (false in paretooptcuts ? [PARETO_NONE]                   : []) )

    # save used settings in folder
    open(joinpath(myfolder, "params.json"), "w") do io
        JSON.print(io, cfg)
    end

    # generate the underlying HNDPwC instances
    if fixnetwork
        hndps = Dict(
            (u, al, nr) => build_random_layer_SiouxFalls(u, al; seed=nr, beta=beta, withweight=true) for
            (u, al, nr) in Base.product(users, alphas, 1:nruns)
        )
    else
        hndps = Dict(
            (u, al, nr) => build_random_SiouxFalls(u, al; seed=nr, two_stage=two_stage, constructioncost=constrac_cost) for
            (u, al, nr) in Base.product(users, alphas, 1:nruns)
        )
    end

    # build GBC instances
    if "GBC" in hsolver
        myfolderGBC = myfolder * "/GBCSolver"
        create_folder_if_not_exists(myfolderGBC)
        for (u, al, nr, partdec, stabopt, wstart, st, LCsub, blc_gbc, trim) in Base.product(users, alphas, 1:nruns, partial_decomposition, pareto_options, warmstartGBC, subsolvertype, LCinGBC, blCandGBCcuts, trim_coeff_opt)
            # check conditions for parameter combinations
            if LCsub 
                if st != "MIP_CYCLE" && st != "MIP_CYCLE"
                    @warn "The chosen GBC subsolver does not guarantee bilevel feasible points. Therefore, discarding the solver setting $st in combination with generating big M with Lagrangian dual."
                    continue 
                end
            end
            if blc_gbc
                if !partdec || !LCsub
                    @warn "We require partial decomposition and automated generation of big M coefficients for the option to add BlC to master in GBCSolver. As partdec=$partdec and LCsub=$LCsub, this is not given, and we scip current setting."
                    continue 
                end
            end
            
            # create output folder
            myfolderrun = myfolderGBC * "/S$(u)_$(al)_$(nr)_P$(partdec)_$(stabopt)_W$(wstart)_$(st)_LC$(LCsub)_BlC$(blc_gbc)_T$(trim)"
            create_folder_if_not_exists(myfolderrun)

            # find out solver type
            if st == "MIP"
                solvertype = SGBC_MIP
            elseif st == "MIP_CYCLE"
                solvertype = SGBC_MIP_CYCLEFREE
            elseif st == "MiBS"
                solvertype = SGBC_MiBS
            else
                throw(ArgumentError("The passed subsolver type $st unknown!"))
            end
            
            # create instance
            hndpt = hndps[u, al, nr]
            inst = to_GBCInstance(
                hndpt,
                GurobiSolver(Gurobi.Env());
                partial_dec=partdec,
                partial_objL2=blc_gbc,
                subtype = solvertype
            )
            gbc_param = GBCparam(
                GurobiSolver(Gurobi.Env()),
                debug_mode,
                myfolderrun,
                "lp",
                stabopt,
                wstart,
                LCsub,
                trim,
                time_limit,
            )

            # set parameter of instance
            new_stat!(get_stats(gbc_param), "fixnetwork", fixnetwork)
            new_stat!(get_stats(gbc_param), "U", u)
            new_stat!(get_stats(gbc_param), "alpha", al)
            if fixnetwork new_stat!(get_stats(gbc_param), "beta", beta) end
            new_stat!(get_stats(gbc_param), "constructioncost", constrac_cost)
            new_stat!(get_stats(gbc_param), "partial_decomposition", partdec)
            new_stat!(get_stats(gbc_param), "seed", nr)
            new_stat!(get_stats(gbc_param), "mip_subsolver", st)
            new_stat!(get_stats(gbc_param), "stabopt", stabopt)
            new_stat!(get_stats(gbc_param), "warmstart", wstart)
            new_stat!(get_stats(gbc_param), "solvertype", solvertype)
            new_stat!(get_stats(gbc_param), "debug_mode", debug_mode)
            new_stat!(get_stats(gbc_param), "bigMwithLC", LCsub)
            new_stat!(get_stats(gbc_param), "blCandGBCcuts", blc_gbc)
            new_stat!(get_stats(gbc_param), "trim_coef", trim)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, gbc_param)
        end
    end

    # build BlC
    if "BlC" in hsolver
        myfolderBlC = myfolder * "/BlCSolver"
        create_folder_if_not_exists(myfolderBlC)
        for (u, al, nr) in Base.product(users, alphas, 1:nruns)
            # create output folder
            myfolderrun = myfolderBlC * "/S$(u)_$(al)_$(nr)"
            create_folder_if_not_exists(myfolderrun)

            # create instance
            hndpt = hndps[u, al, nr]
            inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()); subsolver=SBlC_MIP, fixedBigM=false)
            blc_param = BLCparam(GurobiSolver(Gurobi.Env()), debug_mode, myfolderrun, "lp", time_limit)

            # set parameter of instance
            new_stat!(get_stats(blc_param), "fixnetwork", fixnetwork)
            new_stat!(get_stats(blc_param), "U", u)
            new_stat!(get_stats(blc_param), "alpha", al)
            if fixnetwork new_stat!(get_stats(blc_param), "beta", beta) end
            new_stat!(get_stats(blc_param), "constructioncost", constrac_cost)
            new_stat!(get_stats(blc_param), "seed", nr)
            new_stat!(get_stats(blc_param), "mip_subsolver", SBlC_MIP)
            new_stat!(get_stats(blc_param), "debug_mode", debug_mode)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, blc_param)
        end
    end

    # build BlCLag instances
    if "BlCLag" in hsolver
        myfolderGBC = myfolder * "/BlCLagSolver"
        create_folder_if_not_exists(myfolderGBC)
        for (u, al, nr, subsolver, stabopt, wstart) in Base.product(users, alphas, 1:nruns, subsolvertype, pareto_options, warmstartGBC)
            # create output folder
            myfolderrun = myfolderGBC * "/S$(u)_$(al)_$(nr)_$(subsolver)_$(stabopt)_W$(wstart)"
            create_folder_if_not_exists(myfolderrun)

            # find out solver type
            if subsolver == "MIP"
                @warn "The chosen BlCLag subsolver does not guarantee bilevel feasible points. Therefore, discarding the solver setting $subsolver in combination with generating big M with Lagrangian dual."
                continue # we require bilevel feasible subproblem
            elseif subsolver == "MIP_CYCLE"
                subsolvertype = SBlCLAG_MIP_CYCLEFREE
            elseif subsolver == "MiBS"
                subsolvertype = SBlCLAG_MiBS
            else
                throw(ArgumentError("unsupported subsolver type $(subsolver)"))
            end

            # create instance
            hndpt = hndps[u, al, nr]
            inst = to_BlCInstance(
                hndpt, 
                GurobiSolver(Gurobi.Env()); 
                subsolver = subsolvertype
            )
            blclag_param = BlCLagparam(
                GurobiSolver(Gurobi.Env()),
                debug_mode,
                myfolderrun,
                "lp",
                stabopt,
                wstart,
                time_limit,
            )

            # set parameter of instance
            new_stat!(get_stats(blclag_param), "fixnetwork", fixnetwork)
            new_stat!(get_stats(blclag_param), "U", u)
            new_stat!(get_stats(blclag_param), "alpha", al)
            if fixnetwork new_stat!(get_stats(blclag_param), "beta", beta) end
            new_stat!(get_stats(blclag_param), "constructioncost", constrac_cost)
            new_stat!(get_stats(blclag_param), "seed", nr)
            new_stat!(get_stats(blclag_param), "stabopt", stabopt)
            new_stat!(get_stats(blclag_param), "warmstart", wstart)
            new_stat!(get_stats(blclag_param), "subsolver", subsolvertype)
            new_stat!(get_stats(blclag_param), "debug_mode", debug_mode)

            # save generated and continue
            push!(instances, inst)
            push!(parameters, blclag_param)
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
    - 'betas': A list of beta values employed when generating instances 

    - 'hsolver': The type of solver used for the bilevel problem. Currently supported are "GBC" for GBCSolver, "BlC" for the BlCSolver, "BlCLag" for BlCLagSolver
    where Benders-like cut coeff. are computed by solving a Lagrangian dual,
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
    - 'withweights': If true, add weights to second level. Note that this results in non-convex subproblem. As some solvers do not support them, it will throw an exception if used in combination.
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
    hsolver             = get(params, "hsolver", ["GBC", "BlC", "BlCLag", "CA", "CP", "CH"])
    bigMsetting         = get(params, "bigMsetting", ["F", "C"])
    time_limit          = get(params, "time_limit", 3600)
    time_enum_hybrid    = get(params, "time_enum_hybrid", 300)
    partial_decomposition = get(params, "partial_decomposition", false)
    cycle_free_GBC      = get(params, "cycle_free_GBC", false)
    debug_mode          = get(params, "debug_mode", false)
    withweights          = get(params, "withweights", false)

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
    println("  withweights: ", withweights)
    
    # start main part of function
    instances = []
    parameters = []
    
    # Generate HNDP multi-layer graph instances
    myfolder = init_logging_folder()
    logger, io = JuBiC.new_file_logger(myfolder * "/graph_generation_log.txt", true)

    try
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

                        loggerGBC, io = JuBiC.new_file_logger(myfolderrun * "/setupGBC$(u)_$(nr)_$be.txt", true)
                        try 
                            @info "This debuger only contains information on the generation of the GBCSolver instances used in the function 'test_HNDPfix'."
                            with_logger(loggerGBC) do
                                # create instance
                                hndpt = hndps[u, nr, be]
                                subsolvertype = if cycle_free_GBC SGBC_MIP_CYCLEFREE else SGBC_MIP end
                                inst = to_GBCInstance(
                                    hndpt,
                                    GurobiSolver(Gurobi.Env());
                                    partial_dec=partial_decomposition,
                                    subtype = subsolvertype
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
                                new_stat!(get_stats(gbc_param), "subsolver", subsolvertype)

                                # save generated and continue
                                push!(instances, inst)
                                push!(parameters, gbc_param)
                            end
                        finally
                            close(io)
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

                            loggerBlC, io = JuBiC.new_file_logger(myfolderrun * "/setupBlC$(u)_$(nr)_$be.txt", true)
                            try 
                                @info "This debuger only contains information on the generation of the BlCSolver instances used in the function 'test_HNDPfix'."
                                with_logger(loggerBlC) do
                                    # create instance
                                    hndpt = hndps[u, nr, be]
                                    inst = to_BlCInstance(hndpt, GurobiSolver(Gurobi.Env()); subsolver = SBlC_MIP)
                                    @debug "Time limit is set to $time_limit"
                                    blc_param = BLCparam(GurobiSolver(Gurobi.Env()), debug_mode, myfolderrun, "lp", time_limit)

                                    # set parameter of instance
                                    new_stat!(get_stats(blc_param), "U", u)
                                    new_stat!(get_stats(blc_param), "seed", nr)
                                    new_stat!(get_stats(blc_param), "TypeBigM", bm)
                                    new_stat!(get_stats(blc_param), "beta", be)
                                    new_stat!(get_stats(blc_param), "subsolver", SBlC_MIP)


                                    # save generated and continue
                                    push!(instances, inst)
                                    push!(parameters, blc_param)
                                end
                            finally
                                close(io)
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
                            if withweights
                                throw(ArgumentError("Solver CA does not support weights on the second level as then the second level is non-convex."))
                            end

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

                            loggerCA, io = JuBiC.new_file_logger(myfolderrun * "/setupCA$(u)_$(nr)_$(be)_$bm.txt", true)
                            try 
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
                            finally
                                close(io)
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
                        if withweights
                            throw(ArgumentError("Solver CP does not support weights on the second level as the poath enumeration was not tested for this setting yet."))
                        end

                        # create output folder
                        myfolderrun = myfolderCP * "/S$(u)_$(nr)_$be"
                        create_folder_if_not_exists(myfolderrun)

                        loggerCP, io = JuBiC.new_file_logger(myfolderrun * "/setupCP$(u)_$(nr)_$be.txt", true)
                        try 
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
                        finally
                            close(io)
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
                        if withweights
                            throw(ArgumentError("Solver CH does not support weights on the second level as then the second level is non-convex."))
                        end

                        # create output folder
                        myfolderrun = myfolderCH * "/S$(u)_$(nr)_$(be)"
                        create_folder_if_not_exists(myfolderrun)

                        loggerCH, io = JuBiC.new_file_logger(myfolderrun * "/setupCH$(u)_$(nr).txt", true)
                        try
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
                        finally 
                            close(io)
                        end
                    end
                end
            end
        end # end compact arc-based MIP
    finally
        close(io)
    end

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