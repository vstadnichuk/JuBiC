# Transfer HNDP instances into Instance objects for different solvers
using JuBiC
using JuMP
include("hndp_instances.jl")
include("hdnp_astar_wrapper.jl")


########## Build an Instance for the GBCSolver ##########
"""
    to_GBCInstance(hndp::HNDPwC, solver::SolverWrapper; partial_dec=true, MIPsubsolver=false)

Generate a new Hierarchical Decomposition model from the passed HNDPwC instance. 
# Arguments:
    - 'hndp::HNDPwC': The instance of the HNDPwC.
    - 'solver::SolverWrapper': The MIP solver to use, e.g., Gurobi
    - 'partial_dec=true': If true, employs partial decomposition approach by adding flow constraints to master
    - 'MIPsubsolver=false': If true, use MIP formulation for subproblems instead of Labeling solver. Use this in case first-level objective contains negative values.
    - 'cycle_free_sub=false': If true, additionally add Benders-like cuts to subproblem to ensure that found solutions are cycle free.  
"""
function to_GBCInstance(hndp::HNDPwC, solver::SolverWrapper; partial_dec=true, MIPsubsolver=false, cycle_free_sub=false)
    @info "Starting generation of GBCSolver for passed HNDP instance" 
    A = [(src(e), dst(e)) for e in edges(hndp.mygraph)]
    unames = [user.uname for user in hndp.users]

    # build master
    # master MIP
    mm = Model(() -> get_next_optimizer(solver))

    # create arc binary variables and fix non-decision to 1
    @variable(mm, x[A], Bin)
    for et in A
        if !(et in hndp.edgeA)
            @constraint(mm, x[et] == 1, base_name = "fix$(et))")
        end
    end

    # add variable that represents the construction cost (easier to debug with it)
    constructioncost = sum(hndp.edge_price[a]*x[a] for a in A)
    @variable(mm, construction_cost_var)
    @constraint(mm, constructioncost == construction_cost_var, base_name="constructioncost")

    # create master (obj. has no impact from x-variables)
    @objective(mm, Min, construction_cost_var)  # only first-level cost are the construction cost for the arcs
    xdec_dict = Dict(a => x[a] for a in hndp.edgeA) # Only decision arcs are relevant as all other interdiction constraints are satisfied by default

    

    # if partial decomposition should be applied, create here the function adding L2 vars and constraints
    if partial_dec
        partial(mmodel::JuMP.Model, mobj) = begin
            for user in hndp.users
                oL1, oL2, f = make_user_constraints(user, mmodel, hndp, A, x, false; l2binary=false)  # no need for L2 binary variables in partial decomposition
                @constraint(
                    mmodel,
                    oL1 == mobj[user.uname],
                    base_name = "objLink$(user.uname)_"
                )
            end
        end
        #master = Master(mm, A, xdec_dict, unames, partial)
        # Only decision arcs are relevant as all other interdiction constraints are satisfied by default
        master = Master(mm, hndp.edgeA, xdec_dict, unames, partial)
    else
        #master = Master(mm, A, xdec_dict, unames)
        # Only decision arcs are relevant as all other interdiction constraints are satisfied by default
        master = Master(mm, hndp.edgeA, xdec_dict, unames)
    end

    # build subproblems with A*-search solvers 
    usolvers = []
    for user in hndp.users
        if !MIPsubsolver
            @info "Adding A*search subsolver for user $(user.uname)"
            #solver_user = build_Astar_user(user, hndp, A)
            # Only decision arcs are relevant as all other interdiction constraints are satisfied by default
            solver_user = build_Astar_user(user, hndp, hndp.edgeA)
        else
            @info "Adding MIP subsolver for user $(user.uname)"
            #solver_user = build_MIP_user(user, hndp, A, solver)
            # Only decision arcs are relevant as all other interdiction constraints are satisfied by default
            solver_user = build_MIP_user(user, hndp, A, hndp.edgeA, solver; cycle_free=cycle_free_sub)
        end
        push!(usolvers, solver_user)
    end

    # build instance
    @info "Finished generation of GBCSolver for passed HNDP instance" 
    rV = Instance(master, usolvers)
    return rV
end


########## Build an Instance for the BlCSolver ##########
"""
    to_BlCInstance(hndp::HNDPwC, solver::SolverWrapper; MIPsubsolver = false, fixedbigM=true)

Generate a new Benders-like Decomposition model from the passed HNDPwC instance. 
# Arguments:
    - 'hndp::HNDPwC': The HNDPwC instance.
    - 'solver::SolverWrapper': The MIP solver to use, e.g., Gurobi
    - 'partial_dec=true': If true, employs partial decomposition approach by adding flow constraints to master
    - 'MIPsubsolver=false': If true, use MIP formulation for subproblems instead of Labeling solver. Use this in case first-level objective contains negative values.
    - 'fixedBigM': If 'true', try to generate smart big M based on a path in fixed network. If not possible, use default big M. 
"""
function to_BlCInstance(hndp::HNDPwC, solver::SolverWrapper; MIPsubsolver=false, fixedBigM=true)
    @info "Starting generation of BlCSolver for passed HNDP instance" 
    A = [(src(e), dst(e)) for e in edges(hndp.mygraph)]
    unames = [user.uname for user in hndp.users]

    # build high-point relaxation
    hpr = Model(() -> get_next_optimizer(solver))

    # create arc binary variables and fix non-decision to 1
    @variable(hpr, x[A], Bin)
    for et in A
        if !(et in hndp.edgeA)
            @constraint(hpr, x[et] == 1, base_name = "fix$(et))")
        end
    end

    # create sub_problem constraints and variables
    master_objs = Dict()
    sub_objs = Dict()
    for user in hndp.users
        oL1, oL2, f = make_user_constraints(user, hpr, hndp, A, x, false)

        # a small hack to see easily second level solution later in output
        optL2 = @variable(hpr, base_name = "optL2$(user.uname)")
        @constraint(hpr, optL2 == oL2, base_name = "objUser$(user.uname)")

        # set parameters for master
        master_objs[user.uname] = oL1
        sub_objs[user.uname] = optL2
    end

    # add variable that represents the construction cost (easier to debug with it)
    constructioncost = sum(hndp.edge_price[a]*x[a] for a in A)
    @variable(hpr, construction_cost_var)
    @constraint(hpr, constructioncost == construction_cost_var, base_name="constructioncost")

    # create Benders-like cuts master 
    @objective(hpr, Min, construction_cost_var + sum(values(master_objs)))
    xdec_dict = Dict(a => x[a] for a in hndp.edgeA)
    master = BlCMaster(hpr, hndp.edgeA, xdec_dict, derive_bigM_function(hndp, A, fixedBigM), unames, sub_objs)

    # build subproblems with A*-search solvers 
    usolvers = []
    for user in hndp.users
        if !MIPsubsolver
            @info "Adding A*search subsolver for user $(user.uname)"
            solver_user = build_Astar_user(user, hndp, A)
        else
            @info "Adding MIP subsolver for user $(user.uname)"
            solver_user = build_MIP_user(user, hndp, A, hndp.edgeA, solver)
        end
        push!(usolvers, solver_user)
    end

    # build instance
    @info "Finished generation of BlCSolver for passed HNDP instance" 
    rV = Instance(master, usolvers)
    return rV
end


########## Build an Instance of the compact strong duality model with arc variables ##########
"""
    to_MIPInstance_arc(hndp::HNDPwC, solver::SolverWrapper; fixedBigM=true)

Build the arc-based compact model for the passed HNDP istance. It should not contain weight parameters, as then the reformulation is no longer guaranteed to be correct. 
However, it is up to the user to check this, as there might be cases where you still want to apply the reformulation. Note that there will not be any dual constrains added for any weight related constraints. Arguments
# Arguments
- 'hndp::HNDPwC': The problem instance.
- 'solver::SolverWrapper: The solver to use for the MIP model

# Optional
- 'fixedBigM=true': If 'true', try to generate smart big M based on shortest path in fixed network first. 
- 'indicator=true': If true, use indicator constraints instead of classical big M constraints.
- 'boundsL2vars=true': If true, add additional bounds on second-level variables based on big M values. 
"""
function to_MIPInstance_arc(hndp::HNDPwC, solver::SolverWrapper; fixedBigM=true, indicator=true, boundsL2vars=true)
    @info "Starting generation of compact arc-based MIP model for passed HNDP instance" 
    A = [(src(e), dst(e)) for e in edges(hndp.mygraph)]

    # build MIP model
    mip = Model(() -> get_next_optimizer(solver))

    # create arc binary variables and fix non-decision to 1
    @variable(mip, x[A], Bin)
    for et in A
        if !(et in hndp.edgeA)
            @constraint(mip, x[et] == 1, base_name = "fix$(et))")
        end
    end

    # create sub_problem constraints and variables
    master_objs = Dict()
    sub_objs = Dict()
    for user in hndp.users
        oL1, oL2, f = make_user_constraints(user, mip, hndp, A, x, true; bigM_function=derive_bigM_function(hndp, A, fixedBigM), indicatorconst=indicator, with_var_bounds=boundsL2vars)

        # a small hack to see easily second level solution later in output
        optL2 = @variable(mip, base_name = "optL2$(user.uname)")
        @constraint(mip, optL2 == oL2, base_name = "objUser$(user.uname)")

        # set mip model objective terms
        master_objs[user.uname] = oL1
        sub_objs[user.uname] = optL2
    end

    # add variable that represents the construction cost (easier to debug with it)
    constructioncost = sum(hndp.edge_price[a]*x[a] for a in A)
    @variable(mip, construction_cost_var)
    @constraint(mip, constructioncost == construction_cost_var, base_name="constructioncost")

    # set objective 
    @objective(mip, Min, construction_cost_var + sum(values(master_objs)))

   # build instance
    @info "Finished generation of compact arc-based MIP model for passed HNDP instance" 
    master = MIPMaster(mip)
    rV = Instance(master, nothing)
    return rV
end



########## Build an Instance of the compact strong duality model with path variables ##########
"""
    to_MIPInstance_path(hndp::HNDPwC, solver::SolverWrapper, timelimit)

Build the compact path-based model. Note that the path enumeration conducted can take into accounts weight bounds within the path. 
It is assumed that each user has a feasible path in the fixed network, as otherwise path-based enumeration model reduce to full enumeration. An exception will be thrown if this assumption is violated. 
# Arguments
- 'hndp::HNDPwC': The problem instance.
- 'solver::SolverWrapper: The solver to use for the MIP model
- 'timelimit': The timelimit in seconds

# Returns
The model and also the remaining time for solving taking into account that path enumeration consumed runtime.
For the calculation of the consumed runtime, we assume that the process can be fully parallelized. Therefore, only the runtime of the user with the highest enumeration time is deduced from the original runtime.
    - 'rV': The model
    - 'enummaxtime': The overall runtime required for enumeration. If one user paths could not be enumerated within the time limit, we return 'timelimit'.
    - 'lruntimes_sp': The runtimes messured for computing the shortest path in the fixed network. In case you need it for statistics.
    - 'lruntimes_enum': The runtimes messured for enumeration. In case you need it for statistics.
    - 'npaths': A list of the number of paths enumerated for each user 
"""
function to_MIPInstance_path(hndp::HNDPwC, solver::SolverWrapper, timelimit)
    @info "Starting generation of compact path-based MIP model for passed HNDP instance" 
    A = [(src(e), dst(e)) for e in edges(hndp.mygraph)]

    # build MIP model
    mip = Model(() -> get_next_optimizer(solver))

    # create arc binary variables and fix non-decision to 1
    @variable(mip, x[A], Bin)
    for et in A
        if !(et in hndp.edgeA)
            @constraint(mip, x[et] == 1, base_name = "fix$(et))")
        end
    end

    # enumerate all paths
    upaths = Dict()
    lruntimes_sp = []
    lruntimes_enum = []
    npaths = [] # number of paths enumerated for the user (-1 if no paths enumerated)
    mintimeleft = timelimit # the runtime left if we could fully parallelize the process of enumerating paths 
    for user in hndp.users
        time_left = timelimit

        # compute shortest path in fixed network
        @debug "Started computing shortest path in fixed network for user $(user.uname) with time limit $timelimit"
        shortest_fixed_path, runtime_sp = shortest_unforbiddable_path(user, hndp, time_left)
        if runtime_sp < 0
            # no path found, i.e., no fixed network path for this user
            throw(ArgumentError("User $(user.uname) does not posses a path in the fixed network!"))
        end
        time_left = time_left - runtime_sp
        push!(lruntimes_sp, runtime_sp)
        if time_left <= 0
            # stop if runtime is used up
            @debug "Timeout while computing shortest path in fixed network for user $(user.uname). The timelimit was $timelimit"
            mintimeleft = 0
            push!(lruntimes_enum, 0)  # no time spend on enumeration
            push!(npaths, -1)  # not all paths found
            break
        end
        lengthbound = shortest_fixed_path.cost
        @debug "Computed shortest path in fixed network for user $(user.uname) with cost $lengthbound"

        # enumerate all paths
        @debug "Started enumerating paths in fixed network for user $(user.uname)"
        paths, runtime_enum = enumerate_user_paths(user, hndp, lengthbound, time_left)
        time_left = time_left - runtime_enum
        push!(lruntimes_enum, runtime_enum)
        if time_left <= 0
            # stop if runtime is used up
            @debug "Timeout while enumerating paths in fixed network for user $(user.uname). The timelimit was $timelimit"
            mintimeleft = 0
            push!(npaths, -1)  # not all paths found
            break
        end
        upaths[user] = paths
        push!(npaths, length(paths))  # save number of paths found

        # update 'mintimeleft' to be the minimal amount of runtime left after path enumeration
        @debug "Enumeration finished for user $(user.uname). We found $(length(paths)) paths for user $(user.uname)"
        if time_left < mintimeleft
            if time_left < 0
                throw(ErrorException("While enumerating paths, we used up strictly more runtime that the timelimit. The remaining runtime would be $(mintimeleft)!"))
            end
            mintimeleft = time_left
        end
    end

    # if time is up, we do not build any second level constraints
    if mintimeleft <= 0
        @warn "Enumeration of Path could not be concluded within time limit. Now returning empty model. " 
        master = MIPMaster(mip)
        rV = Instance(master, nothing)
        return rV, timelimit, lruntimes_sp, lruntimes_enum, npaths
    end

    # enumeration was successfull. Now generate constraints 
    @info "Enumeration of Paths was successfull. Worst-case runtime was $(timelimit - mintimeleft)." 
    master_objs = Dict()
    sub_objs = Dict()
    for user in hndp.users
        oL1, oL2 = make_user_constraints_path(user, mip, hndp, upaths[user], A, x, true; l2binary=true)

        # a small hack to see easily second level solution later in output
        optL2 = @variable(mip, base_name = "optL2$(user.uname)")
        @constraint(mip, optL2 == oL2, base_name = "objUser$(user.uname)")

        # set mip objective trms
        master_objs[user.uname] = oL1
        sub_objs[user.uname] = optL2
    end

    # add variable that represents the construction cost (easier to debug with it)
    constructioncost = sum(hndp.edge_price[a]*x[a] for a in A)
    @variable(mip, construction_cost_var)
    @constraint(mip, constructioncost == construction_cost_var, base_name="constructioncost")

    # set objective 
    @objective(mip, Min, construction_cost_var + sum(values(master_objs)))

   # build instance
    @info "Finished generation of compact path-based MIP model for passed HNDP instance" 
    master = MIPMaster(mip)
    rV = Instance(master, nothing)
    return rV, timelimit - mintimeleft, lruntimes_sp, lruntimes_enum, npaths
end

########## Build an Instance of the compact strong duality model with hybrid model for arc and path based constraints ##########
"""
    to_MIPInstance_hybrid(hndp::HNDPwC, solver::SolverWrapper, timelimit_pathenumeration; fixedBigM=true)

Build the compact MIP model for the HNDP. This hybrid MIP combines the path and arc formulation. Up to the passed 'timelimit_pathenumeration', it tries to enumerate all paths of a user. 
If succesfull, the path-based formulation is used for this user. Otherwise, the arc-based formulation is employed for this user. 
It is assumed that each user has a feasible path in the fixed network, as otherwise path-based enumeration model reduce to full enumeration. An exception will be thrown if this assumption is violated. 
# Arguments
    - 'hndp::HNDPwC': The problem instance.
    - 'solver::SolverWrapper: The solver to use for the MIP model
    - 'timelimit': The timelimit in seconds
    - 'timelimit_pathenumeration': The timelimit for the enumeration of paths. If path could not be enuemrated in this time for a user, the arc-based formulation is used for this specific user. 
    - 'fixedBigM': If 'true', try to generate smart big M based on shortest path in fixed network first. 

# Returns
The model and also the remaining time for solving taking into account that path enumeration consumed runtime.
For the calculation of the consumed runtime, we assume that the process can be fully parallelized. Therefore, only the runtime of the user with the highest enumeration time is deduced from the original runtime.
- 'rV': The model
- 'enummaxtime': The overall runtime it required to enumerate the paths. Note that if one user required maximum available time for enumeration, 'timelimit_pathenumeration' is returned.
"""
function to_MIPInstance_hybrid(hndp::HNDPwC, solver::SolverWrapper, timelimit_pathenumeration; fixedBigM=true)
    @info "Starting generation of compact hybrid MIP model that combines path and arc model ideas for passed HNDP instance" 
    A = [(src(e), dst(e)) for e in edges(hndp.mygraph)]

    # build MIP model
    mip = Model(() -> get_next_optimizer(solver))

    # create arc binary variables and fix non-decision to 1
    @variable(mip, x[A], Bin)
    for et in A
        if !(et in hndp.edgeA)
            @constraint(mip, x[et] == 1, base_name = "fix$(et))")
        end
    end

    # enumerate all paths
    upaths = Dict()
    bad_users = [] # list of users we could not enumerate all paths for
    enummaxtime = 0 # the runtime left if we could fully parallelize the process of enumerating paths 
    for user in hndp.users
        time_left = timelimit_pathenumeration

        # compute shortest path in fixed network
        @debug "Started computing shortest path in fixed network within hybrid model for user $(user.uname) with time limit $timelimit_pathenumeration"
        shortest_fixed_path, runtime_sp = shortest_unforbiddable_path(user, hndp, time_left)
        if runtime_sp < 0
            # no path found, i.e., no fixed network path for this user
            throw(ArgumentError("User $(user.uname) does not posses a path in the fixed network! This case could in theory be handled in hybrid model but is not supported currently."))
        end
        time_left = time_left - runtime_sp
        if time_left <= 0
            # we used up all enumeration time for this user
            @debug "Timeout while computing shortest path in fixed network within hybrid model for user $(user.uname). "
            enummaxtime = timelimit_pathenumeration 
            push!(bad_users, user)
            continue
        end
        lengthbound = shortest_fixed_path.cost

        # enumerate all paths
        @debug "Started enumerating paths in fixed network within hybrid model for user $(user.uname)"
        paths, runtime_enum = enumerate_user_paths(user, hndp, lengthbound, time_left)
        time_left = time_left - runtime_enum
        if time_left <= 0
            # we used up all enumeration time for this user
            @debug "Timeout while enumerating paths in fixed network within hybrid model for user $(user.uname)"
            enummaxtime = timelimit_pathenumeration 
            push!(bad_users, user)
            continue
        end
        upaths[user] = paths

        # update 'mintimeleft' to be the minimal amount of runtime left after path enumeration
        enumtime = runtime_sp + runtime_enum
        if  enumtime > enummaxtime
            enummaxtime = enumtime
        end
        @debug "Enumerating paths in fixed network within hybrid model for user $(user.uname) finished. Found $(length(paths)) paths."
    end

    # Generate second level primal, dual, and strong-duality constraints
    @debug "Enumeration part concluded with worst case runtime $enummaxtime. The following users were identified as bad users for which we add arc-based second-level constraints: $bad_users" 
    master_objs = Dict()
    sub_objs = Dict()
    for user in hndp.users
        if !(user in bad_users)
            # we fully enumerated all paths and build path model now
            oL1, oL2 = make_user_constraints_path(user, mip, hndp, upaths[user], A, x, true; l2binary=true)
            optL2 = @variable(mip, base_name = "optL2P$(user.uname)")
        else
            # we could not fully enuemrate all paths and use flow model now
            oL1, oL2, f = make_user_constraints(user, mip, hndp, A, x, true; bigM_function=derive_bigM_function(hndp, A, fixedBigM))
            optL2 = @variable(mip, base_name = "optL2A$(user.uname)")
        end

        # a small hack to see easily second level solution later in output
        @constraint(mip, optL2 == oL2, base_name = "objUser$(user.uname)")

        # set objective terms
        master_objs[user.uname] = oL1
        sub_objs[user.uname] = optL2
    end

    # add variable that represents the construction cost (easier to debug with it)
    constructioncost = sum(hndp.edge_price[a]*x[a] for a in A)
    @variable(mip, construction_cost_var)
    @constraint(mip, constructioncost == construction_cost_var, base_name="constructioncost")

    # set objective 
    @objective(mip, Min, constructioncost + sum(values(master_objs)))

   # build instance
    @info "Finished generation of compact hybrid MIP model that combines path and arc model ideas for passed HNDP instance" 
    master = MIPMaster(mip)
    rV = Instance(master, nothing)
    return rV, enummaxtime, length(bad_users)
end





########## Auxiliary functions for generating subproblem solver instances ##########
function build_Astar_user(user::User, hndp::HNDPwC, A)
    # Check that weight parameters are included. Otherwise, we currently do not support solving the subproblem with labeling algorithms if weigh parameters are not provided!
    if isnothing(hndp.minweights)
        throw(ArgumentError("Currently, solving an instance of the HNDP without weights is not supported by the labeling algorithms implemented. 
        Please use the option to solve the subproblems as MIP (which has okay performance for shortest paths conputations that are needed and can better deal with negative cycles.)"))
    end

    # precompute apsp according to risk and cost 
    spr = floyd_warshall_shortest_paths(hndp.mygraph, user.mrisk)
    spc = floyd_warshall_shortest_paths(hndp.mygraph, user.mcost)
    spw = floyd_warshall_shortest_paths(hndp.mygraph, user.mweight)

    # init structure 
    check_cycles = false  # no negative cycles can exist
    structure = HNDPwC_Structure(hndp.mygraph, user.origin, user.destination, user.weighlimit, check_cycles, user.mrisk, user.mcost, user.mweight, spc.dists, spr.dists, spw.dists, [])

    # init solver for user
    capa = Dict(a => 1 for a in A)
    return AStarSolver(user.uname, A, structure, capa, Inf, false)
end

"""
    build_MIP_user(user::User, hndp::HNDPwC, A, Adec, solver::SolverWrapper; cycle_free=false)

Generate a Subsolver that solves the passed users subproblem by formulating it as an MIP. 
# Arguments
    - 'user': The user whose subproblem we consider
    - 'hndp::HNDPwC': The HNDP instance.
    - 'A': The set of arcs in the network
    - 'Adec': The subset of arcs that are decision arcs. 
    Most algorithms that use subsolvers can take advantage if only a subset of arcs are part of the first-level decision space. 
    - 'solver::SolverWrapper': The MIP solver to use

    - 'cycle_free::Bool': If true, generate Benders-like cuts to ensure that the found solution is cycle free (only usefull for GBCsolver where negative cycles can incur)
"""
function build_MIP_user(user::User, hndp::HNDPwC, A, Adec, solver::SolverWrapper; cycle_free=false)
    # set up sub MIP (min-cost flow for passed user)
    sub = Model()
    set_optimizer(sub, () -> get_next_optimizer(solver))
    @variable(sub, xc[A], Bin)
    oL1, oL2, f = make_user_constraints(user, sub, hndp, A, xc, false; l2binary=true) 
    @objective(sub, Min, oL2)

    # additionally, fix non-decision arc variables in sub to 1
    for et in A
        if !(et in Adec)
            @constraint(sub, xc[et] == 1, base_name = "nondec$(et))")
        end
    end 

    # if requested, add cycle breaking constraints (via callback)
    if cycle_free
        # build up helper model that solves the shortest path for sepration
        subhelper = Model(Gurobi.Optimizer) # TODO: check if same enviroment prevents sub using callback
        #set_optimizer(subhelper, () -> get_next_optimizer(solver))
        @variable(subhelper, x_helper[A], Bin)
        oL1help, oL2help, y_helper = make_user_constraints(user, subhelper, hndp, A, x_helper, false; l2binary=true) # binary vars essential for cycle break
        @objective(subhelper, Min, oL2help)
        set_attribute(subhelper, MOI.NumberOfThreads(), 1) # multi-threads are not fun...

        # additionally, fix non-decision arc variables in sub to 1
        for et in A
            if !(et in Adec)
                @constraint(subhelper, x_helper[et] == 1, base_name = "nondecH$(et))")
            end
        end 

        @debug "helper mip: $subhelper"

        # disable output of helper model (this would be very very confusing)
        set_silent(subhelper) 

        #prepare big M's for Benders-like cuts
        bigMs = derive_bigM_function(hndp, A, true) # always use fixed bigM if available

        # Now use internal JuBiC function for generating Benders-like cuts callback for Subsolver MIP
        saved_cuts = Dict()  # a mapping of master solution to pairs (found lazy constraints, value of second level solution). 
        
        
        # build function that adds cuts to MIP if needed
        #= function extra_blc(timelimit)
            current_L2solution_sub = JuMP.value(oL2)
            x_vals = Dict(a => round(JuMP.value(xc[a])) for a in A)  # round to integer 0 or 1 values

            msolkey = JuBiC.key_master_sol(x_vals, A)
            if haskey(saved_cuts, msolkey)
                # if we already solved this sub_problem, we do not need to solve it again as we already added corresponding cut
                subopt = saved_cuts[msolkey][2]
                need_cut = (current_L2solution_sub > subopt)
                if need_cut
                    @debug "Need to recover cut with rhs $subopt while current master solution is $current_L2solution_sub"
                    cutopt = @constraint(sub, oL2 <= saved_cuts[msolkey][1]) 
                    @debug "Recovered to Subsolver of user $(user.uname) the additional constraint $cutopt"
                else
                    @debug "We recovered cut but no need to add it as cut rhs is $subopt while current master solution is $current_L2solution_sub"
                    return false, 0  # we ignore runtime for the check for internal runtime tracking
                end
            end

            @debug "Starting Benders-like cuts within Subsolver MIP for user $(user.uname) by solving helper MIP for sub solution $x_vals."
            # set time limit for solving
            set_time_limit_sec(subhelper, timelimit)
            # fix values for linking variables
            @constraint(subhelper, fixc[a=A], x_helper[a] == round(x_vals[a]))
            # now, solve the helper MIP
            optimize!(subhelper)

            # handle solution 
            solve_time = JuMP.solve_time(subhelper)
            status = JuMP.termination_status(subhelper)
            @debug "Helper MIP was solved with status $status"
            if status == MOI.OPTIMAL
                y_vals = JuMP.value.(y_helper)
                subopt = JuMP.objective_value(subhelper)
            elseif status == MOI.TIME_LIMIT
                return false, solve_time
            else
                error("Helper MIP of user $(user.uname) termonated with status $status")
            end

            # generate cut coefficients
            bigMterms = 0
            for a in A
                bigMterms +=
                    bigMs(a, user.uname) *
                    y_vals[a] *
                    (1 - xc[a])
            end

            # save found cut
            saved_cuts[msolkey] = (subopt + bigMterms, subopt)

            # clean up for next run by removing constraints enforcing model properties 
            @debug "clean up helper model for next callback call"
            for a in A
                delete(subhelper, fixc[a])
            end
            unregister(subhelper, :fixc)

            # check if cut is needed to be added
            @debug "Solved helper MIP and found solution with value $subopt while current master solution is $current_L2solution_sub"
            need_cut = (current_L2solution_sub > subopt)
            if need_cut
                cutopt = @constraint(sub, oL2 <= subopt + bigMterms) 
                @debug "Added to Subsolver of user $(user.uname) the additional constraint $cutopt"
                
                # finished adding cut. Check if new cut implies resolving and then return
                return true, solve_time # always resolve as we added a constraint. Time to build cut neglected in internal runtime computations
            else
                @debug "No need to add cut as current Subsolver solution is valid"
                return false, solve_time
            end
        end =#
    end 

    # build Subsolver instance
    link_constraints_capacities = Dict(a => 1.0 for a in Adec)
    extra_blc_sub = if cycle_free time -> extra_cuts_benderslike_JuMP(subhelper, sub, A, oL2, xc, x_helper, y_helper, a -> bigMs(a, user.uname), saved_cuts, time) else time -> (false, 0) end
    subS = SubSolverJuMP(
        user.uname,
        sub,
        Adec,
        xc[Adec],
        f[Adec],
        link_constraints_capacities,
        oL1,
        oL2,
        extra_blc_sub
    )
    set_silent(sub) 
    return subS
end



########## Auxiliary functions for generating first and second level constraints ##########
"""
    make_user_constraints(user::User, hpr, hndp::HNDPwC, A, xvars, dualconst; l2binary=true, bigM_function = nothing, indicatorconst = true, with_var_bounds=false)

Adds the variables and constraints for the high-point-relaxation to the model and returns (L1 obj, L2 obj).
The 'dualconst' is a bool. If true, dual constraints and the strong duality constraint are added additionally to primal constraints.
If you want to use dual constraints, please also provide a function 'bigM_function' that takes arc and user and return a valid big M for them.

# Arguments
    - 'user': The user for which we generate the cuts.
    - 'hpr': The JuMP model to which we add cuts.
    - 'hndp': The HNDP instance
    - 'A': The set of arcs (both decision and non-decision arcs)
    - 'xvars': The first-level linking variables.
    - 'dualconst': If true, add also dual constraints to the MIP

# Optional:
    - 'l2binary=true': If true, the second-level flow variables are set to be binary.
    - 'bigM_function = nothing': If you add dual constraints, pass a function here that gives the big M's. It should take (a, u) as arguments where 'a' is the resource 'u' is the user name.
    - 'indicatorconst = true': If true, use indicator constraints for dual constraints. The big M are still used as bounds on second-level variables.
    - 
"""
function make_user_constraints(user::User, hpr, hndp::HNDPwC, A, xvars, dualconst; l2binary=true, bigM_function = nothing, indicatorconst = true, with_var_bounds=true)
    # add flow variables
    if l2binary
        f = @variable(hpr, [i = A], Bin, base_name = "f$(user.uname)_")
    else
        f = @variable(hpr, [i = A], base_name = "f$(user.uname)_")
        for vf in f
            set_lower_bound(vf, 0)
        end 
    end

    # resource constraint (only if weigh parameter provided)
    if !isnothing(hndp.minweights) 
        @info "Added second-level weight capacity constraints to model or MIP submodel of user $(user.uname)"
        @constraint(
            hpr,
            sum(user.mweight[a...] * f[a] for a in A) <= user.weighlimit,
            base_name = "capacity$(user.uname)_"
        )
    end

    # flow conservation
    for n in vertices(hndp.mygraph)
        # set rhs
        demand::Int = 0
        if n == user.origin
            demand = 1
        elseif n == user.destination
            demand = -1
        end

        # build constraint for n
        flowcons =
            reduce(+, [f[(n, o)] for o in outneighbors(hndp.mygraph, n)], init=0.0) -
            reduce(+, [f[(i, n)] for i in inneighbors(hndp.mygraph, n)], init=0.0)
        @constraint(hpr, flowcons == demand, base_name = "flow$(user.uname)_$n")
    end

    # link first and second level
    for a in A
        @constraint(hpr, f[a] <= xvars[a], base_name = "link$(user.uname)_$a")
    end

    # build contribution to master objective
    objL1 = sum(user.mrisk[a...] * f[a] for a in A)
    objL2 = sum(user.mcost[a...] * f[a] for a in A)

    @info "Added second-level arc-flow constraints to model or MIP submodel of user $(user.uname). Overall, $(length(f)) variables were generated."

    # add dual constraints if requested
    if dualconst    
        @debug "adding dual constraints with settings indicatorconst=$indicatorconst and with_var_bounds=$with_var_bounds"
        if isnothing(bigM_function)
            throw(ArgumentError("There was no big M function provided to generate dual constraints of arc-based model."))
        end
        
        # generate variables
        vpi = @variable(hpr, [i = vertices(hndp.mygraph)], base_name = "pi$(user.uname)_")
        if with_var_bounds
            maxbigM = maximum(bigM_function(a, user.uname) for a in A)
            for v in vpi
                set_upper_bound(v, maxbigM)
                set_lower_bound(v, 0)
            end
        end

        # potential conservation constraints
        for a in A
            if indicatorconst
                @constraint(hpr, xvars[a] --> {vpi[a[1]]-vpi[a[2]] <= user.mcost[a...]}, base_name = "potential_ind$(user.uname)_$(a)")
            else
                @constraint(hpr, vpi[a[1]]-vpi[a[2]] <= user.mcost[a...] + bigM_function(a, user.uname)*(1-xvars[a]), base_name = "potential$(user.uname)_$(a)")
            end
        end

        # strong duality constraint
        @constraint(hpr, objL2 == vpi[user.origin] - vpi[user.destination], base_name = "strong_dual$(user.uname)")
        
        @info "Added dual second-level arc-flow and strong-duality constraints to model or MIP submodel of user $(user.uname). Overall, additional $(length(vpi)) variables were generated."

    end

    return objL1, objL2, f
end


"""
    derive_bigM_function(hndp::HNDPwC)

Compute a function that returns, for arc and user, a valid big M.
First, try to compute the shortest path in fixed in network that gives natural big M. 
If does not work, i.e., there is no path in fixed network, use a default big M based on n-1 most expensive arcs (where n is the number of nodes). 
# Parameter
    - 'hndp': The HNDP instance.
    - 'A': the set of all arcs of the graph.
    - 'fixedBigM': If true, first try to compute smart big M based on fixed network. If no path in fixed network found, use default big M.
"""
function derive_bigM_function(hndp::HNDPwC, A, fixedBigM)
    big_m_user = Dict()
    for user in hndp.users
        if fixedBigM
            # first, try to find paths in fixed network that forms natural upper bound
            # TODO: assume that finding such shortest path (or proving none exists) does not consume any runtime. Maybe handle better in future versions
            sp, runtime_sp = shortest_unforbiddable_path(user, hndp, Inf)
            if runtime_sp >= 0 # non-negative runtime indicates we found path
                @debug "For user $(user.uname), we found a shortest path in fixed network and use its length $(sp.cost) as big M."
                big_m_user[user.uname] = sp.cost
                continue
            end
        end

        # No good big M found. Fall back: build big M function: M = sum of n-1 largest elements (n=number nodes)
        arccost = [user.mcost[a...] for a in A]
        costly = arccost[partialsortperm(arccost, 1:(nv(hndp.mygraph)-1); lt = Base.isgreater)]
        big_m_user[user.uname] = sum(costly)
        @debug "For user $(user.uname), we found no shortest path in fixed network. Instead, the big M is approximated by $(big_m_user[user.uname])."
    end
    bigMfunc(a, uname) = begin
        big_m_user[uname]
    end
    return bigMfunc
end


"""
    make_user_constraints_path(user, hpr, hndp::HNDPwC, A, xvars; l2binary = true)

Adds the variables and constraints for the high-point-relaxation to the model and returns (L1 obj, L2 obj).
The 'dualconst' is a bool. If true, dual constraints and the strong duality constraint are added additionally to primal constraints.
"""
function make_user_constraints_path(user, mip, hndp::HNDPwC, P, A, xvars, dualconst; l2binary=true)
    # generate path variables
    if l2binary
        f = @variable(mip, [i = P], Bin, base_name = "f$(user.uname)_")
    else
        f = @variable(mip, [i = P], base_name = "f$(user.uname)_")
        for vf in f
            set_lower_bound(vf, 0)
        end
    end

    # Note: No resource constraints or similar needed as indirectly included in path set

    # select one path
    @constraint(mip, sum(f) == 1, base_name = "Cover$(user.uname)")

    # link first and second level
    for p in P
        Ap = p.arcs
        for a in Ap
            @constraint(mip, f[p] <= xvars[a], base_name = "link$(user.uname)_$(a)_$(string(p))")
        end
    end

    # build contribution to master objective
    objL1 = sum(p.risk * f[p] for p in P)
    objL2 = sum(p.cost * f[p] for p in P)

    @info "Added path-based primal constraints for user $(user.uname). Overall, $(length(f)) variables were generated."

    # add dual constraints if requested
    if dualconst
        # max cost path for big M
        max_subpath = maximum(p -> p.cost, P)
        @debug "Cost of max path for user $(user.uname) is $max_subpath"

        # generate variables
        vpi = @variable(mip, base_name = "pi$(user.uname)")

        # Benders-like constraints 
        for p in P
            @constraint(mip, vpi <= p.cost + sum([(max_subpath-p.cost)*(1-xvars[a]) for a in p.arcs]), base_name = "Benderslike$(user.uname)_$(string(p))")
        end

        # strong duality constraint
        @constraint(mip, objL2 == sum(vpi), base_name = "strong_dual$(user.uname)")

        @info "Added path-based dual and strong-duality constraints for user $(user.uname)"
    end

    return objL1, objL2
end



########## Auxiliary functions for enumerating user paths ##########
struct Subpath
    mynode # the node this paths ends in
    arcs # the list of arcs that this paths consists of
    cost # the computed cost of this path
    risk # the computed risk of the path
    weight # the weight of the paths till now
end

function Base.:(==)(sp1::Subpath, sp2::Subpath)
    # as we have no parallel arcs and unique start node, two paths are equal if they have the same set of arcs they traversed
    return Set(sp1.arcs) == Set(sp2.arcs)
end

function Base.hash(sp::Subpath, h::UInt)
    return hash(Set(sp.arcs), h)
end

function Base.string(paths::Subpath)
    if isempty(paths.arcs)
        return ""
    end

    nodes = [paths.arcs[1][1]]
    for (u, v) in paths.arcs
        push!(nodes, v)
    end
    
    return nodes
end

"""
    shortest_unforbiddable_path(user::User, hndp::HNDPwC, timelimit)

Compute the shortest path in the subgraph given by fixed network, i.e., non-decision arcs. 
If the is no path in fixed network, return an empty path and runtime of -1.
In case of timeout, return empty path and timelimit as runtime.
"""
function shortest_unforbiddable_path(user::User, hndp::HNDPwC, timelimit)
    # compute the shortest path in the fixed network. Gives an error if no path is found. In case of timeout, return empty path. Also alwqays returns runtime needed

    #auxiliary function
    function neighboring_shortest_unforbiddable_path(label::Subpath)
        neighbors_list = []
        for n in neighbors(hndp.mygraph, label.mynode)
            # check if fixed edge. If not, scip
            if ((label.mynode, n) in hndp.edgeA)
                continue
            end

            # check that we do not close a cycle (without it we run into issues with 0 cost cycles)
            has_cycle=false
            for a in label.arcs
                if n in a
                    # we have a cycle
                    has_cycle = true
                    break
                end
            end
            if has_cycle
                continue
            end

            # generate new label 
            narcs = copy(label.arcs) 
            push!(narcs, (label.mynode, n))
            ncost = label.cost + user.mcost[label.mynode, n]
            nrisk = label.risk + user.mrisk[label.mynode, n]
            nweight = label.weight + user.mweight[label.mynode, n]
            newlabel = Subpath(n, narcs, ncost, nrisk, nweight)
            push!(neighbors_list, newlabel)
        end
        return neighbors_list
    end
    

    # TODO: One could implement it smarter but this basis version should do it
    starts = Subpath(user.origin, [], 0, 0, 0)
    ends = Subpath(user.destination, [], 0, 0, 0)
    heuristicf_none(lv::Subpath, lw::Subpath) = 0
    isgoal(label::Subpath) = label.mynode == user.destination
    costf(lv::Subpath, lw::Subpath) = user.mcost[lv.mynode, lw.mynode]
    dominance(lv, lw) = false
    solution = JuBiC.shortest_path_labeling(
            neighboring_shortest_unforbiddable_path,
            starts,
            ends,
            heuristicf_none,
            costf,
            isgoal,
            dominance,
            Inf,
            timelimit,
        )
    if solution.status == JuBiC.LS_NO_SOLUTION
        return solution.goal, -1
        #throw(ArgumentError("User $(user.uname) does not posses a path in the fixed network!"))
    elseif solution.status == JuBiC.LS_TIMEOUT
        return solution.goal, timelimit
    elseif solution.status == JuBiC.LS_OPTIMAL
        return solution.goal, solution.runtime
    else
        throw(ArgumentError("Return status $(solution.status) not known."))
    end
end




"""
    enumerate_user_paths(user::User, hndp::HNDPwC, lengthbound, timelimit)

Enumerates all paths of the 'user' up to the given 'lengthbound'. Note that only cycle-free paths are enumerated (otherwise, we could go in circles for very long time).
Throw exception if no path is found, because we should at least find the path in fixed network
Returns the found paths as list of 'Subpath' objects and the time required for enumeration. In case of timeout, the paths found till now and the timelimit as runtime are returned. 
If no path is found, throws and exception.
"""
function enumerate_user_paths(user::User, hndp::HNDPwC, lengthbound, timelimit)
    # precompute shortest path according to cost for all pairs of nodes to be able to cut of paths that will not reach goal within the lengthbound
    mincostpairs = floyd_warshall_shortest_paths(hndp.mygraph, user.mcost)

    #auxiliary function
    function neighboring_enumerate_user_paths(label::Subpath)
        neighbors_list = []
        for n in neighbors(hndp.mygraph, label.mynode)
            # check weight bound
            nweight = label.weight + user.mweight[label.mynode, n]
            if !isnothing(user.weighlimit) && nweight > user.weighlimit
                continue 
            end

            # check that we do not close a cycle
            has_cycle=false
            for a in label.arcs
                if n in a
                    # we have a cycle
                    has_cycle = true
                    break
                end
            end
            if has_cycle
                continue
            end

            # generate new label 
            narcs = copy(label.arcs) 
            push!(narcs, (label.mynode, n))
            nrisk = label.risk + user.mrisk[label.mynode, n]
            ncost = label.cost + user.mcost[label.mynode, n] # cost for new subpath till neighbour n
            newlabel = Subpath(n, narcs, ncost, nrisk, nweight)

            # check length bound can still be obtained
            d = user.destination
            if n != d && ncost + mincostpairs.dists[n, d] > lengthbound
                #@debug "Label $newlabel cut of because remaining distance is $(ncost + mincostpairs.dists[n, d]) but lengthbound=$lengthbound"
                continue
            end

            # all checks passed. Add label to feasible labels list
            push!(neighbors_list, newlabel)
        end
        return neighbors_list
    end

    @debug "Starting enumeration for user $(user.uname) for lengthbound=$lengthbound"

    # TODO: is definitly not the most efficient implementation but hopefully it is correct
    starts = Subpath(user.origin, [], 0, 0, 0)
    ends = Subpath(user.destination, [], 0, 0, 0)
    isgoal(label::Subpath) = label.mynode == user.destination
    costf(lv::Subpath, lw::Subpath) = user.mcost[lv.mynode, lw.mynode]
    dominance(lv, lw) = false
    solution = JuBiC.enumerate_all_paths(
            neighboring_enumerate_user_paths,
            starts,
            ends,
            costf,
            isgoal,
            dominance,
            lengthbound,
            timelimit,
        )
    if solution.status == JuBiC.LS_NO_SOLUTION
        throw(ArgumentError("User $(user.uname) does not posses a path for length bound $(lengthbound) despite it been reported as length of shortest path in fixed network!"))
    elseif solution.status == JuBiC.LS_TIMEOUT
        @debug "Enumeration reached the time limit of $timelimit for user $(user.uname) "
        return solution.goal, timelimit
    elseif solution.status == JuBiC.LS_OPTIMAL
        @debug "Enumeration concluded for user $(user.uname). We found $(length(solution.goal)) paths. "
        return solution.goal, solution.runtime
    else
        throw(ArgumentError("Return status $(solution.status) not known."))
    end

end

