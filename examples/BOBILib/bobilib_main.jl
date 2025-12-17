using JuBiC, JuMP, Gurobi, BilevelJuMP

include("../logging.jl")

optimizer = Gurobi.Optimizer
const instances_dir = "./examples/BOBILib/instances/"
const instances_to_solve = []
const test_gbc = true
const test_gbc_without_preprocessing = false
const test_gbc_without_partial_decomposition = false
const test_mibs = true


function _solve_GBC(mps_data, aux_data, logging_gbc, partial_decomposition, preprocessing)
    gbc_parameters = GBCparam(GurobiSolver(Gurobi.Env()), false, logging_gbc, "lp", PARETO_OPTIMALITY_ONLY)
    gbc_instance = get_GBC_instance(mps_data, aux_data, optimizer; partial_decomposition=partial_decomposition, preprocessing=preprocessing, stats=gbc_parameters.stats)
    gbc_statistics = solve_instance!(gbc_instance, gbc_parameters)
    return gbc_statistics
end

function _solve_MibS(mps_data, aux_data, logging_mibs)
    bilevel_parameters = MibSparam(false, logging_mibs)
    bilevel_instance = get_MibS_instance(mps_data, aux_data; stats=bilevel_parameters.stats)
    bilevel_statistics = solve_instance!(bilevel_instance, bilevel_parameters)
    return bilevel_statistics
end


function main()
    logging_file = init_logging_folder()
    logging_gbc = joinpath(logging_file, "GBC_Solver")
    logging_gbc_without_preprocessing = joinpath(logging_file, "GBC_Solver(Without_Preprocessing)")
    logging_gbc_without_partial_decomposition = joinpath(logging_file, "GBC_Solver(Without_Partial_Decomposition)")
    logging_mibs = joinpath(logging_file, "MibS_Solver")

    if test_gbc
        mkpath(logging_gbc)
    end

    if test_gbc_without_preprocessing
        mkpath(logging_gbc_without_preprocessing)
    end

    if test_gbc_without_partial_decomposition
        mkpath(logging_gbc_without_partial_decomposition)
    end

    if test_mibs
        mkpath(logging_mibs)
    end

    statistics_list = []

    for instance in instances_to_solve
        mps_path = joinpath(instances_dir, instance * ".mps")
        aux_path = joinpath(instances_dir, instance * ".aux")

        if !isfile(mps_path) || !isfile(aux_path)
            @warn "Instance files for $instance not found. Skipping."
            continue
        end

        if test_gbc
            println("Solving instance $instance with GBC")
            statistic = _solve_GBC(mps_path, aux_path, logging_gbc, true, true)
            push!(statistics_list, statistic)
        end

        if test_gbc_without_preprocessing
            println("Solving instance $instance with GBC (without preprocessing)")
            statistic = _solve_GBC(mps_path, aux_path, logging_gbc_without_preprocessing, true, false)
            push!(statistics_list, statistic)
        end

        if test_gbc_without_partial_decomposition
            println("Solving instance $instance with GBC (without partial decomposition)")
            statistic = _solve_GBC(mps_path, aux_path, logging_gbc_without_partial_decomposition, false, true)
            push!(statistics_list, statistic)
        end

        if test_mibs
            println("Solving instance $instance with MibS")
            statistic = _solve_MibS(mps_path, aux_path, logging_mibs)
            push!(statistics_list, statistic)
        end
    end

    output_csv = joinpath(logging_file, "statistics.csv")
    print_stats_to_csv(statistics_list, output_csv)
    println("Done with all tests. Statistics written to $output_csv")
end

main()
