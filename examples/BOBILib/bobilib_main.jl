using JuBiC, JuMP, Gurobi, BilevelJuMP

include("../logging.jl")

optimizer = Gurobi.Optimizer
const instances_dir = "C:\\Users\\stadnichuk\\OneDrive - rwth-aachen.de\\Documents\\Arbeit\\Julia Bilevel\\L_Shape_tests\\JuBiC\\examples\\BOBILib\\instances"
const instances_to_solve = ["K5010W07.KNP"]
const test_gbc = true
const test_mibs = true
const test_mibs_transform = false


function _solve_GBC(mps_data, aux_data, logging_gbc)
    gbc_instance = get_GBC_instance(mps_data, aux_data, optimizer)
    gbc_parameters = GBCparam(GurobiSolver(Gurobi.Env()), false, logging_gbc, "lp", PARETO_OPTIMALITY_ONLY)
    gbc_statistics = solve_instance!(gbc_instance, gbc_parameters)
    return gbc_statistics
end

function _solve_MibS(mps_data, aux_data, logging_mibs)
    bilevel_instance = get_MibS_instance(mps_data, aux_data)
    bilevel_parameters = MibSparam(false, logging_mibs)
    bilevel_statistics = solve_instance!(bilevel_instance, bilevel_parameters)
    return bilevel_statistics
end

function _solve_MibS_transformation(mps_data, aux_data, logging_mibs)
    gbc_instance = get_GBC_instance(mps_data, aux_data, optimizer)
    bilevel_instance = transform_GBC_to_MibS(gbc_instance)
    bilevel_parameters = MibSparam(false, logging_mibs)
    bilevel_statistics = solve_instance!(bilevel_instance, bilevel_parameters)
    return bilevel_statistics
end

function main()
    logging_file = init_logging_folder()
    logging_gbc = joinpath(logging_file, "GBC_Solver")
    logging_mibs = joinpath(logging_file, "MibS_Solver")
    logging_mibs_transform = joinpath(logging_file, "MibS_Solver(Transformation)")

    if test_gbc
        mkpath(logging_gbc)
    end

    if test_mibs
        mkpath(logging_mibs)
    end

    if test_mibs_transform
        mkpath(logging_mibs_transform)
    end

    statistics_list = []

    for instance in instances_to_solve
        mps_path = joinpath(instances_dir, instance * ".mps")
        aux_path = joinpath(instances_dir, instance * ".aux")

        if !isfile(mps_path)
            @warn "Instance mps files for $instance not found under path $mps_path. Skipping."
            continue
        end
        if !isfile(aux_path)
            @warn "Instance aux files for $instance not found under path $aux_path. Skipping."
            continue
        end
        

        if test_gbc
            println("Solving instance $instance with GBC")
            statistic = _solve_GBC(mps_path, aux_path, logging_gbc)
            new_stat!(statistic, "instance", instance)
            push!(statistics_list, statistic)
        end

        if test_mibs
            println("Solving instance $instance with MibS")
            statistic = _solve_MibS(mps_path, aux_path, logging_mibs)
            new_stat!(statistic, "instance", instance)
            push!(statistics_list, statistic)
        end

        if test_mibs_transform
            println("Solving instance $instance with MibS (Transformation)")
            statistic = _solve_MibS_transformation(mps_path, aux_path, logging_mibs_transform)
            new_stat!(statistic, "instance", instance)
            statistic.data["Solver"] = "MibSSolver(Transformation)"
            push!(statistics_list, statistic)
        end
    end

    output_csv = joinpath(logging_file, "statistics.csv")
    print_stats_to_csv(statistics_list, output_csv)
    println("Done with all tests. Statistics written to $output_csv")
end

main()
