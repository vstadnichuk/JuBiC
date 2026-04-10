using JuBiC, Test

"""
    test_batch_solve_simple_bilevel()

Test the JSON-based batch runner on one simple bilevel example with two JuBiC
solvers. The batch configuration runs one GBC, one BlC, and one BlCLag
experiment, each with its own named instance factory. The test checks that all
experiments solve successfully, that the reported objective values are correct,
and that the batch runner writes the aggregated CSV and resolved per-experiment
config files.
"""
function test_batch_solve_simple_bilevel()
    batch_output_folder = "test/logs/batch_runs"
    if isdir(batch_output_folder)
        rm(batch_output_folder; recursive=true, force=true)
    end

    results = solve_batch!(
        [
            "gbc_instance" => generate_gbc_simple_bilevel_instance,
            "blc_instance" => generate_blc_simple_bilevel_instance,
            "blclag_instance" => generate_blclag_simple_bilevel_instance,
        ],
        "test/batch_config.json",
    )

    @test haskey(results, "gbc_instance")
    @test haskey(results["gbc_instance"], "gbc_simple")
    @test results["gbc_instance"]["gbc_simple"].data["Opt"] ≈ 1

    @test haskey(results, "blc_instance")
    @test haskey(results["blc_instance"], "blc_simple")
    @test results["blc_instance"]["blc_simple"].data["Opt"] ≈ 1

    @test haskey(results, "blclag_instance")
    @test haskey(results["blclag_instance"], "blclag_simple")
    @test results["blclag_instance"]["blclag_simple"].data["Opt"] ≈ 1

    @test isfile(joinpath(batch_output_folder, "batch_summary.csv"))
    @test isfile(joinpath(batch_output_folder, "gbc_simple", "gbc_instance", "resolved_experiment.json"))
    @test isfile(joinpath(batch_output_folder, "blc_simple", "blc_instance", "resolved_experiment.json"))
    @test isfile(joinpath(batch_output_folder, "blclag_simple", "blclag_instance", "resolved_experiment.json"))
end

test_batch_solve_simple_bilevel()
