using JuBiC, Test, Logging, Dates
using Gurobi

# Build each test model with a silent Gurobi environment so the suite reports
# test results instead of solver console output.
optimizer = () -> JuBiC.get_next_optimizer(JuBiC.GurobiSolver())


# create logging folder
logging_folder = JuBiC.init_logging_folder()

logger, io = JuBiC.new_file_logger(logging_folder * "/debuglog.txt", true)
with_logger(logger) do
    @testset "JuBiC Tests" begin
        @testset "Labeling Tests" begin
            include("labeling.jl")
        end
    end

    @testset "MibS Solver Tests" begin
        include("mibs.jl")
    end

    @testset "GBC Solver Tests" begin
        include("gbc.jl")
    end

    @testset "Bilevel SubSolver Tests" begin
        include("bilevel_subsolver.jl")
    end

    @testset "Batch Solver Tests" begin
        include("batch.jl")
    end
end
