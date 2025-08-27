using JuBiC, Test, Logging, Dates
using Gurobi

optimizer = Gurobi.Optimizer


# create logging folder
include("../examples/logging.jl")
logging_folder = init_logging_folder()

logger = JuBiC.new_file_logger(logging_folder * "/debuglog.txt", true)
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
end
