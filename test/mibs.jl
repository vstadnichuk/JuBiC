using JuBiC, BilevelJuMP, MibS_jll


function test_mibs()
    # Create a simple bilevel model
    model = BilevelModel()
    @variable(Upper(model), x[1:2], Bin)
    @variable(Lower(model), y[1:2], Bin)

    @constraints(Lower(model), begin
        l1, y[1] <= x[1]
        l2, y[2] <= x[2]
        l3, y[1] >= 1  # somehow using y[1] == 1 does not work
        l4, y[1] <= 1
    end)

    @objective(Upper(model), Min, x[1] - x[2] + 10 * y[2])
    @objective(Lower(model), Min, -y[1] - y[2])

    # Construct the Master and Parameter
    master = MibSMaster(model)
    inst = Instance(master, nothing)
    """folder_path = logging_folder * "/MibS_Solver"
    mkpath(folder_path)"""
    params = MibSparam(false, logging_folder, JuBiC.RunStats())

    # Solve the instance
    stats = solve_instance!(inst, params)

    @test stats.data["Opt"] ≈ 1
end


test_mibs()
