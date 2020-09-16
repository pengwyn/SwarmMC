using Test
using SwarmMC

pushfirst!(LOAD_PATH, joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks"))

thorough_run = parse(Bool, get(ENV, "SWARMMC_THOROUGH", "false"))

function SpeedUp!(params)
    @assert !thorough_run
    params.t_grid = LinRange(0*SwarmMC.uT, min(params.t_grid[end]/100, 1e3*SwarmMC.uT), length(params.t_grid))
end

@testset "All benchmarks" begin
    if "MOD_LIST" in keys(ENV)
        mod_list = split(ENV["MOD_LIST"], ',') .|>  Symbol
    else
        mod_list =[:ConstIon, :LucasSaelee, :MagneticField, :Maxwell, :NessRobsonLoss, :NessRobsonSharing, :ReidRamp, :SuperModel, :Argon, :DaleTanhTest, :Hardsphere]
    end

    @testset "Benchmark $mod" for mod in mod_list
        filename = "test_$(string(mod)).jl"
        if isfile(filename)
            include(filename)
        else
            @eval import $mod
            params = @eval $mod.SetupParams()
            thorough_run || SpeedUp!(params)
            # Assuming that 2 particles will always be enough to test all features
            props = BunchedPropagate(params, 2)
        end
    end
end
