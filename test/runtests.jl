using Test
using SwarmMC

pushfirst!(LOAD_PATH, joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks"))

@testset "All benchmarks" begin
    if "MOD_LIST" in keys(ENV)
        mod_list = split(ENV["MOD_LIST"], ',') .|>  Symbol
    else
        mod_list =[:ConstIon, :LucasSaelee, :MagneticField, :Maxwell, :NessRobsonLoss, :NessRobsonSharing, :PercusYevick, :ReidRamp, :SuperModel, :Argon, :DaleTanhTest, :Hardsphere]
    end

    @testset "Benchmark $mod" for mod in mod_list
        filename = "test_$(string(mod)).jl"
        if isfile(filename)
            include(filename)
        else
            @eval import $mod
            params = @eval $mod.SetupParams()
            # Forcing small time so that tests complete quick enough
            params.t_grid = LinRange(zero(eltype(params.t_grid)), params.t_grid[end]/100, 101)
            # Assuming that 2 particles will always be enough to test all features
            props = BunchedPropagate(params, 2)
        end
    end
end
