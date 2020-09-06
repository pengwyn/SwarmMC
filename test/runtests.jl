using Test
using SwarmMC

pushfirst!(LOAD_PATH, joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks"))

@testset begin
    if "MOD_LIST" in keys(ENV)
        mod_list = split(ENV["MOD_LIST"], ',') .|>  Symbol
    else
        mod_list =[:Hardsphere, :ConstIon, :LucasSaelee, :MagneticField, :Maxwell, :NessRobsonLoss, :NessRobsonSharing, :PercusYevick, :ReidAniso, :ReidRamp, :SuperModel, :Argon, :DaleTanhTest]
    end

    @testset "Benchmark $mod" for mod in mod_list
        filename = "test_$(string(mod)).jl"
        if isfile(filename)
            include(filename)
        else
            @eval import $mod
            params = @eval $mod.SetupParams()
            # Assuming that 2 particles will always be enough to test all features
            props = BunchedPropagate(params, 2)
        end
    end
end
