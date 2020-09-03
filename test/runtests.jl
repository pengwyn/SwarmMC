using Test
using SwarmMC

pushfirst!(LOAD_PATH, joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks"))

@testset begin
    @testset "Benchmark $mod" for mod in [:Hardsphere, :ConstIon, :LucasSaelee, :MagneticField, :Maxwell, :NessRobsonLoss, :NessRobsonSharing, :PercusYevick, :ReidAniso, :ReidRamp, :SuperModel, :Argon]
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
