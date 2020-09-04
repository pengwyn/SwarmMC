import Argon

cd(joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks")) do
    params = Finalise(Argon.SetupParams())
    props = BunchedPropagate(params,2)
end
