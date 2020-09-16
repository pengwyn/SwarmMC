import Argon

cd(joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks")) do
    params = Argon.SetupParams()
    thorough_run || SpeedUp!(params)
    props = BunchedPropagate(params,2)
end
