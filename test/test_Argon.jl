import Argon

cd(joinpath(dirname(pathof(SwarmMC)), "..", "examples","SwarmMCBenchmarks")) do
    params = Argon.SetupParams()
    params.t_grid = LinRange(0*SwarmMC.uT, params.t_grid[end]/100, 10)
    props = BunchedPropagate(params,2)
end
