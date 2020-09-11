import ReidAniso
using NonNegLeastSquares

@testset "Model $model, method $method" for model in [nothing ; collect("ABCD")],
    method in [:func,:grid,:legendre]

    params = Finalise(ReidAniso.SetupParams(model, method))
    params.t_grid = LinRange(0*SwarmMC.uT, params.t_grid[end]/100, 101)
    props = BunchedPropagate(params,2)
end
