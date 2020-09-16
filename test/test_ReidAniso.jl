import ReidAniso
using NonNegLeastSquares

@testset "Model $model, method $method" for model in [nothing ; collect("ABCD")],
    method in [:func,:grid,:legendre]

    params = ReidAniso.SetupParams(model, method)
    thorough_run || SpeedUp!(params)
    props = BunchedPropagate(params,2)
end
