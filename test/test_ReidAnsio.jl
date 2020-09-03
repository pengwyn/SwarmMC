import ReidAniso

@testset "Model $model, method $method" for model in [nothing ; collect("ABCD")],
    method in [:func,:grid,:legendre]

    params = Finalise(ReidAniso.SetupParams(model, method))
    props = BunchedPropagate(params,2)
end
