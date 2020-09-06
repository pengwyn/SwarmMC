
using ForwardDiff, NonNegLeastSquares
include("ReidAniso.jl")
using .ReidAniso

@CheckTurns for model in ['A', 'B', 'C', 'D', nothing], method in [:legendre, :func, :grid]
    model == nothing && (method == :func || continue)

    p = SetupParams(model, method)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
