

include("ReidRamp.jl")
using .ReidRamp

@CheckTurns for E in [1., 12., 24.]*u"Td"
    p = SetupParams(E)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
