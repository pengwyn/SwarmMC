
include("Hardsphere.jl")
using .Hardsphere

@CheckTurns for Temp in [nothing, 293u"K"],
                ETd in [1., 0.1, 10.]*u"Td"
    p = SetupParams(Temp, ETd=ETd)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
