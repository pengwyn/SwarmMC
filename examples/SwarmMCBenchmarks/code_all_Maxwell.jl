
include("Maxwell.jl")
using .Maxwell

@CheckTurns for tmtr in [nothing, 293u"K"],
                ETd in [nothing ; [1., 0.1, 10.]u"Td"]

    tmtr == 0. && ETd == nothing && continue
    
    p = SetupParams(tmtr, ETd)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
