
union!(LOAD_PATH, ["."])
using Maxwell

@CheckTurns for tmtr in [0., 293.]u"K",
                ETd in [nothing ; [1., 0.1, 10.]u"Td"]

    tmtr == 0. && ETd == nothing && continue
    
    p = SetupParams(tmtr, ETd)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
