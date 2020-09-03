

using ReidRamp

using SwarmMC
using DanUtils

@CheckTurns for E in [1., 12., 24.]*u"Td"
    p = SetupParams(E)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
