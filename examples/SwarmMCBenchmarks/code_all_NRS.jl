
using SwarmMC
using DanUtils

using NessRobsonSharing

@CheckTurns for E in [300., 500., 800.]*u"Td"
        for S in [-1., 0., 1/4., 1/3., 1/2.]
            p = SetupParams(E, S)
            props = LoopMaxTime(p, 1000)

            Save(p, props)
        end
    end
end
