
using SwarmMC
using Generic

using ReidRamp

@CheckTurns for BHx in [50.0, 200.0]
        ETd = 12.0
        p = SetupParams(ETd, BHx)

        props = LoopMaxTime(p, 1)

        Save(p, props)
    end
end
