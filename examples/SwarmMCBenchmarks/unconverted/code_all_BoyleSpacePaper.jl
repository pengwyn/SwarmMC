using BoyleSpacePaper

using Generic, SwarmMC

@CheckTurns for phi in reverse([0.0, 0.2, 0.3, 0.4])
    #p = SetupParams(phi, :drifted_maxwellian, :gaussian)
    p = SetupParams(phi, :delta, :delta)

    props = LoopMaxTime(p, 1)

    Save(p, props)
end
