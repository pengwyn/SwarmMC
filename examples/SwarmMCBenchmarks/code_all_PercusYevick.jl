
using SwarmMC
using DanUtils

using PercusYevick

@CheckTurns for E in [0.03, 0.3, 3., 10., 30.]*u"Td"
    #for E in [a*10^c for a in [1,3,6] for c in -3.:2.]
    #for E in [round(a*10^c,Int(-c+1)) for a in 1.:9. for c in -3.:2.]
    #for E in [round(a*10^c,Int(-c+1)) for a in 1.:9. for c in -1.:0.]
        for phi in [0., 0.2, 0.3, 0.4]
            p = SetupParams(E, phi)

            #props = LoopBunchedPropagate(p, 1*nworkers(), 10, final_time/100, final_time)
            props = LoopMaxTime(p, 1)

            #p2 = SwarmMC.FinaliseParams(p)
            #Save(p2, props)
            Save(p, props)
        end
    # end
end
