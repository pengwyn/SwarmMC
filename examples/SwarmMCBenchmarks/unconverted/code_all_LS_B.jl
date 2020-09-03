
using Generic 
using LucasSaelee

@CheckTurns for F in [0., 0.5, 1.0], BHx in [100., 200., 1000.]
    p = SetupParams(F, BHx)

    props = LoopMaxTime(p, 1000)

    Save(p, props)
end
