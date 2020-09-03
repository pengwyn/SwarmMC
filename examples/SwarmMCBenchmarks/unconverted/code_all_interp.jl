
using SwarmMC

@everywhere include("InterpTest.jl")

#while true
for iter = 1:20
    for E in [0.1, 0.5, 1., 2., 10.]
        for SetupFunc in [SetupParamsInterp, SetupParamsFunc]
            p = SetupFunc(E)
            props = LoopBunchedPropagate(p, 1*nworkers(), 1000, 1e8, 1e11)

            Save(p, props)
        end
    end
end
