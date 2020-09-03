
using SwarmMC

@everywhere include("SSF_Maxwell.jl")

#while true
for iter = 1:20
    for S in [0.1, 0.5, 1.0, 1.5, 2.0, 5.0]
        p = SetupParams(S)

        final_time = 1e10

        props = LoopBunchedPropagate(p, 1*nworkers(), 1000, final_time/100, final_time)

        Save(p, props)
    end
end
