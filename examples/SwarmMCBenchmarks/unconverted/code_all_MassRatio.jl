
using SwarmMC
using Generic

using MassRatio

#while true
#for iter = 1:20

@CheckTurns while true
    for ratio in [1e-4, 0.1, 1.0, 2.0, 10.0]
        for Temp in [0., 293.]
            if ratio == 1e-4
                ETd_list = [0.1, 12.]
            else
                ETd_list = [2.0, 12.]
            end

            for ETd in ETd_list
                p = SetupParams(Temp, ETd, ratio)

                props = LoopMaxTime(p, 1)

                Save(p, props)
            end
        end
    end
end
