
using SwarmMC
using Generic

@everywhere using SuperModel

@CheckTurns begin
    Tlist = [0., 77., 293.]

    E = 1e-7
    for elastic_T = Tlist, inelastic_T = elastic_T, inelastic_recoil = true
        p = SetupParams(E, elastic_T, inelastic_T, inelastic_recoil)

        if inelastic_T == 0
            if elastic_T == 0
                if E < 1e-6
                    final_time = 1e15
                else
                    final_time = 5e13
                end
            else
                final_time = 5e12
            end
        else
            final_time = 5e11
        end

        p.t_grid = GRID(GRID_LOG(), final_time*1e-6, final_time)

        props = LoopMaxTime(p, 1)

        Save(p, props)
    end
end
