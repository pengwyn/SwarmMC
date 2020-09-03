
using SuperModelSingle
using Generic, SwarmMC

Eprelist = 1.0:0.5:10.0
Eprelist = Eprelist[1:end-1]
#Epowlist = -5:2
#Epowlist = -5:-2
Epowlist = -1:2

Tlist = [0., 77., 293.]

if false
    @CheckTurns for Epre = Eprelist, Epow = Epowlist
        E = Epre * 10.0^Epow
        E = round(E, sigdigits=5)

        for Temp = Tlist

            params = SetupParams(E, Temp)

            if E < 1e-4 && Temp == 0
                final_time = 5e12
            elseif E < 10.
                final_time = 5e11
            else
                final_time = 5e10
            end

            # if E < 1e-1
            #     npart = 10
            # else
            #     npart = 1
            # end

            #final_time = 1e9

            params.t_grid = GRID(final_time)

            props = LoopMaxTime(params, 1)

            Save(params, props)
        end
    end
end

# Some high-accuracy runs
@CheckTurns for E in [0.01, 0.1, 0.5, 1.0, 5.0, 10.0]
    for Temp = Tlist

        params = SetupParams(E, Temp)

        if E < 1e-4 && Temp == 0
            final_time = 5e12
        elseif E < 10.
            final_time = 5e11
        else
            final_time = 5e10
        end

        params.t_grid = GRID(final_time)

        props = LoopMaxTime(params, 10)

        Save(params, props)
    end
end
