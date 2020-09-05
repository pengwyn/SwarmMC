
union!(LOAD_PATH, ["."])
using SuperModel

    #Eprelist = 1.0:0.5:10.0
    #Eprelist = 1.0:2.0:10.0
    #Eprelist = Eprelist[1:end-1]
    #Epowlist = -5:2
    #Epowlist = -5:-2
    #Epowlist = -5:-5
    #Epowlist = -3:2

    #Tlist = [0., 77., 293.]
    #Tlist = [0]

# For double checking errors
# Eprelist = 1.0
# Epowlist = -3:0
Eprelist = Float64[1, 3, 5, 7, 9]
Epowlist = -3:2
#Tlist = 293.
Tlist = [0., 77., 293.]*u"K"
recoil = [true,false]
# recoil = false

@CheckTurns for Epre = Eprelist, Epow = Epowlist
    E = Epre * 10.0^Epow

    for elastic_T = Tlist, inelastic_T = elastic_T, inelastic_recoil = recoil
        # (elastic_T == 77u"K" && 0.03 <= E <= 0.1) || (elastic_T == 0u"K" && E ≈ 1e-3) || continue

        p = SetupParams(E*u"Td", elastic_T, inelastic_T, inelastic_recoil)

        if E < 1e-4 && inelastic_T == 0u"K"
            if elastic_T == 0u"K"
                final_time = 5e13
            else
                final_time = 5e12
            end
        elseif E < 10
            final_time = 5e11
        else
            final_time = 5e10
        end

        final_time *= SwarmMC.uT

        p.t_grid = LogRange(final_time*1e-6, final_time, 101, as_val=true)

        props = LoopMaxTime(p, 1)

        Save(p, props)
    end
end
