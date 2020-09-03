
using Generic
using SwarmMC
using Argon



@CheckTurns for just_el in [true, :allexc, false]
        for E in [1., 5., 10.]

            p = SetupParams(E, "Biagi", just_el)
            if just_el === true
                p.save_name = "ArForAriJustEl_E$E"
            elseif just_el === :allexc
                p.save_name = "ArForAriWithExc_E$E"
            elseif just_el === false
                p.save_name = "ArForAriAll_E$E"
            else
                error("Unknown just_el $just_el")
            end

            p.init_style = PIS_SEPARABLE(PIVS_ISOTROPIC(1.))
            p.ptype_list[1].cf_list[1].gas = GAS("argon", Argon.Ar_mass, 1e-6, 0.)

            p.measure_time_dists = :yes

            props = LoopMaxTime(p, 1)

            Save(p, props)
        end
    end
end
