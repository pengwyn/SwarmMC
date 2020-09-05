
module NessRobsonSharing

using Reexport ; @reexport using SwarmMC

NRCrossSection(eps) = 1u"Å^2"

function CollFreqs(params, gas, ptype, iss, isratio)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(10u"Å^2"))
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_MANUAL(NRCrossSection), threshold=10eV)
    ionisation = CreateCollFreq(params, gas, ptype, "ionisation", CFS_IONISATION(), GENCF_MANUAL(NRCrossSection), threshold=15eV, ion_sharing_style=iss, ion_sharing_ratio=isratio, new_ptype=ptype)

    [elastic, inelastic, ionisation]
end


@xport function SetupParams(ETd=300u"Td", S=-1.)
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1eV)))

    m0 = 25amu
    ρ = 1e-6u"Å^-3"
    p.E_field = EFromETd(ETd, ρ) * [0,0,1]

	p.eps_grid = LogRange(-5, 3, 1001)*eV

	gas = GAS(; m0, ρ)
    electron = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    
    if S == -1
        iss = ISS_AFE()
        isratio = 0.
        Sstr = "AFE"
    else
        iss = ISS_FRACTION()
        isratio = S
        Sstr = string(S)
    end

    p.gen_cf = (args...) -> CollFreqs(args..., iss, isratio)

    p.t_grid = LogRange(5, 8, 101) * SwarmMC.uT

    # p.max_substepsize = ???

    p.save_name = MakeSaveName("NRS", S=Sstr, ETd=ETd)
    
    p
end

end
