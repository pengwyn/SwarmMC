module ReidRamp

using SwarmMC

ReidInelastic(eps) = 10u"Å^2"*(eps/eV - 0.2)

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(6.0u"Å^2"))
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_MANUAL(ReidInelastic), threshold=0.2eV)

    [elastic, inelastic]
end

@xport function SetupParams(ETd=12.0Td, BHx=nothing)
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1eV)))

    m0 = 4*amu
    ρ = 1e-6u"Å^-3"
    p.E_field = EFromETd(ETd, ρ) * [0,0,1]

	gas = GAS(;m0, ρ)
    electron = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    
    p.gen_cf = CollFreqs

    p.t_grid = LinRange(0,1e9,101) * p.time_unit

    if BHx === nothing
        B = 0u"T"
        p.save_name = MakeSaveName("ReidRamp",E=ETd)
    else
        B = BFromBHx(BHx, n0, p)
        p.save_name = MakeSaveName("ReidRampB", E=ETd, B=BHx)
    end

    p
end

end
