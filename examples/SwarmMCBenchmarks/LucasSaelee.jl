module LucasSaelee

using SwarmMC

LSFunc(eps) = 0.1u"Å^2"*(eps/eV - 15.6)

function CollFreqs(params, gas, ptype, F)
    gencf = GENCF_MAXWELL(4.0u"Å^2" * sqrt(eV))
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), gencf)
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_MANUAL(LSFunc, 1-F), threshold=15.6eV)
    ionisation = CreateCollFreq(params, gas, ptype, "ionisation", CFS_IONISATION(), GENCF_MANUAL(LSFunc, F), threshold=15.6eV, new_ptype=ptype, ion_sharing_style=ISS_AFE())

    [elastic,inelastic,ionisation]
end

@xport function SetupParams( ; F=0.5, ETd=10u"Td", BHx=nothing, Bθ=0, ion_weight=false)
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))

    m0 = 1000u"mₑ"
    n0 = 1e-6u"Å^-3"

	gas = GAS(name=:gas,m0=m0,ρ=n0)
	ptype = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, ptype)

    p.E_field = EFromETd(ETd, n0) * [0,0,1]
    if BHx != nothing
        p.B_field = BFromBHx(BHx, n0) * [sind(Bθ),0,cosd(Bθ)]
    end
    p.gen_cf = (args...) -> CollFreqs(args..., F)

    p.t_grid = LogRange(1e5, 1e9, 101, prefactor=p.time_unit, as_val=true)
    p.max_substepsize = 1e7*p.time_unit

    p.ionisation_as_weight = TypeBool(ion_weight)

    p.meas_bins = (MEAS_BIN(:cum, :t), MEAS_BIN(:t))

    p.save_name = MakeSaveName("LucasSaelee" ; F, ETd, BHx, Bθ, ion_weight)

    p
end

end
