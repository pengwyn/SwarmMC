module ConstIon

using Reexport ; @reexport using SwarmMC

function CollFreqs(params, gas, ptype, ion_style, ion_frac)
    elastic =    CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MAXWELL(4.0u"Å^2" * sqrt(eV)))
    ionisation = CreateCollFreq(params, gas, ptype, "ionisation", CFS_IONISATION(), GENCF_MAXWELL(0.01u"Å^2" * sqrt(eV)), new_ptype=ptype, ion_sharing_style=ion_style, ion_sharing_ratio=ion_frac)

    [elastic,ionisation]
end

@xport function SetupParams(; ionweight=false, gns=GNS_REGEN_ALL(), ion_style=ISS_AFE(), ion_frac=0.5)
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))

    m0 = 1000u"mₑ"
    n0 = 1e-6u"Å^-3"
    p.E_field = EFromETd(10u"Td", n0) * [0,0,1]

	gas = GAS(name=:gas,m0=m0,ρ=n0)
	ptype = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, ptype)

    p.gen_cf = (args...) -> CollFreqs(args..., ion_style, ion_frac)

    p.save_name = MakeSaveName("ConstIon" ; ionweight, gns, ion_style, ion_frac)

    p.meas_bins = (MEAS_BIN(:cum, :t), MEAS_BIN(:t))

    # p.t_grid = LogRange(1e5, 1e9, 101, prefactor=p.time_unit, as_val=true)
    p.t_grid = LogRange(1e5, 1e8, 101, prefactor=p.time_unit, as_val=true)
    p.max_substepsize = 1e7*p.time_unit

    p.ionisation_as_weight = convert(TypeBool, ionweight)
    p.generate_next_step_style = gns

    p
end

end
