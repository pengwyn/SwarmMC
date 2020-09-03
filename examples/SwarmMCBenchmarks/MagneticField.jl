module MagneticField

using SwarmMC

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MAXWELL(6.0u"Å^2*eV^(1/2)"))

    [elastic]
end
@xport function SetupParams(; BHx=1000Hx, Bθ=90, ETd=1u"Td", tmtr=293u"K")
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(1.0eV)))
               
    m0 = 4*amu
    ρ = 1e-6u"Å^-3"
    p.E_field = EFromETd(ETd, ρ) * [0,0,1]
    p.B_field = BFromBHx(BHx, ρ) * [sind(Bθ),0,cosd(Bθ)]

	gas = GAS(; name=:gas, m0, ρ, tmtr)
    electron = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)

    p.gen_cf = CollFreqs

    p.save_name = MakeSaveName("MagneticField" ; BHx, ETd, Bθ, tmtr)

    p.t_grid = LinRange(0,1e11,101)*SwarmMC.uT

    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:t), MEAS_BIN(:ss))

    p
end

end
