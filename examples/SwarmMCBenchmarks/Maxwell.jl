
module Maxwell

using DanUtils, Constants
@reexport using SwarmMC

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MAXWELL(5.0u"Å^2*eV^(1/2)"))

    [elastic]
end

@xport function SetupParams(tmtr=293u"K", ETd=1u"Td")
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))
    m0 = 4amu
    n0 = 1e-6u"Å^-3"
    if ETd !== nothing
        p.E_field = EFromETd(ETd, n0) * [0,0,1]
    end

	gas = GAS(; name=:gas,m0,ρ=n0,tmtr)
    electron = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    
    p.gen_cf = CollFreqs

    p.save_name = MakeSaveName("Maxwellmodel",T=tmtr,E=ETd)

    p.t_grid = LinRange(0, 1e11, 101) * SwarmMC.uT

    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:t), MEAS_BIN(:cum, :ss, :eps))

    p
end

end
