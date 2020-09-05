
module Hardsphere

using Reexport ; @reexport using SwarmMC

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(5u"Å^2"))

    [elastic]
end

@xport function SetupParams(tmtr ; m0=4amu, ETd=1Td)
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))
    n0 = 1e-6u"Å^-3"
    E = EFromETd(ETd, n0)

	gas = GAS(name=:gas,m0=m0,ρ=n0,tmtr=tmtr)
    electron = PARTTYPE(:electron)

    push!(p.ptype_list, electron)
    push!(p.gas_list, gas)

    p.E_field = E * [0,0,1]

    p.gen_cf = CollFreqs

    p.save_name = MakeSaveName("Hardspheremodel", T=tmtr, ETd=ETd)

    # p.t_grid = LinRange(0, 1e11, 101) * SwarmMC.uT
    p.t_grid = LinRange(0, 1e10, 101) * SwarmMC.uT

    # p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:t), MEAS_BIN(:ss,:cum,:eps), MEAS_BIN(:cum, :t, :eps))
    # p.meas_bins = (MEAS_BIN(:cum,:t),)
    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:cum, :ss, :eps))
    # p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:t))
    # p.meas_bins = (MEAS_BIN(:t),)
    # p.meas_bins = (MEAS_BIN(:t), MEAS_BIN(:ss,:eps),)

    # p.meas_quants = SwarmMC.Q_quick_list
    p.meas_quants = SwarmMC.Q_basic_list

    p
end

end
