
module DaleTanhTest

using Reexport ; @reexport using SwarmMC

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(5u"Å^2"))

    [elastic]
end

@xport function SetupParams()
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1eV)))

    ρ_left = 1e-7u"Å^-3"
    ρ_right = 1e-6u"Å^-3"

    zstar = 1 / (ρ_left * u"Å^2")
    @show zstar

    z_final = 300 * zstar
    p.z_grid = LinRange(0u"Å", z_final, 101)

    tanh_mid = 100. * zstar
    tanh_spread = 25. * zstar

    ρ_func(pos) = ρ_left + (ρ_right - ρ_left)/2 * (1 + tanh(2*(pos[3] - tanh_mid) / tanh_spread))

    # p.E_field = EFromETd(1Td,ρ_func,SwarmMC.XYZ([0,0,1]))
    # p.E_field = EFromETd(1Td,ρ_left) * [0,0,1]
    p.E_field = EFromETd(1Td,ρ_right) * [0,0,1]
    
    m0 = 4*amu

	gas = GAS(; m0, ρ=ρ_func)
    electron = PARTTYPE(:electron)

    push!(p.ptype_list, electron)
    push!(p.gas_list, gas)

    p.gen_cf = CollFreqs

    # What was the original benchmark? Same ETd everywhere?
    # EFromN(N) = EFromETd(ETd, N, p.len_unit)
    
    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:cum,:z))

    p.t_grid = LinRange(0,5e11,101) * p.time_unit
    # p.t_grid = LinRange(0,5e1,101) * p.time_unit
    
    p.save_name = MakeSaveName("DaleTanh")

    p
end

end
