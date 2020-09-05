
module PYStepModel

using Reexport ; @reexport using SwarmMC
using Generic, Constants

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(6.))
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_HARDSPHERE(0.1), threshold=2.0)

    [elastic, inelastic]
end

@xport function SetupParams(;ETd=3., phi=0.)
    @assert phi == 0.

    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1)))

    m0 = 4*amu/p.mass_unit
    n0 = 1e-6
    E = EFromETd(ETd, n0, p.len_unit)

	gas = GAS{:gas}(m0)
    electron = PARTTYPE{:electron}()
    region = REGION{:region}(E=E)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    push!(p.region_list, region)
    
    p.chars[:dens] = n0
    p.chars[:temp] = 0.
    p.chars[:cflist] = CollFreqs

    tn0 = 1e-42 / p.time_unit / p.len_unit^3
    final_t = tn0 / n0

    p.t_grid = GRID(final_t)

    p.meas_types = (MEAS_TYPE(:t, :cum), MEAS_TYPE(:eps, :ss, :cum))

    p.save_name = MakeSaveName("PYStepModel", E=ETd, phi=phi)

    p
end

end
