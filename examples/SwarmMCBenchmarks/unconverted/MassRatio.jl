
module MassRatio

export SetupParams

using SwarmMC
using Constants

function SetupParams(Temp, ETd, mratio=1e-4)
    p = PARAMS(ETd=ETd,
               init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1)))

    m0 = 4 * amu/p.mass_unit
    m = m0 * mratio

	gas = GAS("gas", m0, 1e-6, Temp)
    ptype = PARTTYPE("mratio", m)

    elastic = COLLFREQ(CFS_ELASTIC(), "elastic", gas, ptype, ("Hardsphere", 6.))

    push!(ptype.cf_list, elastic)
    push!(p.ptype_list, ptype)

    p.measure_time_dists = :yes
    
    p.save_name = "MassRatio_T$(Temp)_mratio$(round(mratio,sigdigits=2))_E$ETd"

    p.final_time = 1e10

    p
end

end
