
module MaxwellDaleTest

export SetupParams

using SwarmMC
using Constants

function SetupParams(ETd=1., final_time=1e-6)
    p = PARAMS(ETd=ETd,
               init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1)),
               num_times=1001)
    m0 = 4*amu/p.mass_unit / 100.

	gas = GAS("gas", m0, 1e-7, 0.)
    ptype = PARTTYPE()

    elastic = COLLFREQ(CFS_ELASTIC(), "elastic", gas, ptype, ("Maxwell", 5.))

    push!(ptype.cf_list, elastic)
    push!(p.ptype_list, ptype)

    p.save_name = "MaxwellDaleTest"

    #p.final_time = 1e11
    p.final_time = final_time/p.time_unit

    p
end

end
