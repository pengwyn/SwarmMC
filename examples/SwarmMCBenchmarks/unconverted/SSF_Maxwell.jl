
using Reexport ; @reexport using SwarmMC

import BisectInterp

function SetupParams(S=1.0, ETd=1.0)
    p = PARAMS()
    m0 = 4*SwarmMC.units_amu/p.mass_unit
    p.ETd = ETd
    p.init_eps = 0.1
    p.max_eps = 1000.

    a = 6.
    n0 = 1e-6
	gas = GAS("gas", m0, n0)

    elastic = COLLFREQ(CFS_ELASTIC(), "elastic", gas, ("Maxwell", a))
    push!(p.cf_list, elastic)

    p.static_structure_factor = @eval eps -> $S

    #if phi > 0
    #    eps_grid = logspace(-5,3,10001)
    #    #eps_grid = linspace(0, 100, 1001)
    #    SSF = PercusYevickIntegrate.(eps_grid, phi, n0)
    #    SSF_interp = BisectInterp.INTERP1D(eps_grid, SSF)

    #    p.static_structure_factor = SSF_interp
    #end

    p.save_name = "SSFMaxwell_S$(S)_E$(ETd)"

    p
end
