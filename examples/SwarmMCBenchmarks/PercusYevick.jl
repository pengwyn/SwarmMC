
module PercusYevick

using Reexport ; @reexport using SwarmMC
using Dierckx#, UnitfulDierckxInternal

using PercusYevickSSF

const a = 6u"Å^2"
function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(a))
    [elastic]
end

@xport function SetupParams(ETd=3u"Td", phi=0.4)
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1eV)))

    m0 = 4amu

    # This is missing a 4 but this is what the benchmark paper (Tattersall et al 2015) used.
    # Actually I think this is correct because at high energy it is cs = pi*r^2 not cs = 4pi*r^2 as it is at zero energy.
    r,ρ = PercusYevickSSF.HardsphereRadiusDensity(a, phi, factor_4pi=false)
    p.E_field = ETd * ρ * [0,0,1]
    
	gas = GAS(; m0,ρ)
	ptype = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, ptype)

    if phi > 0
        eps_grid = LogRange(-5,3,10001)*eV
        #eps_grid = LogRange(-5,3,101)*eV
        SSF = PercusYevickIntegrate.(eps_grid, phi, ρ)
        SSF_interp = Spline1D(eps_grid, SSF)

        p.static_structure_factor = SSF_interp
    end

    p.gen_cf = CollFreqs

    if ETd > 1Td
        final_time = 1e6
    elseif ETd > 0.1Td
        final_time = 1e7
    else
        final_time = 1e7
    end

    p.t_grid = LinRange(0, final_time, 101) * SwarmMC.uT

    p.save_name = MakeSaveName("PercusYevick", Phi=phi, E=ETd)

    p
end

end
