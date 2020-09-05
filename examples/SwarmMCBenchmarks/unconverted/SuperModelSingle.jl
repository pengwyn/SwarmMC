
module SuperModelSingle

using Generic, Constants
using Reexport ; @reexport using SwarmMC

SuperModelElastic(eps) = (1. + 5*eps)
SuperModelInelastic(eps) = 1.

function CollFreqs(params, gas, ptype)
    Temp = params.chars[:temp, :default, :default]
    
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MANUAL(SuperModelElastic), temperature=Temp)
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_MANUAL(SuperModelInelastic), threshold=0.002, temperature=Temp)

    [elastic, inelastic]
end

@xport function SetupParams(ETd, Temp, grid_count=nothing)
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1)))

    m0 = 28 * amu / p.mass_unit
    n0 = 1e-6
    E = EFromETd(ETd, n0, p.len_unit)

	gas = GAS{:gas}(m0)
    electron = PARTTYPE{:electron}()
    region = REGION{:region}(E=E)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    push!(p.region_list, region)
    
    p.chars[:dens, :default, :default] = n0
    p.chars[:temp, :default, :default] = Temp
    p.chars[:cflist, :default, :default] = CollFreqs

    p.save_name = "SuperModelSingle:T=$(Temp):E=$(ETd)"

    if grid_count == nothing
        p.eps_grid = GRID(GRID_LOG(), 1e-5, 1e4, 1001)
    else
        p.eps_grid = GRID(GRID_LOG(), 1e-5, 1e4, grid_count)
        p.save_name *= ":gridsize=$(grid_count)"
    end
    
    #p.meas_types = (MEAS_TYPE(:cum,:t), MEAS_TYPE(:ss,:cum,:eps))

    p
end

end
