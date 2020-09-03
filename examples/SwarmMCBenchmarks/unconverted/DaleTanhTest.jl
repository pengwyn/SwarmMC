
module DaleTanhTest

using Generic
@reexport using SwarmMC
using Constants

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(5e-20 / params.len_unit^2))

    [elastic]
end

@xport function SetupParams()
    # zstar_unit = 1 / (1e23 * 1e-20)
    # p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.)), len_unit=zstar_unit)

    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.)))

    ETd = 1. #Td
    
    m0 = 4*amu/p.mass_unit

	gas = GAS{:gas}(m0)
    electron = PARTTYPE{:electron}()

    push!(p.ptype_list, electron)
    push!(p.gas_list, gas)

    p.chars[:temp] = 0.
    p.chars[:cflist, :gas, :electron] = CollFreqs

    dens_left = 1e23 * p.len_unit^3
    dens_right = 1e24 * p.len_unit^3

    zstar_conv = 1 / (dens_left * (1e-20 / p.len_unit^2))
    # zstar_conv = 1.

    z_final = 300 * zstar_conv
    p.z_grid = GRID(GRID_LINEAR(), 0, z_final, 101)

    tanh_mid = 100. * zstar_conv
    tanh_spread = 25. * zstar_conv
    
    z_mid = RunningAvg(p.z_grid)
    dens_grid = @. dens_left + (dens_right-dens_left)/2*(1 + tanh(2*(z_mid-tanh_mid) / tanh_spread))
    @show @. 1/2*(1 + tanh(2*(z_mid-tanh_mid) / tanh_spread))

    EFromN(N) = EFromETd(ETd, N, p.len_unit)
    
    push!(p.region_list, REGION{:left}(E=EFromN(dens_left), right=p.z_grid[1]))
    p.chars[:dens, :left, :gas] = dens_left
    for ind = 1:length(z_mid)
        sym = Symbol("mid" * string(ind))
        push!(p.region_list, REGION{sym}(E=EFromN(dens_grid[ind]), left=p.z_grid[ind], right=p.z_grid[ind+1]))
        p.chars[:dens, sym, :gas] = dens_grid[ind]
    end
    push!(p.region_list, REGION{:right}(E=EFromN(dens_right), left=p.z_grid[end]))
    p.chars[:dens, :right, :gas] = dens_right

    p.meas_types = (MEAS_TYPE(:cum,:t), MEAS_TYPE(:cum,:z))

    p.t_grid = GRID(5e-5 / p.time_unit)
    
    p.save_name = MakeSaveName("DaleTanh")

    p
end

end
