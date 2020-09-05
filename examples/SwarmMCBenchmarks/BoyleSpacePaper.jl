
module BoyleSpacePaper

using Generic, Constants, Reexport
@reexport using Reexport ; @reexport using SwarmMC
using BisectInterp
using PercusYevickSSF

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(CSmag))
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_HARDSPHERE(0.1), threshold=2.0)

    [elastic, inelastic]
end

const CSmag = 6.0

@xport function SetupParams(phi, init_eps_style=:delta, init_z_style=:delta)
    r,n0 = PercusYevickSSF.HardsphereRadiusDensity(CSmag, phi, default_phi=0.1, factor_4pi=false)
    println("Chose density of ",n0)

    zstar_conv = n0*1

    p = PARAMS()

    m0 = 4 * amu / p.mass_unit
    ETd = 3.
    E = EFromETd(ETd, n0, p.len_unit)
    
    tn0 = 1e-42 / p.time_unit / p.len_unit^3
    final_t = tn0 / n0
    p.t_grid = GRID(final_t)

    p.z_grid = GRID(-5 / zstar_conv, 40 / zstar_conv, 1001)
    z_kill = 80 / zstar_conv

	gas = GAS{:gas}(m0)
    electron = PARTTYPE{:electron}()
    region = REGION{:region}(E=E, right=z_kill, leave_event=REGION_EVENT_DIE())

    # TODO: This needs to changed into a drifted Maxwellian.
    if init_eps_style == :delta
        init_eps = PIVS_ISOTROPIC(1.0)
    elseif init_eps_style == :drifted_maxwellian
        init_eps = PIVS_MAXWELLIAN(p, electron.mass, 1e4, 1e5)
    else
        error("Not implemented: $init_eps_style.")
    end
    # TODO: The start spatial profile needs to be the Gaussian from the paper.
    if init_z_style == :delta
        init_z = PISS_DELTA()
    elseif init_z_style == :gaussian
        z_spread = 0.1 / zstar_conv
        init_z = PISS_GAUSSIAN(z_spread)
    else
        error("Not implemented: $init_z_style.")
    end

    p.init_style = PIS_SEPARABLE(init_eps, init_z)
        

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    push!(p.region_list, region)
    
    p.chars[:dens] = n0
    p.chars[:temp] = 0.
    p.chars[:cflist] = CollFreqs

    if phi > 0
        eps_grid = 10 .^ linspace(-5,3,10001)
        # Note: the Boyle paper used a wrong definition of the cross section in
        # r. It used σ=πr² not σ=4πr². To correct for that, we leave the cross
        # section at σ=6 but modify the apparent density in the SF.
        #SSF = PercusYevickIntegrate.(eps_grid, phi, rho0 / 4^(3/2) )
        # This has now been changed in my function HardsphereRadiusDensity above.
        SSF = PercusYevickIntegrate.(eps_grid, phi, n0)
        SSF_interp = BisectInterp.INTERP1D(eps_grid, SSF)

        p.static_structure_factor = SSF_interp
    end

    
    p.generate_next_step_style = GNS_NOTHING()
    p.meas_types = (p.meas_types..., MEAS_TYPE(:z, :cum))

    p.save_name = MakeSaveName("BoyleSpacePaper", Phi=phi, init_eps=init_eps_style, init_z=init_z_style)

    p
end

end
