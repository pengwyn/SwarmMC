
module RegionTest

using DanUtils, Constants
@reexport using Reexport ; @reexport using SwarmMC

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(1.))

    [elastic]
end

@xport function SetupParams(;ETd=0.1, init_eps=0., bias=10., ETd_right=ETd, n0_right_factor=1., have_absorb=true)
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(init_eps)))

	m0 = ustrip(u"kg", 1amu) / p.mass_unit
    n0 = 1e-6
    Eleft = EFromETd(ETd, n0, p.len_unit)
    Eright = EFromETd(ETd_right, n0 * n0_right_factor, p.len_unit)

	gas = GAS{:gas}(m0)
    electron = PARTTYPE{:electron}()

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)

    region_left = REGION{:left}(E=Eleft, V=bias, right=6e7)
    push!(p.region_list, region_left)
    if have_absorb
        region_right = REGION{:right}(E=Eright, V=0., left=6e7, right=10e7)
        region_absorb = REGION{:absorb}(E=0., V=0., left=10e7, enter_event=REGION_EVENT_MEASURE((:eps,:t)))

        push!(p.region_list, region_right)
        push!(p.region_list, region_absorb)
    else
        region_right = REGION{:right}(E=Eright, V=0., left=6e7)

        push!(p.region_list, region_right)
    end
    
    p.chars[:dens, :default, :default] = n0
    p.chars[:dens, :right, :default] = n0 * n0_right_factor
    p.chars[:temp, :default, :default] = 0.
    p.chars[:cflist, :gas, :electron] = CollFreqs

    p.save_name = MakeSaveName("RegionTest", ETd=ETd, init_eps=init_eps, separate_dir=false)

    p.meas_types = (MEAS_TYPE(:t,:cum), MEAS_TYPE(:cum,:z))

    p.t_grid = GRID(GRID_LOG(), 1e5, 1e11)
    p.z_grid = GRID(-5e7, 11e7, 1001)

    p.adapt_noncons_style = ANS_NOTHING()
    
    p
end


using Plots

function DoATimeProfileRun(N=10000)

    params = SetupParams(ETd=0.5, bias=1.)
    params.meas_types = (MEAS_TYPE(:cum,:z), MEAS_TYPE(:cum, :z, :t))

    params = FinaliseParams(params)
    props = BunchedPropagate(params, N)

    profiles = props[:duration, (:cum,:z,:t)]
    dt = diff(params.t_grid.grid)
    profiles ./= dt

    z = ZGrid(params)
    ymax = maximum(profiles[50:end,:])
    anim = @animate for i = 1:length(TGrid(params))
       plot(z[2:end-1], profiles[i,2:end-1], ylim=(0,ymax))
    end

    #gif(anim, "test.gif")
    mp4(anim, "test.mp4")

    return anim
end

end
