__precompile__()

module CurrentTrace

using Generic, Constants, Reexport
@reexport using Reexport ; @reexport using SwarmMC

const Temp = 27.



#TrappingCollFreq(eps, eps_threshold, mag) = (eps < 1.0 * kB*Temp / eV) * 20.0
TrappingCollFreq(eps, eps_threshold, mag) = (eps < eps_threshold) * mag

function CollFreqs(params, gas, ptype, solvation_as_loss, temperature, trap_mag)
    if GetSymbol(ptype) == :trapped
        return []
    end
    
    neon_cs = readcsv("Neon_LandoltBornstein.csv", skipstart=1)
    neon_cs[:,1] *= eV / params.eps_unit
    neon_cs[:,2] *= (1e-10 / params.len_unit)^2

    cf_list = COLLFREQ[]
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_INTERP(neon_cs[:,1], neon_cs[:,2]), temperature=temperature)
    push!(cf_list, elastic)

    # eps_threshold = kB*temperature/eV
    # trapfunc = let eps_threshold=eps_threshold, trap_mag=trap_mag
    #     eps -> TrappingCollFreq(eps, eps_threshold, trap_mag)
    # end

    epsstar = 37.29 * kB
    trap_x = epsstar * [0., 0.00010081, 0.00020922, 0.0005721, 0.003461, 0.010102, 0.04541, 0.11103, 0.22551, 0.2743, 0.4139]
    trap_y = 1e12 * [3000., 2960, 2674, 1761, 650.8, 215.5, 27.58, 4.123, 0.9915, 0.12596, 0.011061]

    trap_x /= params.eps_unit
    trap_y /= 1/params.time_unit

    trap_x *= trap_mag
    trap_y /= TotDensity(params)

    gencf_trap = GENCF_INTERPLOG(trap_x, trap_y)

    if solvation_as_loss
        loss = CreateCollFreq(params, gas, ptype, "loss", CFS_LOSS(), gencf_trap, is_cross_section=false)
        push!(cf_list, loss)
    else
        capture = CreateCollFreq(params, gas, ptype, "capture", CFS_CHANGETYPE(), gencf_trap, new_ptype=FindPtypeByName(params, :trapped), is_cross_section=false)
        push!(cf_list, capture)
    end

    return cf_list
end

@xport function SetupParams(; solvation_as_loss=false, delay_type=:fixed_dt, delay_median=1e6, trap_mag=20., Eorig=92.3, cell_width_SI=2e-3, return_eps=0.)
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_PUREZ(0.3)))

    # Using a Sakai value
    #E = 92.3 # kV/cm
    E = Eorig * 1000 * 100
    E *= p.len_unit
    
    rho0 = 0.036019 # For 1.207 g/cm

    m0 = 20.1797*amu / p.mass_unit

	neon = GAS{:neon}(m0)
    quasifree = PARTTYPE{:quasifree}()

    push!(p.gas_list, neon)
    push!(p.ptype_list, quasifree)

    if solvation_as_loss
        p.save_name = MakeSaveName("CurrentTrace", AsLoss=solvation_as_loss, E=Eorig, width=cell_width_SI, trap_mag=trap_mag)
    else
        trapped = PARTTYPE{:trapped}(Inf)
        push!(p.ptype_list, trapped)

        let median=delay_median, freq=1/delay_median

            if delay_type == :fixed_dt
                trapped.delay_lookup = () -> median
            elseif delay_type == :exp
                trapped.delay_lookup = () -> median * -log(rand())
            elseif delay_type == :fractional
                # This is an infinite mean but a median of 3/freq
                # P(t) = 1/(freq*t + 1)^3/2
                trapped.delay_lookup = () -> 1/freq * (1/rand()^2 - 1)
            else
                error("Unknown delay type $delay_type")
            end
        end
        
        if return_eps == :thermal
            return_eps_val = let w = wFromT(quasifree, Temp, p.eps_unit)
                () -> SampleFromMaxwellian(w)
            end
        else
            return_eps_val = return_eps
        end
        trapped.delay_event = DELAY_EVENT_CHANGETYPE(quasifree, return_eps_val)

        p.save_name = MakeSaveName("CurrentTrace", AsLoss=solvation_as_loss, E=Eorig, width=cell_width_SI, delay_type=string(delay_type), delay_median=delay_median, trap_mag=trap_mag, return_eps=return_eps)
    end

    cell_width = cell_width_SI / p.len_unit
    
    region_left = REGION{:left}(right=-0.1, enter_event=REGION_EVENT_DIE())
    region_right = REGION{:right}(left=cell_width, enter_event=REGION_EVENT_DIE())
    region_absorb = REGION{:cell}(E=E, left=-0.1, right=cell_width)

    push!(p.region_list, region_left)
    push!(p.region_list, region_right)
    push!(p.region_list, region_absorb)

    p.chars[:dens, :default, :default] = rho0
    p.chars[:temp, :default, :default] = Temp
    p.chars[:cflist, :neon, :default] = (args...) -> CollFreqs(args..., solvation_as_loss, Temp, trap_mag)

    p.meas_types = (MEAS_TYPE(:t,:cum), MEAS_TYPE(:cum,:z))

    #final_time = 1e-6 / p.time_unit
    #p.t_grid = GRID(GRID_LOG(), final_time*1e-4, final_time, 1001)
    p.t_grid = GRID(GRID_LOG(), 1e-9 / p.time_unit, 1e-5 / p.time_unit, 1001)
    p.z_grid = GRID(GRID_LINEAR(), 0., cell_width, 1001)

    p.generate_next_step_style = GNS_DOUBLE()
    p.adapt_noncons_style = ANS_NOTHING()

    p
end


function ShowPlotsSlowly(;plot_type=:current, mode=:pausing, restrict=nothing, num_bins=1)
    prefixes = PrefixSets(r"^CurrentTrace") |> keys |> collect |> sort

    if restrict != nothing
        prefixes = filter!(x -> ismatch(restrict, x), prefixes)
    end

    @eval using PyPlot
    figure(figsize=(4,3))
    for prefix in prefixes
        params,props = ReadAll(prefix)

        dt = diff(params.t_grid.grid)
        if plot_type == :current
            yval = props[(:vel,:z)].*props[:duration]./dt ./ props.num_particles
        elseif plot_type == :duration
            yval = props[:duration] ./ dt ./ props.num_particles
        else
            yval = props[plot_type]
        end
        time = TGrid(params) * params.time_unit / 1e-6

        if num_bins != 1
            time = Binned(time, num_bins)
            yval = Binned(yval, num_bins)
        end

        plot(time, yval)
        xlabel(raw"Time ($\mu$s)")
        ylabel("Current")


        if mode == :pausing
            readline(STDIN)
        elseif mode == :saving
            savefig(prefix * "_fig.png")
            xscale("log")
            yscale("log")
            plot_val[isnan(plot_val)] = 0.
            ymax = maximum(plot_val)
            ymin = minimum(plot_val)
            ymin = max(ymin, ymax*1e-8)
            ylim(ymin, ymax)
            
            savefig(prefix * "_fig_loglog.png")
            close()
        elseif mode == :separate
            title(prefix, fontsize=4)
            #ylim(ymin=0)
            xscale("log")
            yscale("log")
            tight_layout()
            figure(figsize=(4,3))
        else
            error("Shouldn't get here")
        end
    end
end

end
