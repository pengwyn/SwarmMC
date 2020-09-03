
module CoulombPositron

# Basing this off the values in Bobylev and Potapenko (JCP 246, 123, 2013)

using DanUtils, Constants, Reexport
@reexport using SwarmMC

function CollFreqs(params, gas, ptype ; ϵ, logL)
    m = ptype.mass
    m0 = gas.m0
    μ = m*m0 / (m + m0)

    gencf,angle_func = SetupCoulombGenCF(ϵ, logL, μ, 1, 1)
    elastic = CreateCollFreq(params, gas, ptype, "coulomb", CFS_ELASTIC(), gencf, angle_dist_cum=angle_func, is_cross_section=false)

    [elastic]
end

@xport function SetupParams(Temp, ϵ, logL)
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(3000.0)))

    # Hydrogen??? This doesn't seem right....
    m0 = 64
    # This is Yaniss' choice (scaled out in Bobylev)
    n0 = ustrip(u"m^-3", 1e-5*amagat) / p.len_unit^-3

	electron = GAS{:electron}(1.0)
    ion = GAS{:ion}(m0)
    positron = PARTTYPE{:positron}()
    region = REGION{:region}()

    push!(p.ptype_list, positron)
    push!(p.gas_list, electron)
    push!(p.gas_list, ion)
    push!(p.region_list, region)

    p.chars[:dens] = n0
    p.chars[:temp] = Temp
    p.chars[:cflist] = (args...) -> CollFreqs(args..., ϵ=ϵ, logL=logL)

    p.save_name = MakeSaveName("CoulombPositron", T=Temp)

    p.t_grid = GRID(GRID_LOG(), 1e7, 1e11, 1001)
    # p.t_grid = GRID(0., 1e11, 1001)

    p
end

end
