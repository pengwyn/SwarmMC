
module NessRobsonLoss

using DanUtils, Constants
@reexport using SwarmMC

function CollFreqs(params, gas, ptype, lossa, lossp)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MAXWELL(10u"Å^2*eV^(1/2)"))
    loss = CreateCollFreq(params, gas, ptype, "loss", CFS_LOSS(), GENCF_ATTACH(lossa, lossp))

    [elastic, loss]
end

@xport function SetupParams(;lossa=1e-3u"Å^2/sqrt(eV)",
                            lossp=0.5,
                            ans_style=ANS_NOTHING(),
                            gns_style=GNS_UPDATE_LOG2FAC(),
                            weight_reduction=0.5,
                            # split_fake=false,
                            tmtr=293.0u"K")

    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))

    m0 = 16*amu
    ρ = 1e-6u"Å^-3"
    p.E_field = EFromETd(0.4u"Td", ρ) * [0,0,1]

	gas = GAS(; m0, ρ, tmtr)
    electron = PARTTYPE(:electron)

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    
    p.gen_cf = (args...) -> CollFreqs(args..., lossa, lossp)

    p.generate_next_step_style = gns_style
    p.adapt_noncons_style = ans_style
    p.weight_reduction = weight_reduction
    # p.fake_noncons_split_weight = split_fake

    p.eps_grid = 1eV * LogRange(-5, 3, 1001)
    p.t_grid = LinRange(0, 5e10, 1001) * SwarmMC.uT
    p.min_log2_weight = -Inf

    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:ss, :cum), MEAS_BIN(:t))

    ansstr = string(typeof(ans_style).name.name)
    ansstr = ansstr[5:end]
    gnsstr = string(typeof(gns_style).name.name)
    gnsstr = gnsstr[5:end]

    p.save_name = MakeSaveName("NRL",
                               p=lossp,
                               a=ustrip(lossa),
                               ans=ansstr,
                               gns=gnsstr,
                               red=weight_reduction,
                               # split=split_fake,
                               T=ustrip(tmtr))

    p
end

end
