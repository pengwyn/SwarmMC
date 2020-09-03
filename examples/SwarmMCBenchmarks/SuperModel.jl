
module SuperModel

using DanUtils, Constants
@reexport using SwarmMC

SuperModelElastic(eps) = (1u"Å^2" + 5u"Å^2/eV"*eps)
SuperModelInelastic(eps) = 1u"Å^2"

function CollFreqs(params, gas, ptype)
    if gas.name == :elastic
        elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MANUAL(SuperModelElastic))
        [elastic]
    elseif gas.name == :inelastic
        inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_MANUAL(SuperModelInelastic), threshold=0.002eV)
        [inelastic]
    end
end

@xport function SetupParams(ETd=0.0u"Td", elastic_T=293.0u"K", inelastic_T=293.0u"K", inelastic_recoil=true)
    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(0.1eV)))

	m0 = 28amu
    ρ = 1e-6u"Å^-3"
    p.E_field = ETd * ρ * [0,0,1]

	elastic_gas = GAS(; name=:elastic, m0, ρ, tmtr=elastic_T)

    if inelastic_recoil
        inel_m0 = m0
    else
        inel_m0 = Inf*mₑ
    end
    inelastic_gas = GAS(; name=:inelastic, m0=inel_m0, ρ, tmtr=inelastic_T)

    electron = PARTTYPE(:electron)

    push!(p.gas_list, elastic_gas)
    push!(p.gas_list, inelastic_gas)
    push!(p.ptype_list, electron)
    
    p.gen_cf = CollFreqs

    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:ss,:cum,:eps))
    p.t_grid = LinRange(0,1e11,101)*SwarmMC.uT
    
    p.save_name = MakeSaveName("SuperModel", elT=elastic_T, inelT=inelastic_T, E=ETd, inelrecoil=inelastic_recoil)

    p
end

end
