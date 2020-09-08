
module Hardsphere_stoptest

using Reexport ; @reexport using SwarmMC

function CollFreqs(params, gas, ptype)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(5u"Å^2"))

    [elastic]
end

@xport function SetupParams(tmtr ; m0=4amu, ETd=1Td)
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))

    n0 = 1e-6u"Å^-3"
    p.E_field = EFromETd(ETd, n0) * [0,0,1]

	gas = GAS(; m0, ρ=n0, tmtr)
    electron = PARTTYPE(:electron)

    push!(p.ptype_list, electron)
    push!(p.gas_list, gas)

    p.gen_cf = CollFreqs

    p.save_name = MakeSaveName("Hardspheremodel_asldjflskdjf", T=tmtr, ETd=ETd)

    p.t_grid = LinRange(0, 1e10, 101) * SwarmMC.uT

    p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:t), MEAS_BIN(:ss,:cum,:eps), MEAS_BIN(:cum, :t, :eps))
    # p.meas_bins = (MEAS_BIN(:cum,:t), MEAS_BIN(:cum, :ss, :eps))

    p.user_callbacks = [MyCallback]
    
    p
end



using UnPack

function TestCond(u,t,int)
    # @unpack pos = UnpackVector(u)
    t - 1e9
end
function TestAffect!(int)
    @unpack part,params,extra_part_list = int.p
    part.weight = 0*SwarmMC.uN
    terminate!(int)
end
using OrdinaryDiffEq
MyCallback = ContinuousCallback(TestCond, TestAffect!, nothing, abstol=0.0)



end
