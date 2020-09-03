
using SwarmMC
using Generic

@everywhere using LucasSaelee

@CheckTurns for F in [0., 0.5, 1.0]
    p = SetupParams(F)
    p.meas_types = (MEAS_TYPE(:t), MEAS_TYPE(:cum,:t))
    p.init_style = PIS_SEPARABLE(PIVS_PUREZ(100.))
    p.t_grid = GRID(unique([0:1e-13:1e-12 ; 0:1e-11:2e-10]) / SwarmMC.CalcTimeUnit(p.mass_unit, p.len_unit, p.eps_unit))
    p.do_extra_meas = TypeTrue()
    p.save_name = MakeSaveName("LucasSaeleeTransient", F=F)

    p2 = FinaliseParams(p)
    props = LoopMaxTime(p2, 1000)

    Save(p2, props)
end
