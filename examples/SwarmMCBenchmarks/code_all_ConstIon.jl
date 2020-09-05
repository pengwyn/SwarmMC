
union!(LOAD_PATH, ["."])
using ConstIon

@CheckTurns for gns in [GNS_UPDATE_LOG2FAC(), GNS_NOTHING(), GNS_REGEN_ALL(), GNS_DOUBLE()],
                ionweight in [true,false],
                ion_style in [ISS_AFE(), ISS_FRACTION()],
                ion_frac in [0.5, 0.8]

    ionweight && (gns ∈ [GNS_UPDATE_LOG2FAC(), GNS_NOTHING(), GNS_REGEN_ALL()] || continue)
    ion_frac != 0.5 && (ion_style == ISS_FRACTION() || continue)

	p = SetupParams(; ionweight, ion_style, ion_frac, gns)

	props = LoopMaxTime(p, 10)

	Save(p, props)
end
