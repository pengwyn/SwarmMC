
include("LucasSaelee.jl")
using .LucasSaelee

@CheckTurns for F in [0., 0.25, 0.5, 0.75, 1.0],
    ion_weight = [false,true]

	p = SetupParams(; F, ion_weight)

	props = LoopMaxTime(p, 100)

	Save(p, props)
end


@CheckTurns for F in [0., 0.5, 1.0],
    BHx in [100,200,1000]*Hx,
    Bθ in [0, 30, 60, 90],
    ion_weight = [false,true]

	p = SetupParams(; F, BHx, Bθ)

	props = LoopMaxTime(p, 100)

	Save(p, props)
end
