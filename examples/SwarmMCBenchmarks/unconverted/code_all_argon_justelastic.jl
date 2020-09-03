
using DanUtils
using Argon

@CheckTurns for exponent in -4:3, coeff in 1:9
# @CheckTurns for exponent in [-1], coeff in [1]
	E = coeff * 10.0^exponent
	p = SetupParams(E, just_elastic=true, Temp=0.)

	props = LoopMaxTime(p, 1)

	Save(p, props)
end
