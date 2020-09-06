
include("Argon.jl")
using .Argon

@CheckTurns for exponent in -4:3, coeff in 1:9
# @CheckTurns for exponent in [-1], coeff in [1]
	E = coeff * 10.0^exponent * u"Td"
	p = SetupParams(E)

	props = LoopMaxTime(p, 10)

	Save(p, props)
end
