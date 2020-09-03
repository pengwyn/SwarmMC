
module ReidAniso

using DanUtils, Constants
@reexport using SwarmMC
using BisectInterp

RiedInelastic(eps) = 0.4*(eps - 0.516eV)*u"Å^2/eV"

AnisoA(eps, R) = R*2 - 1
AnisoB(eps, R) = sign(2*R-1) * abs(2*R - 1.)^(1/5.)
AnisoC(eps, R) = -1/1.5 * log(1 - R*(1-exp(-3))) - 1
function AnisoD(eps, Rin)
    pointA = cos(3*pi/4.)
    pointB = cos(0.134*pi)

    widthA = pointA - -1
    widthB = 1 - pointB
    Rmax = widthA + widthB

    R = Rin * Rmax

    if R < widthA
        return (R/widthA) * widthA - 1
    else
        return (R - widthA)/widthB * widthB + pointB
    end
end

Aniso = Dict()
Aniso['A'] = AnisoA
Aniso['B'] = AnisoB
Aniso['C'] = AnisoC
Aniso['D'] = AnisoD

AnisoA2(eps, costh) = costh + 1
AnisoB2(eps, costh) = costh^5 + 1
AnisoC2(eps, costh) = -1/1.5 * (exp(-1.5*(costh+1)) - 1)
function AnisoD2(eps, costh)
    pointA = cos(3*pi/4.)
    pointB = cos(0.134*pi)

    widthA = pointA - -1
    widthB = 1 - pointB

    if costh < pointA
        return costh + 1
    elseif costh < pointB
        return pointA + 1
    else
        return pointA + (costh - pointB) + 1
    end
end

Aniso2 = Dict()
Aniso2['A'] = AnisoA2
Aniso2['B'] = AnisoB2
Aniso2['C'] = AnisoC2
Aniso2['D'] = AnisoD2






function CollFreqs(params, gas, ptype, cs, aniso_var)
    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_HARDSPHERE(cs), angle_dist_cum=aniso_var)
    inelastic = CreateCollFreq(params, gas, ptype, "inelastic", CFS_INELASTIC(), GENCF_MANUAL(RiedInelastic), threshold=0.516eV, angle_dist_cum=aniso_var)

    [elastic,inelastic]
end


using LegendrePolys
using QuadGK
using ForwardDiff

@xport function SetupParams(model::Union{Nothing,Char}, method=:func)
    @info "Preparing" model method
    p = PARAMS(init_style=PIS_SEPARABLE(PIVS_ISOTROPIC(0.1eV)))

	m0 = 2amu
    ρ = 1e17u"m^-3"
    ETd = 25u"Td"
    p.E_field = ETd * ρ * [0,0,1]

	gas = GAS(; m0, ρ)
    ptype = PARTTYPE(:ptype)

    push!(p.ptype_list, ptype)
    push!(p.gas_list, gas)

    if model === nothing
        aniso_func = nothing
    elseif model ∉ "ABCD"
        error("Unknown model '$model'.")
    else
        aniso_func = Aniso[model]
    end

    if model ∈ "CD"
        cs = 6.954u"Å^2"
    else
        cs = 10u"Å^2"
    end

    if method == :func
        aniso_var = aniso_func
    elseif model == nothing
        aniso_var = nothing
		method = :func
    elseif method == :grid
        epsrange = linspace(0,1,11)*eV
        costhrange = linspace(-1,1,1001)

        func = Aniso2[model]
        mat = [func(eps, costh) for costh=costhrange, eps=epsrange]

        #mat .+= mat[1:1,:]
        mat ./= mat[end:end,:]

        aniso_var = BisectInterp.INTERP2D(epsrange, costhrange, mat)
    elseif method == :legendre
        func = Aniso2[model]
        l_list = 0:2
        legendre_poly = LegendrePoly.(l_list)
        #sigma_l = [(2l+1)/4pi * quadgk(x -> ForwardDiff.derivative(x->func(1.,x), x)*legendre_poly[lind](x), -1.,1.)[1] for (lind,l) in enumerate(l_list)]
        sigma_l = [quadgk(x -> ForwardDiff.derivative(x->func(1.,x), x)*legendre_poly[lind](x), -1.,1.)[1] for (lind,l) in enumerate(l_list)]
        aniso_var = SwarmMC.MatchWithDeltas([1eV], sigma_l[newaxis,:])
    end

    p.gen_cf = (args...) -> CollFreqs(args..., cs, aniso_var)

    p.t_grid = LogRange(13, 17, 101) * SwarmMC.uT
    
    p.save_name = MakeSaveName("ReidAniso",model=model,method=method)

    p
end

end
