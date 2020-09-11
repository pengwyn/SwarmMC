
function SetupAngleDistCum(costh, val)
    dx = diff(costh)
    dy = diff(val)

    cumval = [zeros(size(val, 2)) ; cumsum(val[1:end-1].*dx, 2)]
    m = dy./dx
    cumval[2:end] += cumsum(0.5*m.*dx.^2, 2)

    cumval
end

abstract type GENCF_PARAMS end

##############################
# * GENCF_MANUAL
#----------------------------
@xport struct GENCF_MANUAL <: GENCF_PARAMS
    func::Function
    prefac::Float64
    GENCF_MANUAL(func, prefac=1.) = new(func, prefac)
end
    
function GenCF(gencf::GENCF_MANUAL)
    if gencf.prefac != 1.
        let prefac=gencf.prefac, func=gencf.func
            thisfunc = eps -> prefac * func(eps)
        end
    else
        thisfunc = gencf.func
    end

    return thisfunc
end

##############################
# * GENCF_ATTACH
#----------------------------
@xport struct GENCF_ATTACH <: GENCF_PARAMS
    a
    p
end
    
function GenCF(gencf::GENCF_ATTACH)
    #p = gencf.p + 0.5
    p = gencf.p
    if p == 1
        sub_expr = eps -> eps
    elseif p == 0.5
        sub_expr = eps -> sqrt(eps)
    elseif p == 0.
        sub_expr = eps -> 1.0
    elseif p == -0.5
        sub_expr = eps -> 1/sqrt(eps)
    elseif p == -1.0
        sub_expr = eps -> 1/eps
    else
        println("Warning! This is a generic version for the attach cf!")
        sub_expr = eps -> eps^p
    end

    local func
    let a=gencf.a, sub_expr=sub_expr
        func = eps -> a * sub_expr(eps)
    end
end

##############################
# * GENCF_MAXWELL
#----------------------------

@xport struct GENCF_MAXWELL <: GENCF_PARAMS
    a::Q(uσ * sqrt(uE))
end

function GenCF(gencf::GENCF_MAXWELL)
    local func
    let a=gencf.a
        func = eps -> a / sqrt(eps)
    end

    func
end

##############################
# * GENCF_HARDSPHERE
#----------------------------

@xport struct GENCF_HARDSPHERE <: GENCF_PARAMS
    a::Q(uσ)
end

function GenCF(gencf::GENCF_HARDSPHERE)
    local func
    let a=gencf.a
        func = eps -> a
    end

    func
end

##############################
# * GENCF_FAKENONCONS
#----------------------------

@xport struct GENCF_FAKENONCONS <: GENCF_PARAMS
    params::PARAMS
end

function GenCF(gencf::GENCF_FAKENONCONS)
    let p=gencf.params
        func = eps -> abs(p.fake_noncons_rate)
    end
end

##############################
# * GENCF_INTERP
#----------------------------
@xport struct GENCF_INTERP{T} <: GENCF_PARAMS
    x::Vector{Q(uE)}
    y::Vector{T}
    extrap
end
GENCF_INTERP(x,y) = GENCF_INTERP(x,y,"extrapolate")

function GenCF(gencf::GENCF_INTERP)
    interpobj = Spline1D(gencf.x, gencf.y, k=1, bc=gencf.extrap)
end

@xport struct GENCF_INTERPLOG <: GENCF_PARAMS
    x::Vector{Float64}
    y::Vector{Float64}
    extrap
end
GENCF_INTERPLOG(x,y) = GENCF_INTERPLOG(x,y,"extrapolate")

function GenCF(gencf::GENCF_INTERPLOG)
    @unpack x,y,extrap = gencf

    if any(y .<= 0)
        error("Can't have any points zero or negative for a InterpLog!")
    end

    interpobj = Spline1D(log.(x),log.(y), k=1, bc=extrap)

    local func
    if extrap == "zero"
        let interpobj=interpobj, x0=x[1], xend=x[end]
            func = eps -> (eps < x0 || eps > xend) ? 0. : exp(interpobj(log(eps)))
        end
    else
        let interpobj=interpobj
            func = eps -> exp(interpobj(log(eps)))
        end
    end

    func
end


################################
# * Coulomb collisions - not a new GENCF
#------------------------------

function CoulombLogarithm(ρ_e, T_e, u=:temperature, μ=(u == :temperature) ? nothing : error("Need a reduced mass!"), q1=1, q2=1)
    λ_D = sqrt(ustrip(u"q^2/eV/Å",eps0) * ustrip(u"eV/K",kB) * T_e / ρ_e / ustrip(u"q",1echarge)^2)

    b0 = q1*q2 * ustrip(u"q", 1echarge)^2  / (4π*ustrip(u"q^2/eV/Å",eps0))
    if u == :temperature
        b0 /= 3/2 * ustrip(u"eV/K",kB)*T_e
    else
        b0 /= 0.5 * μ * u^2
    end

    logL = log(sqrt(1 + λ_D^2/b0^2))
end

@xport function SetupCoulombCFFixedLogL(params, gas, ptype, ϵ, logL, q1=1, q2=1 ; name="coulomb")
    m = ptype.mass
    m0 = gas.m0
    μ = m*m0 / (m + m0)

    # This needs to be used with is_cross_section=false.
    gencf = GENCF_HARDSPHERE(1/ϵ)

    # In units of eV, Å, e we have:
    eps0 = 0.005526349361053035
    # Assuming the above units and that μ is in electron masses
    C = 4π * q1^2*q2^2 * logL / ( (4π*eps0)^2 * μ^2 )
    # Also going to factor out the terms for v, assuming we are an electron
    C /= (2/1)^(3/2)

    angle_change = COULOMB_ANGLE(ϵ * C)

    cf = CreateCollFreq(params, gas, ptype, name, CFS_ELASTIC(), gencf, angle_dist_cum=angle_change, is_cross_section=false)
end

struct COULOMB_ANGLE_FIXED_LOGL
    C::Float64
end

function SelectCollCosTheta(angle_change::COULOMB_ANGLE_FIXED_LOGL, eps)
    # TODO: I should really recalc the logL at each energy here. Particularly relevant when we're at high energies.

    v_like = sqrt(eps)

    temp = angle_change.C / v_like^3

    temp = ClampMax(temp, 2.0)

    costheta = 1 - temp
end

@xport function SetupCoulombCF(params, gas, ptype, ϵ, T_e, ρ_e, q1=1, q2=1 ; name="coulomb")
    m = ptype.mass
    m0 = gas.m0
    μ = m*m0 / (m + m0)

    # This needs to be used with is_cross_section=false.
    gencf = GENCF_HARDSPHERE(1/ϵ)

    angle_change = COULOMB_ANGLE(ϵ, q1, q2, μ, T_e, ρ_e)

    cf = CreateCollFreq(params, gas, ptype, name, CFS_ELASTIC(), gencf, angle_dist_cum=angle_change, is_cross_section=false)
end

struct COULOMB_ANGLE
    ϵ::Float64
    q1::Float64
    q2::Float64
    μ::Float64
    T_e::Float64
    ρ_e::Float64
end

function SelectCollCosTheta(angle_change::COULOMB_ANGLE, eps)
    @unpack ϵ,ρ_e,T_e,μ,q1,q2 = angle_change

    # In units of eV, Å, e we have:
    eps0 = 0.005526349361053035
    # Assuming the above units and that μ is in electron masses
    C = 4π * q1^2*q2^2 / ( (4π*eps0)^2 * μ^2 )

    u = sqrt(2*eps)
    logL = CoulombLogarithm(ρ_e, T_e, u, μ)

    temp = ϵ * C * logL / u^3

    temp = ClampMax(temp, 2.0)

    costheta = 1 - temp
end

@xport function SetupCoulombCFVarying(params, gas, ptype, target_angle_change, T_e, ρ_e, q1=1, q2=1 ; name="coulomb")
    m = ptype.mass
    m0 = gas.m0
    μ = m*m0 / (m + m0)

    eps0 = 0.005526349361053035
    C = 4π * q1^2*q2^2 / ( (4π*eps0)^2 * μ^2 )

    logL = CoulombLogarithm(ρ_e, T_e)
    cf = C * logL / 2^(3/2)

    # This would make the angle change 1. Let's fix it up a bit.
    cf *= 1/target_angle_change

    min_en = ustrip(eV, 3/2 * kB * T_e*u"K")

    # This needs to be used with is_cross_section=false.
    # gencf = GENCF_HARDSPHERE(1/ϵ)
    gencf = let cf=cf, min_en=min_en
        GENCF_MANUAL(en -> max(min_en,en)^(-3/2) * cf)
    end

    angle_change = COULOMB_ANGLE_VARYING(C, μ, T_e, ρ_e, cf, min_en)

    cf = CreateCollFreq(params, gas, ptype, name, CFS_ELASTIC(), gencf, angle_dist_cum=angle_change, is_cross_section=false)
end

struct COULOMB_ANGLE_VARYING
    C::Float64
    μ::Float64
    T_e::Float64
    ρ_e::Float64
    cf::Float64
    min_en::Float64
end

function SelectCollCosTheta(angle_change::COULOMB_ANGLE_VARYING, eps)
    @unpack ρ_e,T_e,μ,C,cf,min_en = angle_change

    u = sqrt(2*eps)
    logL = CoulombLogarithm(ρ_e, T_e, u, μ)

    collfreq = max(min_en,eps)^(-3/2) * cf

    temp = 1/collfreq * C * logL / u^3
    # @show 1/collfreq C logL u^3 temp

    temp = ClampMax(temp, 2.0)

    costheta = 1 - temp
end


##########################################
# * Generic stuff here
#----------------------------------------

@inline function WrapThreshold(eps, func, threshold)
    if eps > threshold
        func(eps)
    else
        0.0*uL^3 * ufreq
    end
end
@inline function WrapSSF(eps, func, SSF::Function)
    func(eps) * max(1., SSF(eps))
end
@inline function WrapSSF(eps, func, SSF::INTERP1D)
    func(eps) * max(1., Interp(SSF, eps))
end
    


"""
CreateCollFreq(params, gas, ptype,
               name,
               colltype::CF_STYLE_ENUM, 
               gencf::GENCF_PARAMS
               ;
               threshold=0.0uE,
               ion_sharing_style=ISS_AFE(),
               ion_sharing_ratio=0.5,
               is_cross_section=true)
               
Generate a `COLLFREQ` object, should be used for creating the list to be given
to `gen_cf` of the params object.

The collision type should be chosen from one of those in
[`CF_STYLE_ENUM`](@ref). The details of the form of the collision frequency
should be provided with `gencf`. Options are:

- GENCF_MANUAL - an arbitrary functional form of the cross section
- GENCF_INTERP - interpolate a given list of energy and cross section
- GENCF_INTERPLOG - as with GENCF_INTERP but interpolate on the log values

Some convenience structs:
- GENCF_HARDSPHERE
- GENCF_MAXWELL
- GENCF_ATTACH - for the Ness-Robson attachment benchmark

By default, these are interpreted as a cross section. If a collision frequency
is directly required, then set `is_cross_section` to false.

For inelastic cross sections, `threshold` is important.
For ionisation cross sections, `ion_sharing_style` and `ion_sharing_ratio` is important.

DCS values can be provided using `angle_dist_cum`.
"""               
@xport function CreateCollFreq(params::PARAMS, gas::GAS, ptype::PARTTYPE, name::String, colltype::CF_STYLE_ENUM, gencf::GENCF_PARAMS ; threshold=0.0uE, ion_sharing_style=ISS_AFE(), ion_sharing_ratio=0.5, is_cross_section=true, new_ptype=nothing, angle_dist_cum=nothing, manual_degen=0.)

    manual_degen != 0. && @assert colltype == CFS_INELASTIC_MANUAL_DEGEN()
    
    func = GenCF(gencf)
    m = ptype.mass

    if is_cross_section
        prefac = sqrt(2/m)

        func = let oldfunc=func, prefac=prefac
            eps -> oldfunc(eps) * prefac * sqrt(eps)
        end
    end
    
    if threshold > 0.0uE
        func = let oldfunc=func, threshold=threshold
            eps -> WrapThreshold(eps, oldfunc, threshold)
        end
    end

    if isa(colltype, CFS_ELASTIC)
        if params.static_structure_factor != nothing
            func = let ssf=params.static_structure_factor,oldfunc=func
                eps -> WrapSSF(eps, oldfunc, ssf)
            end
        end
    end

    if isa(colltype, CFS_LOSS)
        if params.weight_reduction != 1.
            func = let oldfunc=func, factor=1/(1-params.weight_reduction)
                eps -> factor * oldfunc(eps)
            end
        end
    end
        
            
    com_func = func

    if gas.tmtr !== nothing
        gas.tmtr == 0uTmtr && error("Should not set temperature to zero, but rather `nothing`.")
        if colltype == CFS_ELASTIC()
            # TODO I should really include a lookup for different regions.
            # Some special cases
            # Disabling these for testing
            # if false #p.static_structure_factor == nothing && name == "Maxwell"
            #     temperature_wrap = SSF_wrap
            # elseif false #p.static_structure_factor == nothing && name == "Hardsphere"
            #     println("Thermal cross section comes from Hardsphere formula.")
            #     hardsphere_a = args[1]

            #     function temperature_wrap(eps)
            #         v = sqrt(2/m)*sqrt(eps)
            #         w = cf.gas.w

            #         R = v/w

            #         #return cf.gas.n0 * hardsphere_a * w * ((R + 1/(2R))*erf(R) + 1/sqrt(pi)*exp(-R^2))
            #         return hardsphere_a * w * ((R + 1/(2R))*erf(R) + 1/sqrt(pi)*exp(-R^2))
            #     end
            # else
            if true
                # Need to do this integral manually.  Could think about
                # implementing the interpolated version in Ristivojevic
                # here. But it seems like it isn't worth the effort.

                nuintlist = @msgwrap "Creating new thermal cross section by integrating manually." IntegrateTemperatureCF(collect(params.eps_grid), wFromT(gas, gas.tmtr), m, func)

                gencf = GENCF_INTERP(collect(params.eps_grid), nuintlist)
                func = GenCF(gencf)
            end
        elseif colltype isa CFS_INELASTIC_ABSTRACT
            if colltype == CFS_INELASTIC()
                kT = kB * gas.tmtr
                boltzfac = 1 / (1 + exp(-threshold / kT))
            else
                boltzfac = manual_degen
            end

            func = let boltzfac=boltzfac, oldfunc=func
                eps -> boltzfac * oldfunc(eps)
            end
        end
    end

    if new_ptype == nothing
        new_ptype_ind = ptype.ind
    elseif isa(new_ptype, PARTTYPE)
        new_ptype_ind = new_ptype.ind
    else
        new_ptype_ind = new_ptype
    end
    new_ptype_ind::Int
        
    COLLFREQ(colltype, name, ptype.ind, gas.ind, threshold, manual_degen, angle_dist_cum, ion_sharing_ratio, ion_sharing_style, new_ptype_ind, func, com_func)
end

function IntegrateTemperatureCF(eps_grid, w, mass, func)
    nuintlist = similar(eps_grid, Q(ufreq/uρ*uVel^2))
    vlist = sqrt.(2/mass * eps_grid)

    intone(x,v,w=w) = x * func(0.5*mass*x^2) * exp(-(x - v)^2 / (w^2))
    inttwo(x,v,w=w) = x * func(0.5*mass*x^2) * exp(-(x + v)^2 / (w^2))

    for i = 1:length(eps_grid)
        v = vlist[i]
        # Quad version - need invokelatest because the cross sections may have just been generated.
        # FIXME: Is this true anymore?
        val1 = Base.invokelatest(quadgk, x -> intone(x,v), max(0uVel, v - 6w), v + 6w)
        val2 = Base.invokelatest(quadgk, x -> inttwo(x,v), 0uVel, v + 6w)
        val = val1[1] - val2[1]

        nuintlist[i] = val
    end

    nuintlist = nuintlist ./ (sqrt(pi) * w * vlist)
end
                          



@xport function PlotAllCollFreqs(params::PARAMS ; in_SI=true, eps=EpsGrid(params), weight=nothing, do_cs=false, plot_kwds=[])
    @eval using Plots

    params = Finalise(params)

    #eps = EpsGrid(params)

    plot_list = []

    for (r_ind,region) in enumerate(params.region_list), (p_ind,ptype) in enumerate(params.ptype_list)
        p = plot(title="Region $(GetName(region)), ptype $(GetName(ptype))" ; plot_kwds...)
        push!(plot_list, p)

        # This is mimicing EvalCollFreq. Perhaps a better way to do this.
        cf_list = COLLFREQ[]
        for gas_ind in 1:length(params.gas_list)
            append!(cf_list, params.cf_matrix[gas_ind, p_ind])
        end
        cf_names = getfield.(cf_list, :name)

        cf_vals = EvalCollFreqs.(eps, Ref(params), r_ind, p_ind) |> matrix

        if in_SI
            cf_vals /= params.time_unit
        end

        if do_cs
            cf_vals ./= (TotDensity(params) / params.len_unit^3)
            # Assume all the same mass
            m = params.ptype_list[1].mass
            cf_vals ./= sqrt(2/m) * sqrt.(eps)' * params.len_unit / params.time_unit
        end


        if weight == nothing
            # nothing!
        elseif weight == :eps
            cf_vals .*= [AverageEpsLoss(cf,params,eps) for cf = cf_list, eps = eps]
        end

        loglog!(eps, permutedims(cf_vals), label=reshape(cf_names, 1, :) ; plot_kwds...)
    end

    display(plot(plot_list...))

    return plot_list
end

AverageEpsLoss(cf::COLLFREQ, args...) = AverageEpsLoss(cf.colltype, cf, args...)

function AverageEpsLoss(::CFS_ELASTIC, cf::COLLFREQ, params::PARAMS, eps::Float64)
    m = params.ptype_list[cf.ptype_ind].mass
    m0 = params.gas_list[cf.gas_ind].m0

    2*m*m0 / (m + m0)^2 * eps
end

function AverageEpsLoss(::CFS_INELASTIC_ABSTRACT, cf::COLLFREQ, params::PARAMS, eps::Float64)
    if eps < cf.threshold
        return 0.
    end
    
    loss_recoil = AverageEpsLoss(CFS_ELASTIC(), cf, params, eps - cf.threshold)

    return cf.threshold + loss_recoil
end

function AverageEpsLoss(::CFS_LOSS, cf::COLLFREQ, params::PARAMS, eps::Float64)
    return 0.
end

function AverageEpsLoss(::CFS_IONISATION, cf::COLLFREQ, params::PARAMS, eps::Float64)
    if eps < cf.threshold
        return 0.
    end
    
    # Ignoring recoil here.
    eps_f = eps - cf.threshold
    if cf.ionisation_sharing_style isa ISS_FRACTION
        eps_share = (1 - cf.ionisation_sharing_ratio) * eps_f
        loss_recoil = AverageEpsLoss(CFS_ELASTIC(), cf, params, eps - cf.threshold - eps_share)
    elseif cf.ionisation_sharing_style isa ISS_AFE
        eps_share = 0.5 * eps_f
        loss_recoil = AverageEpsLoss(CFS_ELASTIC(), cf, params, eps - cf.threshold - eps_share)
    elseif cf.ionisation_sharing_style isa ISS_FUNC
        error("Not implemented")
    end

    return eps_share + cf.threshold + loss_recoil
end

@xport function SSFIntegrateFromSK(en, K, SK, mass=1mₑ)
    sinchion2(cosχ) = sqrt((1-cosχ)/2)
    deltaK(cosχ) = 2*sqrt(2*en * mass)/hbar * sinchion2(cosχ)

    interp = Spline1D(K, SK, k=1)
    integrand(cosχ) = @. 1/2. * (1 - cosχ) * interp(deltaK(cosχ))

    val,err = quadgk(integrand, -1.0, 1.0)
    val
end
