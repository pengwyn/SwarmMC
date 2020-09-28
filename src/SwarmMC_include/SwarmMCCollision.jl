
@xport function SampleFromMaxwellian(w)
    u = @SVector(randn(3))*w / sqrt(2)
end

function SelectNeutralU(cf::COLLFREQ, gas::GAS, part::PARTICLE, w=0uVel)
    if w == 0uVel || cf.colltype != CFS_ELASTIC()
        return XYZ(0,0,0)*uVel
    end
    
    m = part.cur_mass
    
    vmag = norm(part.vel)

    max_w_mult = 5

    high_speed = vmag + w*max_w_mult
    min_cf_eps = 1e-5eV
    lowest_speed = sqrt(2*min_cf_eps / m)
    low_speed = max(vmag - w*max_w_mult, lowest_speed)
    # low_speed = max(vmag - w*max_w_mult, zero(w))

    high_eps = m/2 * high_speed^2
    low_eps = m/2 * low_speed^2

    innercffunc = cf.com_func

    # high_cf = innercffunc(high_eps)::Q(ufreq/uρ)
    # low_cf = innercffunc(low_eps)::Q(ufreq/uρ)
    high_cf = innercffunc(high_eps)
    low_cf = innercffunc(low_eps)

    max_cf = max(high_cf, low_cf)

    iter = 0
    local u
    while true
        iter += 1
        if iter >= 100000
            eps = EpsFromVel(part)
            @show u w high_cf low_cf max_cf eps part.vel
            throw("Rejected a lot of choices for neutral vel. Something must be wrong!")
        end

        u = SampleFromMaxwellian(w)

        relspeed = norm(u-part.vel)
        thiseps = m/2 * relspeed^2

        # thiscf = innercffunc(thiseps)::Q(ufreq/uρ)
        thiscf = innercffunc(thiseps)
        dist_prob = thiscf / max_cf

        # TODO: Measurements of choice success rates.

        if rand() < dist_prob
            break
        end
    end

    # TODO: Measurements of neutral choice
    u
end

@xport struct DELTA_ANGLE_DIST
    eps_list::Vector{Q(uE)}
    x_locs::Vector{Float64}
    Rvals::Matrix{Float64}
end

SelectCollCosTheta(angle_dist_cum::Nothing, eps) = rand()*2 - 1
SelectCollCosTheta(angle_dist_cum::Function, eps) = angle_dist_cum(eps, rand())
function SelectCollCosTheta(angle_dist_cum::INTERP2D, eps)
    BisectInterpInternal.InterpRevY(angle_dist_cum, eps, rand())
end
SelectCollCosTheta(angle_dist_cum::DELTA_ANGLE_DIST, eps) = AngleFromDeltaRvals(eps, rand(), angle_dist_cum.eps_list, angle_dist_cum.x_locs, angle_dist_cum.Rvals)
    
function Collision(params::PARAMS, cf::COLLFREQ, part::PARTICLE, energy_subtraction=0.0)
    m0 = params.gas_list[cf.gas_ind].m0
    m = part.cur_mass

    #mu = m*m0 / (m0 + m)
    # This version also works if m0 is inf
    mu = m / (1. + m/m0)
    mu_m0 = mu / m0
    mu_m = mu / m

    gas = params.gas_list[cf.gas_ind]
    w = wFromT(gas, part.pos)
    if params.handle_temperature_bias
        u = SelectNeutralU(cf, params.gas_list[cf.gas_ind], part, w)
    else
        u = SampleFromMaxwellian(w)
    end


    comvel = part.vel*mu_m0 + u*mu_m
    relvel = part.vel - u

    relvel_mag = norm(relvel)
    if relvel_mag == 0
        orig_costheta = 1.
        orig_sintheta = 0.
        orig_cosphi = 1.
        orig_sinphi = 0.
    else
        orig_costheta = relvel[3] / relvel_mag
        orig_sintheta = sqrt(1 - orig_costheta^2)
        if orig_sintheta == 0
            orig_cosphi = 1.
            orig_sinphi = 0.
        else
            orig_cosphi = relvel[1] / relvel_mag / orig_sintheta
            orig_sinphi = relvel[2] / relvel_mag / orig_sintheta
        end
    end

    if energy_subtraction != 0
        # Need to make sure it's larger than the energy available.
        K = EpsFromVel(relvel, mu)

        if energy_subtraction <= K
            factor = sqrt(1 - energy_subtraction / K)

            relvel_mag *= factor
        end
    end


    eps = EpsFromVel(part)
    coll_costheta = SelectCollCosTheta(cf.angle_dist_cum, eps)::Float64
    coll_sintheta = sqrt(1 - coll_costheta^2)
    coll_phi = rand()*2*pi

    thex = coll_sintheta * cos(coll_phi)
    they = coll_sintheta * sin(coll_phi)
    thez = coll_costheta
    coll_vel = relvel_mag * XYZ(thex, they, thez)

    Ry = @SMatrix [orig_costheta 0. orig_sintheta ; 0. 1. 0. ; -orig_sintheta 0. orig_costheta]
    Rz = @SMatrix [orig_cosphi -orig_sinphi 0. ; orig_sinphi orig_cosphi 0. ; 0. 0. 1.]

    new_relvel = Rz * Ry * coll_vel

    part.vel = mu_m*new_relvel + comvel
    UpdateEpsBin(part, params.eps_bin_grid)

    part
end


function ReduceEnergy(part::PARTICLE, deps::Q(uE), eps_bin_grid::Vector)
    eps = EpsFromVel(part)
    final_eps = eps - deps
    velscale = sqrt(final_eps / eps)
    part.vel *= velscale
    UpdateEpsBin(part, eps_bin_grid)
end

# function MeasureNonCons(props_set::PROPS_SET, params::PARAMS, part::PARTICLE, weight::Float64)
#     prefac = exp2(part.log2fac - props_set.log2fac) * weight * part.weight
#     #prefac = exp2(part.log2fac) * weight * part.weight
#     props_set.time[part.ptype.ind].noncons        += prefac
#     props_set.time[part.ptype.ind].noncons_pos    += prefac * part.pos
#     props_set.time[part.ptype.ind].noncons_possqr += prefac * part.pos.^2

#     if params.measure_z_quants != :no
#         props_set.z[part.z_bin, part.ptype.ind].noncons        += prefac
#         props_set.z[part.z_bin, part.ptype.ind].noncons_pos    += prefac * part.pos
#         props_set.z[part.z_bin, part.ptype.ind].noncons_possqr += prefac * part.pos.^2
#     end
#     if params.measure_r_quants != :no
#         props_set.r[part.r_bin, part.ptype.ind].noncons        += prefac
#         props_set.r[part.r_bin, part.ptype.ind].noncons_pos    += prefac * part.pos
#         props_set.r[part.r_bin, part.ptype.ind].noncons_possqr += prefac * part.pos.^2
#     end
# end

############################################
# * Particular collision types

function HandleCollFreq(cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list::Vector{PARTICLE})
    HandleCollFreq(cf.colltype, cf, p, part, extra_part_list)
    nothing
end


function HandleCollFreq(::CFS_ELASTIC, cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list)
    SSF = p.static_structure_factor
    eps = EpsFromVel(part)
    
    if SSF === nothing
        Collision(p, cf, part)
    else
        local SSF_val
        if isa(SSF, Function)
            SSF_val = SSF(eps)
        elseif isa(SSF, INTERP1D)
            SSF_val = Interp(SSF, eps)
        end

        R = rand()*max(SSF_val, 1.)

        debug_before = eps
        
        if R > 1.
            before_eps = eps
            before_vel = part.vel
            Collision(p, cf, part)
            eps = EpsFromVel(part)
            part.vel *= sqrt(before_eps / eps)
            UpdateEpsBin(part, p.eps_bin_grid)

        elseif R < SSF_val
            Collision(p, cf, part)
        else
            before_eps = eps
            before_vel = part.vel
            Collision(p, cf, part)
            eps = EpsFromVel(part)
            part.vel = before_vel * sqrt(eps / before_eps)
        end
    end
    nothing
end

function HandleCollFreq(::CFS_LOSS, cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list)
    #MeasureNonCons(props_set, p, part, -p.weight_reduction)
    part.weight *= (1-p.weight_reduction)
    nothing
end

# function HandleCollFreq(::CFS_NULL, cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list)
#     meas.num += part.weight
#     return part
# end

# function HandleCollFreq(::CFS_FAKE_NONCONS, cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list)
#     if p.fake_noncons_rate > 0
#         if p.fake_noncons_split_weight
#             part.weight /= 2.
#         end
#         push!(extra_part_list, copy(part))
#     else
#         part.weight = 0.
#     end

#     return part
# end

function HandleCollFreq(::CFS_INELASTIC_ABSTRACT, cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list)
    Collision(p, cf, part, cf.threshold)
    nothing
end

function HandleCollFreq(::CFS_SUPERELASTIC, cf::COLLFREQ, p::PARAMS, part::PARTICLE, extra_part_list)
    Collision(p, cf, part, cf.threshold)
    nothing
end

function HandleCollFreq(::CFS_IONISATION, cf::COLLFREQ, params::PARAMS, part::PARTICLE, extra_part_list)
    ReduceEnergy(part, cf.threshold, params.eps_bin_grid)

    other_eps = ChooseIonisationEnergy(part, cf)

    ReduceEnergy(part, other_eps, params.eps_bin_grid)

    Collision(params, cf, part)

    if cf.new_ptype_ind != -1
        if IsTrue(params.ionisation_as_weight) && cf.new_ptype_ind == part.ptype_ind
            if rand(Bool)
                # Update the particle energy - using dodgy "reduce" but it may
                # actually increase it.
                # other_eps = EpsFromVel(part) - other_eps
                ReduceEnergy(part, -(other_eps - EpsFromVel(part)), params.eps_bin_grid)
            end

            # Double the weight
            # Can't change log2fac in the ode version
            part.weight *= 2
        else
            # TODO: Should still allow for directional choice of ionisation.
            other_part = GenerateInitParticle(PIS_SEPARABLE(PIVS_ISOTROPIC(other_eps), PISS_DELTA(part.pos), PITS_DELTA(part.time), params.ptype_list[cf.new_ptype_ind]), params)
            other_part.weight = part.weight
            other_part.log2fac = part.log2fac
            # I don't know why this has to be done this way... but I think there's a good reason for it.
            other_part.t_bin = part.t_bin

            # TODO: Perhaps need this to be other_part??
            #MeasureNonCons(props_set, params, part, 1.)
            push!(extra_part_list, other_part)
        end
    end

    return part
end

function ChooseIonisationEnergy(part::PARTICLE, cf::COLLFREQ)
    local other_eps
    eps = EpsFromVel(part)
    if isa(cf.ionisation_sharing_style, ISS_FRACTION)
        other_eps = eps * (1-cf.ionisation_sharing_ratio)
    elseif isa(cf.ionisation_sharing_style, ISS_AFE)
        other_eps = eps * rand()
    elseif isa(cf.ionisation_sharing_style, ISS_FUNC)
        other_eps = cf.ionisation_sharing_style.func(eps)
    else
        # Lookups need to go here.
        error("Not implemented")
    end

    # For safety, max the other_eps to be at most part.eps
    if other_eps > eps
        @warn "The ionisation eps was larger than the input eps" eps other_eps
        other_eps = eps
    end

    other_eps
end

function HandleCollFreq(::CFS_CHANGETYPE, cf::COLLFREQ, params::PARAMS, part::PARTICLE, extra_part_list)
    ChangePType(part, params.ptype_list[cf.new_ptype_ind], params)

    return part
end
