
################################################################################
################################################################################
# * Particle initialisation

@xport struct PIS_UNSET <: PARTICLE_INIT_STYLE end
@xport struct PIS_SEPARABLE <: PARTICLE_INIT_STYLE
    vel_style::PARTICLE_INIT_VEL_STYLE
    spatial_style::PARTICLE_INIT_SPATIAL_STYLE
    time_style::PARTICLE_INIT_TIME_STYLE
    ptype::Union{PARTTYPE,Int}

    PIS_SEPARABLE(vel, spatial=PISS_DELTA(), time=PITS_DELTA(), ptype=1) = new(vel,spatial,time,ptype)
end

@xport struct PIVS_ISOTROPIC <: PARTICLE_INIT_VEL_STYLE
    init_eps::Q(uE)
end
@xport struct PIVS_PUREZ <: PARTICLE_INIT_VEL_STYLE
    init_eps::Q(uE)
end
@xport struct PIVS_GAUSSIAN_ENERGY <: PARTICLE_INIT_VEL_STYLE
    centre::Q(uE)
    width::Q(uE)
end
@xport struct PIVS_MAXWELLIAN <: PARTICLE_INIT_VEL_STYLE
    init_w::Q(uVel)
    init_drift::Q(uVel)
end
function PIVS_MAXWELLIAN(params::PARAMS, mass, temperature, drift=0.0*uVel)
    init_w = wFromT(mass, temperature)
    PIVS_MAXWELLIAN(init_w, drift)
end
@xport struct PIVS_UNIFORM_ENERGY <: PARTICLE_INIT_VEL_STYLE
    min::Q(uE)
    max::Q(uE)
end

@xport struct PISS_DELTA <: PARTICLE_INIT_SPATIAL_STYLE
    pos::PosVec
end
PISS_DELTA() = PISS_DELTA(zero(PosVec))

@xport struct PISS_GAUSSIAN <: PARTICLE_INIT_SPATIAL_STYLE
    spread::Q(uL)
    pos::PosVec
    PISS_GAUSSIAN(spread, pos=zero(PosVec)) = new(spread,pos)
end

@xport struct PITS_DELTA <: PARTICLE_INIT_TIME_STYLE
    time::Q(uT)
    function PITS_DELTA(time)
        @assert(time >= 0uT, "Time for PITS_DELTA can't be negative, it would stuff up sampling in first time bin.")
        new(time)
    end
end
PITS_DELTA() = PITS_DELTA(0.0*uT)


GenerateInitParticle(params::PARAMS) = GenerateInitParticle(params.init_style, params)

function GenerateInitParticle(style::PIS_UNSET, params::PARAMS)
    throw("Init style is unset!")
end

function GenerateInitParticle(style::PIS_SEPARABLE, params::PARAMS)
    if isa(style.ptype, Int)
        ptype = params.ptype_list[style.ptype]
    else
        ptype = style.ptype
    end

    vel = GenerateInitVel(style.vel_style, ptype)
    pos = GenerateInitPos(style.spatial_style, ptype)
    time = GenerateInitTime(style.time_style, ptype)

    # region = FindRegion(pos, params.region_list)

    PARTICLE(pos, vel, time, ptype, params)
end
      
function GenerateInitVel(style::PIVS_ISOTROPIC, ptype::PARTTYPE)
    coschi = rand()*2 - 1
    sinchi = sqrt(1 - coschi^2)
    phi = rand()*2*pi

    velmag = sqrt(2/ptype.mass*style.init_eps)

    vel = VelVec(velmag * [sinchi*cos(phi), sinchi*sin(phi), coschi])
end
function GenerateInitVel(style::PIVS_PUREZ, ptype::PARTTYPE)
    velmag = sqrt(2/ptype.mass*style.init_eps)

    vel = VelVec(0uVel,0uVel,velmag)
end
function GenerateInitVel(style::PIVS_MAXWELLIAN, ptype::PARTTYPE)
    coschi = rand()*2 - 1
    sinchi = sqrt(1 - coschi^2)
    phi = rand()*2*pi

    vel = SampleFromMaxwellian(style.init_w) + VelVec(0.,0.,style.init_drift)

    return vel
end
function GenerateInitVel(style::PIVS_GAUSSIAN_ENERGY, ptype::PARTTYPE)
    coschi = rand()*2 - 1
    sinchi = sqrt(1 - coschi^2)
    phi = rand()*2*pi

    eps = randn() * style.width + style.centre

    velmag = sqrt(2/ptype.mass*eps)
    vel = VelVec(velmag * [sinchi*cos(phi), sinchi*sin(phi), coschi])

    return vel
end
using Distributions
function GenerateInitVel(style::PIVS_UNIFORM_ENERGY, ptype::PARTTYPE)
    coschi = rand()*2 - 1
    sinchi = sqrt(1 - coschi^2)
    phi = rand()*2*pi

    eps = rand(Uniform(style.min,style.max))

    velmag = sqrt(2/ptype.mass*eps)
    vel = velmag * VelVec([sinchi*cos(phi), sinchi*sin(phi), coschi])

    return vel
end

function GenerateInitPos(style::PISS_DELTA, ptype::PARTTYPE)
    style.pos
end
function GenerateInitPos(style::PISS_GAUSSIAN, ptype::PARTTYPE)
    z = randn()*style.spread
    style.pos + PosVec(0uL,0uL,z)
end

function GenerateInitTime(style::PITS_DELTA, ptype::PARTTYPE)
    style.time
end
