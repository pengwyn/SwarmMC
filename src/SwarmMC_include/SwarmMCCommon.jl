# Introduce a special unit for number of particles, which can also be used as
# the particle weight.
# New version where I just use the existing dimension and define it in terms of moles.
# Should only be an issue when using upreferred
@unit uN "ñ" Particle Unitful.mol/6.02214076e23 false
const localunits = Unitful.basefactors

const uL = u"Å"
const uE = u"eV"
const uM = u"mₑ"
const uT = sqrt(1.0 * uM * uL^2 / uE)
const uFonN = u"Td"
const uF = u"V/m"
const uV = u"V"
const uC = u"q"
const uB = u"T"
const uTmtr = u"K"

const uVel = uL / uT
const uAcc = uL / uT^2
const uρ = uL^-3
const uσ = uL^2
const ufreq = uT^-1

using Unitful: AbstractQuantity, Unitlike, Units

Base.convert(::Type{Union{Union,Function,T}}, x::T2) where {T<:Union{AbstractQuantity,Units}, T2<:AbstractQuantity} = convert(T, x)

const UnitQuants = Union{AbstractQuantity,Number}
Q(unit) = typeof(1.0 * unit)
Q(unit::AbstractQuantity) = typeof(unit)
Q(unit::Type{T}) where {T<:UnitQuants} = T
QComb(unit1::Type{<:UnitQuants}, unit2::Type{<:UnitQuants}) = typeof(zero(Q(unit1)) * zero(Q(unit2)))
QComb(unit1, unit2) = QComb(Q(unit1), Q(unit2))

const XYZ{T} = SVector{3,T}
const PosVec = XYZ{Q(uL)}
const VelVec = XYZ{Q(uVel)}

Base.string(x::Type{T}) where {EL,T <: XYZ{EL}} = isconcretetype(x) ? "XYZ{$EL}" : "XYZ"
Base.show_datatype(io::IO,x::Type{<:XYZ}) = print(io, string(x))

# This is a convenient conversion for constructors when I have a possible constant/function form
Base.convert(::Type{Union{Nothing,Function,T}}, x::T2) where {T<:XYZ, T2<:AbstractArray} = convert(T, x)
Base.convert(::Type{Union{Function,T}}, x::T2) where {T<:XYZ, T2<:AbstractArray} = convert(T, x)
Base.convert(::Type{Union{Nothing,Function,T}}, x::T2) where {T<:Number, T2<:Number} = convert(T, x)
Base.convert(::Type{Union{Function,T}}, x::T2) where {T<:Number, T2<:Number} = convert(T, x)

################################################################################
# ** Measurements

"""
MEAS_BIN
    
Define the times, positions or energies in which measurements are to be divided
into. It also selects whether measurements should be cumulative or
instantaneous.

The constructor takes symbols, e.g. `MEAS_BIN(:cum,:t)` represents cumulative
measurements in time bins and `MEAS_BIN(:ss,:eps,:r)` are instantaneous
measurements in energy and radial bins, only for the steady-state times, defined
by `steady_state_frac` in `PARAMS`.
"""
@xport @AutoParm struct MEAS_BIN
    steady_state::AUTO <: TypeBool
    cumulative::AUTO <: TypeBool
    # user::AUTO <: TypeBool

    with_t::AUTO <: TypeBool
    with_r::AUTO <: TypeBool
    with_z::AUTO <: TypeBool
    with_eps::AUTO <: TypeBool
    with_costh::AUTO <: TypeBool
end

# Convenience constructors allowing for e.g. MEAS_BIN(:t, :cum, :ss)
MEAS_BIN(syms::NTuple{N,Symbol} where N) = MEAS_BIN(syms...)
function MEAS_BIN(syms::Symbol...)
    steady_state = cumulative = with_cf = with_t = with_r = with_z = with_eps = with_costh = TypeFalse()

    for sym in syms
        if sym == :ss
            steady_state = TypeTrue()
            cumulative = TypeTrue()
        elseif sym == :cum
            cumulative = TypeTrue()
        elseif sym == :t
            with_t = TypeTrue()
        elseif sym == :r
            with_r = TypeTrue()
        elseif sym == :z
            with_z = TypeTrue()
        elseif sym == :eps
            with_eps = TypeTrue()
        elseif sym == :costh
            with_costh = TypeTrue()
        else
            error("Unknown meas type requested with symbol $sym.")
        end
    end

    if Bool(with_costh)
        @assert !Bool(cumulative) && !Bool(steady_state)
    end

    MEAS_BIN(steady_state, cumulative, with_t, with_r, with_z, with_eps, with_costh)
end

function Base.show(io::IO, x::MEAS_BIN)
    print(io, "MEAS_BIN(" * GetNames(x) * ")")
end

# GetSymbols(meas_type::MEAS_BIN) = GetSymbols(typeof(meas_type))
# WARNING! This relies on the types being in an expected order
function GetSymbols(x::MEAS_BIN)
    thelist = []
    
    IsTrue(x.steady_state) && push!(thelist, :ss)
    IsTrue(x.cumulative) && push!(thelist, :cum)

    for sym in [:t, :r, :z, :eps, :costh]
        if IsTrue(getproperty(x, Symbol(:with_, sym)))
            push!(thelist, sym)
        end
    end

    thelist
end

function GetNames(meas_type::Union{MEAS_BIN, Type{T} where T<: MEAS_BIN})
    thelist = GetSymbols(meas_type)

    join(string.(thelist), ",")
end

##############################
# * COMB_MEAS
#----------------------------

"""
MEAS_QUANT

References the quantity that is measured and will be placed into `MEAS_BIN`.
A quantity defined this way should then extend `InstMeas` to implement the
measure itself.
"""
struct MEAS_QUANT{label, unit, kind} end
MEAS_QUANT(label, unit, kind=Float64) = MEAS_QUANT{label,unit,kind}()

"""
InstValue(::Val{sym}, weight, log2fac, pos, vel, eps)

The value of the quantity (`sym`) being measured, given the current particle
properties are `pos`, `vel`, `eps` and the particle's weight is given by
`weight + 2^log2fac`.
"""
InstValue(::Val{T}, weight, log2fac, pos, vel, eps) where T = error("InstValue not defined for quantity $T.")

"""
Log2Fac(::Val{sym}, log2fac)

What log2fac should be assigned to this quantity. Only relevant for measurements
that are considering the weights.
"""
Log2Fac(::Val, log2fac) = log2fac

FullType(::Type{Float64}, unit) = QComb(unit,uN)
FullType(::Type{XYZ}, unit) = XYZ{QComb(unit,uN)}
FullType(::Type{Vector}, unit) = Vector{QComb(unit,uN)}

FullTypeCumul(kind, unit) = FullType(kind, QComb(unit,uT))

Q_energy = MEAS_QUANT(:energy, Q(uE))
Q_energysqr = MEAS_QUANT(:energy_pow2, Q(uE^2))
Q_vel = MEAS_QUANT(:vel, Q(uVel), XYZ)
Q_velsqr = MEAS_QUANT(:vel_pow2, Q(uVel^2), XYZ)
Q_pos = MEAS_QUANT(:pos, Q(uL), XYZ)
Q_possqr = MEAS_QUANT(:pos_pow2, Q(uL^2), XYZ)
Q_posvel = MEAS_QUANT(:pos_vel, Q(uL*uVel), XYZ)
Q_epsflux = MEAS_QUANT(:energy_vel, Q(uE*uVel), XYZ)
Q_epspressure = MEAS_QUANT(:energy_vel_vel, Q(uE*uVel^2), XYZ)

Q_denom = MEAS_QUANT(:denom, Q(1.))
Q_sqrdenom = MEAS_QUANT(:sqrdenom, Q(uN))
Q_invdenom = MEAS_QUANT(:invdenom, Q(uN^-2))
Q_weight = MEAS_QUANT(:weight, Q(uN))
Q_preweight = MEAS_QUANT(:preweight, Q(uN))
Q_log2weight = MEAS_QUANT(:log2weight, Q(1.))

Q_basic_list = (Q_denom, Q_invdenom, Q_sqrdenom,
                Q_weight, Q_preweight, Q_log2weight,
                Q_energy, Q_vel, Q_velsqr,
                Q_pos, Q_possqr,
                Q_posvel,
                Q_epsflux, Q_epspressure)

Q_quick_list = (Q_denom, Q_sqrdenom,
                Q_weight,
                Q_energy, Q_vel, Q_pos,
                Q_possqr, Q_posvel)

@AutoParm struct MEAS_LOG2{LABEL, QUANT, BINS, N}
    label::Val{LABEL}

    mbin::BINS <: MEAS_BIN

    ptype_ind::Int
    cf_ind::Int
    
    vals::Array{QUANT,N}
    log2w::Array{Int,N}
end

################################################################################
# * PROPS_OUT

"""
PROPS_OUT

Structure to store measurements taken during a simulation. This struct supports
convenient querying of different quantities for different bins.
"""
@xport @AutoParm mutable struct PROPS_OUT
    meas::Array{MEAS_LOG2,3}
    walltime::Q(u"s")
    num_particles::Q(uN)
end

function Base.show(io::IO, props_out::PROPS_OUT)
    bins = unique(getproperty.(props_out.meas, :mbin))
    meastypes = join("(".*GetNames.(bins).*")", ",")
    print(io, "PROPS_OUT with N=$(props_out.num_particles), meas=$meastypes, walltime=$(round(props_out.walltime ; sigdigits=3))")
end


######################################################
# * Optional function/values
#----------------------------------------------------

GetVal(val::Number, args...) = val
GetVal(val::XYZ, args...) = val
GetVal(func::Function, args...) = func(args...)
GetVal(val::Nothing, args...) = 0.0 # Not so useful - no units! Should at least cause an error.

################################################################################
# ** Gas
const last_gas_ind = Ref(0)
NextGasInd() = (last_gas_ind[] += 1)

"""
GAS(; name::Symbol=:gas, m0, ρ, tmtr)

The definition of a gas. The temperature `tmtr` defaults to zero, which itself
must be specified as `nothing`.

Both `ρ` and `tmtr` can accept a number or a function which itself takes one
argument: `pos`.
"""
@AutoParm @xport struct GAS
    name::Symbol = :gas
	m0::Q(uM)
    ρ::AUTO <: Union{Function,Q(uρ)}
    tmtr::AUTO <: Union{Nothing,Function,Q(uTmtr)} = nothing
    ind::Int = NextGasInd()
end


################################################################################
# ** PARTTYPE

# TODO: Need to make these type params not slow somehow. I suspect this means
# not including them as parameters, because that at least makes the outer type a
# concrete type.
# TODO: However the same thing applies to NAME. I think I need to
# ditch this in favour of faster calculation (but less flexibility in the
# definition).

const last_ptype_ind = Ref(0)
NextPTypeInd() = (last_ptype_ind[] += 1)
# @xport mutable struct PARTTYPE{NAME,
#                         T_DELAY_LOOKUP,
#                         T_DELAY_EVENT}

"""
PARTTYPE(;name=:particle, mass=1mₑ, charge=1q)

The definition of the type of particle, giving a mass and charge.
"""
@AutoParm @xport mutable struct PARTTYPE
    name::Symbol = :particle
    mass::Q(uM) = 1.0u"mₑ"
    charge::Q(uC) = 1.0u"q"

    # delay_lookup::T_DELAY_LOOKUP
    # delay_event::T_DELAY_EVENT
    # delay_lookup = ()->Inf*uT
    # delay_event = DELAY_EVENT_NULL()

    ind::Int = NextPTypeInd()
end

################################################################################
# ** COLLFREQ

@TypeEnumExport CF_STYLE_ENUM CFS_ELASTIC CFS_LOSS CFS_IONISATION CFS_NULL CFS_FAKE_NONCONS CFS_SUPERELASTIC CFS_CHANGETYPE CFS_UNSET
export CFS_INELASTIC,
    CFS_INELASTIC_MANUAL_DEGEN
abstract type CFS_INELASTIC_ABSTRACT <: CF_STYLE_ENUM end
struct CFS_INELASTIC <: CFS_INELASTIC_ABSTRACT end
struct CFS_INELASTIC_MANUAL_DEGEN <: CFS_INELASTIC_ABSTRACT end

@doc """
CF_STYLE_ENUM

The type of collision frequency process.

These types are hopefully obvious:
- CFS_ELASTIC
- CFS_INELASTIC
- CFS_LOSS
- CFS_IONISATION

The types that might need some more explanation are:
- CFS_NULL - for introducing manually configured null collisions.
- CFS_INELASTIC_MANUAL_DEGEN - an inelastic cross section with manually set ground/excitation populations, as opposed to using the gas temperature.
- CFS_CHANGETYPE - changes the type of particle.

""" CF_STYLE_ENUM

@TypeEnumExport ISS_ENUM ISS_FRACTION ISS_AFE #ISS_FUNC
@xport struct ISS_FUNC{T} <: ISS_ENUM
    func::T
end

@doc """
ISS_ENUM

Options for ionisation energy sharing

- ISS_FRACTION - specify an exact constant fraction
- ISS_AFE - all fractions equiprobable
- ISS_FUNC - the fraction is chosen at the time of collision.
""" ISS_ENUM

@xport mutable struct COLLFREQ
    colltype::CF_STYLE_ENUM
    name::String

    ptype_ind::Int
    gas_ind::Int

    threshold::Q(uE)

    manual_degen::Float64

    angle_dist_cum

    ionisation_sharing_ratio::Float64
    ionisation_sharing_style::ISS_ENUM

    new_ptype_ind::Int

    func
    com_func
end

# @xport abstract type DELAY_EVENT end
# @xport struct DELAY_EVENT_NULL <: DELAY_EVENT end
# @xport struct DELAY_EVENT_CHANGETYPE <: DELAY_EVENT
#     new_ptype::PARTTYPE
#     return_vel
#     DELAY_EVENT_CHANGETYPE(new_ptype, return_vel=0.) = new(new_ptype, return_vel)
# end
# @xport struct DELAY_EVENT_LOSS <: DELAY_EVENT end


################################################################################
# ** PARTICLE

@xport @AutoParm mutable struct PARTICLE
    pos::PosVec = zero(PosVec)
    vel::VelVec = zero(VelVec)
    time::Q(uT) = zero(uT)
    coll_C::Float64 = GenerateCollC()

    ptype_ind::Int = 1
    # region_ind::Int

    weight::Q(uN) = 1uN
    log2fac::Int = 0

    eps_bin::Int = 0
    z_bin::Int = 0
    r_bin::Int = 0
    t_bin::Int = 0

    cur_mass::Q(uM) = zero(Q(uM))
    cur_charge::Q(uC) = zero(Q(uC))
end

Base.copy(x::PARTICLE) = GenericCopy(x)



################################################################################
# ** PARAMS

abstract type PARTICLE_INIT_STYLE end
abstract type PARTICLE_INIT_VEL_STYLE end
abstract type PARTICLE_INIT_SPATIAL_STYLE end
abstract type PARTICLE_INIT_TIME_STYLE end


@TypeEnumExport GNS_ENUM GNS_DOUBLE GNS_REGEN_ALL GNS_NOTHING GNS_UPDATE_LOG2FAC

@doc """
GNS_ENUM

The style for generating or removing particles for the next time bin.

- GNS_NOTHING - don't ever generate or remove new particles
- GNS_UPDATE_LOG2FAC - as with GNS_NOTHING but update the log2fac of the particles.
- GNS_DOUBLE - when there are too few particles, double those that are left.
- GNS_REGEN_ALL - when the effective number of particles gets too high or low,
                  regenerate by sampling the same weights of particles.
""" GNS_ENUM

@enum StepEvents Event_none Event_coll Event_epsbin Event_zbin Event_region Event_rbin Event_delay
Base.to_index(x::StepEvents) = Int(x)
Base.length(x::Type{StepEvents}) = length(instances(StepEvents))


"""
PARAMS(; kwds...)

The main parameters object. The parameters can be set as keywords or by mutating
the structure afterwards. When finished setting the parameters, call
`Finalise(...)` on the object which will set the structure type parameters and
fill in the necessary interal-use-only fields.

Keywords that must be set (no defaults):
- `init_style` - how particles are generated
- `t_grid` - the time grid used for the simulation and for measurements
- `z_grid`, `r_grid` - the grids used for measurements
- `gas_list` - a list of `GAS` objects
- `ptype_list` - a list of `PARTTYPE` objects
- `gen_cf` - a function that returns the list of collision frequencies for each
  gas and particle type combination. Function signature `(params, gas, ptype) -> ...`

Keywords that have defaults:
- `eps_grid` - the energy grid used for measurements
- `meas_bins` - the bins for measurements
- `meas_quants` - the quantities measured
- `steady_state_timefrac` - defines the time at which steady-state measurements
  are taken from.

- `E_field` - the electric field. Can be `nothing`, an XYZ vector, or a function
    of signature `(pos,time) -> ...`.
- `B_field` - the magnetic field. Can be `nothing`, an XYZ vector, or a function
    of signature `(pos,time) -> ...`.

- `min_log2_weight` - the point at which if a particle's weight is less than `2^min_log2_weight`
- `weight_reduction` - the percentage of the weight that is lost in a loss collision.
- `ionisation_as_weight` - whether to either create a secondary particle in
  ionisation, or to double the weight. 

- `generate_next_step_style` - one of [`GNS_ENUM`](@ref) types for refreshing the
    particle list for each timestep.
- `only_regen_for_too_many` - when true, don't ever generate extra particles
    when using `GNS_REGEN_ALL`.
- `max_substepsize` - the maximum timestep. Use this when the time grid is too
    coarse to handle changing numbers of particles.
- `runaway_threshold` - the multiplicative factor to judge when a runaway of
  particles is occuring, which will stop the simulation.

- `static_structure_factor` - the static structure factor used for elastic collisions

- `save_name` - the name for saving, see [`MakeSaveName`](@ref)

- `do_event_meas` - if true, store information in `event_counts`
- `show_events` - log whenever an event occurs.

- `handle_temperature_bias` - set this to false to force sampling gas velocities
  directly (and incorrectly) from a Maxwellian distribtion.

- `user_callbacks` - a list of additional callbacks to provide to the ODE
  solver. These callbacks will directly interface with
  `DifferentialEquations.jl`. For reference, the `p` in the integrator object
  will be of type [`INTEGRATOR`](@ref). 
"""
@xport @AutoParm mutable struct PARAMS
    len_unit::Q(uL) = 1.0uL
    eps_unit::Q(uE) = 1.0uE
    mass_unit::Q(uM) = 1.0uM
    time_unit::Q(uT) = 1.0uT

    init_style::AUTO = PIS_UNSET()

    t_grid::AUTO <: AbstractVector{Q(uT)} = Q(uT)[]
    eps_grid::AUTO <: AbstractVector{Q(uE)} = LogRange(1e-5, 100.0, 1001, as_val=true, prefactor=1eV)
    z_grid::AUTO <: AbstractVector{Q(uL)} = Q(uL)[]
    r_grid::AUTO <: AbstractVector{Q(uL)} = Q(uL)[]

    # Internal use only
    eps_bin_grid::AUTO <: AbstractVector{Q(uE)} = Q(uE)[]
    z_bin_grid::AUTO <: AbstractVector{Q(uL)} = Q(uL)[]
    r_bin_grid::AUTO <: AbstractVector{Q(uL)} = Q(uL)[]
    costh_bin_grid::AUTO <: AbstractVector{Float64} = range(-1.0,1.0, length=51)

    meas_bins::AUTO <: Tuple{Vararg{MEAS_BIN}} = (MEAS_BIN(:t,:cum),)
    meas_quants::AUTO <: Tuple{Vararg{MEAS_QUANT}} = Q_basic_list

    # Internal use only
    has_cum_meas::AUTO = nothing

    # meas_num_l::Int = 0

    steady_state_timefrac::Float64 = 0.5

	gas_list::AUTO = GAS[]
    ptype_list::AUTO = PARTTYPE[]

    gen_cf::Function = (params, gas, ptype) -> []

    # Internal use only
    cf_matrix::AUTO = nothing
    cf_list::AUTO = nothing

    E_field::AUTO <: Union{Nothing,Function,XYZ{typeof(1.0uF)}} = nothing
    B_field::AUTO <: Union{Nothing,Function,XYZ{typeof(1.0uB)}} = nothing

    min_log2_weight::Float64 = -200.0
    weight_reduction::Float64 = 0.5
    ionisation_as_weight::AUTO <: TypeBool = TypeTrue()

    generate_next_step_style::AUTO = GNS_REGEN_ALL()
    only_regen_for_too_many::Bool = false
    max_substepsize::Q(uT) = Inf*uT
    runaway_threshold::Int = 100

    static_structure_factor::AUTO = nothing

    save_name::String = "DefaultSave"

    do_event_meas::AUTO = TypeFalse()
    event_counts::Vector{Base.Enums.basetype(StepEvents)} = zeros(Int, length(StepEvents))
    show_events::AUTO <: TypeBool = TypeFalse()

    handle_temperature_bias::Bool = true

    user_callbacks::AUTO = []

    finalised::Bool = false
end

# For convenience
Base.broadcastable(x::PARAMS) = Ref(x)

# Annoying with all the finalisation in here, but very convenient otherwise.
"""
TGrid(params)

Return the midpoint of all time bins, which is the value relevant for cumulative measurements.
"""
@xport TGrid(params::PARAMS) = RunningAvg(Finalise(params).t_grid)
"""
TInstGrid(params)

Return the points at which instantaneous measurements correspond to.
"""
@xport TInstGrid(params::PARAMS) = Finalise(params).t_grid[begin+1:end]
"""
ΔT(params)

Return the size of the cumulative time measurement bins.
"""
@xport ΔT(params::PARAMS) = diff(Finalise(params).t_grid)
"""
RGrid(params)

Return the midpoint of the radial measurement bins.
"""
@xport RGrid(params::PARAMS) = RunningAvg(Finalise(params).r_bin_grid)
"""
ΔR(params)

Return the size of the radial measurement bins.
"""
@xport ΔR(params::PARAMS) = diff(Finalise(params).r_bin_grid)
"""
ZGrid(params)

Return the midpoint of the spatial z measurement bins.
"""
@xport ZGrid(params::PARAMS) = RunningAvg(Finalise(params).z_bin_grid)
"""
ΔZ(params)

Return the size of the spatial z measurement bins.
"""
@xport ΔZ(params::PARAMS) = diff(Finalise(params).z_bin_grid)
"""
EpsGrid(params)

Return the midpoint of the energy measurement bins.
"""
@xport EpsGrid(params::PARAMS) = RunningAvg(Finalise(params).eps_bin_grid)
"""
ΔEps(params)

Return the size of the energy measurement bins.
"""
@xport ΔEps(params::PARAMS) = diff(Finalise(params).eps_bin_grid)
"""
CosthGrid(params)

Return the midpoint of the costheta measurement bins.
"""
@xport CosthGrid(params::PARAMS) = RunningAvg(Finalise(params).costh_bin_grid)
"""
ΔCosth(params)

Return the size of the costheta measurement bins.
"""
@xport ΔCosth(params::PARAMS) = diff(Finalise(params).costh_bin_grid)

function Base.show(io::IO, params::PARAMS)
    print(io, "PARAMS with Gas=[" * join(getproperty.(params.gas_list, :name), ",") * "], "
                      * "Ptype=[" * join(getproperty.(params.ptype_list, :name), ",") * "]")
                     #* "Region=[" * join(GetName.(params.region_list), ",") * "]")
end
################################################################################
################################################################################
################################################################################
# ** Constructors

function MEAS_LOG2(params::PARAMS, quant, mbin::MEAS_BIN, ptype_ind=1, cf_ind=0)
    label,unit,kind = typeof(quant).parameters

    if kind <: MVector || kind == Vector
        error("Not implemented")
    end
    if IsTrue(mbin.cumulative)
        # Need to special case the sqrdenom variant
        # This is silly... I'm not going to bother since the results are useless
        # anyway (a single "measurement" is ill defined in the cumulative case).
        # if label == :sqrdenom
        #     quant = FullTypeCumul(Float64, Q(uN*uT))
        # elseif label == :invdenom
        #     quant = FullTypeCumul(Float64, Q(uN^-2*uT^-2))
        # else
            quant = FullTypeCumul(kind, unit)
        # end
    else
        quant = FullType(kind, unit)
    end
    dims = Dims(params, mbin)
    if dims == ()
        dims = 1
    end

    vals = zeros(quant, dims...)
    log2w = fill(typemin(Int), dims...)

    MEAS_LOG2(Val(label), mbin, ptype_ind, cf_ind, vals, log2w)
end
    

GenerateCollC() = -log(rand())

function PARTICLE(pos, vel, time, ptype, params)
    eps = EpsFromVel(vel, ptype.mass)
    eps_bin = EpsBin(eps, params.eps_bin_grid)

    z_bin = ZBin(pos, params.z_bin_grid)
    r_bin = RBin(pos, params.r_bin_grid)
    t_bin = 0

    weight = 1uN
    log2fac = 0

    part = PARTICLE(; pos, vel, time, ptype_ind=ptype.ind, eps_bin, z_bin, r_bin)
    UpdateParticleLookups(part, params)

    return part
end

function UpdateParticleLookups(part::PARTICLE, params::PARAMS)
    # region = params.region_list[part.region_ind]
    ptype = params.ptype_list[part.ptype_ind]

    part.cur_mass = ptype.mass
    part.cur_charge = ptype.charge
end


function ResetParamsIndices()
    last_ptype_ind[] = 0
    # last_region_ind = 0
    last_gas_ind[] = 0
end

# Overwrite the basic PARAMS creation function because it breaks type inference
# for some reason!
function PARAMS(; kwds...)
    ResetParamsIndices()

    thetype = PARAMS
    while !isconcretetype(thetype)
        thetype = thetype{thetype.var.ub}
    end

    p = thetype(; kwds...)
end

################################################################################
# ** Simple maths

"""
EpsFromVel(vel, mass)
EpsFromVel(vel, ::PARTTYPE)
EpsFromVel(::PARTICLE)

Calculate the energy of a particle from the velocity and mass.
"""
function EpsFromVel(vel::VelVec, mass::Q(uM))
    mass == Inf && return Inf

    @inbounds 0.5 * mass * (vel[1]^2 + vel[2]^2 + vel[3]^2)
end
EpsFromVel(vel, ptype::PARTTYPE) = EpsFromVel(vel, ptype.mass)
EpsFromVel(part::PARTICLE) = EpsFromVel(part.vel, part.cur_mass)

export EFromETd, ETdFromE
EFromETd(ETd::Number, tot_dens::Number) = ETd == 0 ? nothing : ETd * tot_dens
EFromETd(ETd::Number, tot_dens::Function, dir) = ETd == 0 ? nothing : (pos,time) -> ETd * tot_dens(pos) * dir
EFromETd(ETd::Symbol, tot_dens) = :inhomogeneous
EFromETd(ETd::Nothing, tot_dens) = nothing

ETdFromE(E::AbstractVector, tot_dens::Number) = norm(E) / tot_dens
ETdFromE(E::Symbol, tot_dens) = :inhomogeneous
ETdFromE(E::Nothing, tot_dens) = nothing

export BFromBHx
BFromBHx(BHx::Nothing, n0) = nothing
BFromBHx(BHx, n0) = (BHx == 0) ? nothing : BHx * n0

BHxFromB(B::AbstractVector, n0) = norm(B) / n0
BHxFromB(B::Symbol, n0) = :inhomogeneous
BHxFromB(B::Nothing, n0) = nothing

# TODO This logic is wrong - really should be returning ETd, even if it has to
# be an anonymous function. The "inhomogeneous" thing should only be for the
# analysis printing.
ETdFromE(params::PARAMS) = ETdFromE(HomogeneousProp(params.E_field), TotDensity(params))
BHxFromB(params::PARAMS) = BHxFromB(HomogeneousProp(params.B_field), TotDensity(params))

export wFromT
function wFromT(gas::GAS, pos)
    gas.tmtr === nothing && return 0.0uVel
    return wFromT(gas.m0, GetVal(gas.tmtr, pos))
end
wFromT(mass, tmtr) = sqrt(2*kB*tmtr / mass)

################################################################################
# ** Bin conversions

function EpsBin(eps, eps_bin_grid, cur_bin=1)
    # Short circuit check - as this is the most common from Collision()
    if cur_bin > 0 && eps_bin_grid[cur_bin] < eps < eps_bin_grid[cur_bin+1]
        return cur_bin
    else
        ind = BisectInterpInternal.BisectInsert(eps, eps_bin_grid)
        return ind - 1
    end
end

function UpdateEpsBin(part::PARTICLE, eps_bin_grid::Vector)
    eps = EpsFromVel(part)
    part.eps_bin = EpsBin(eps, eps_bin_grid, part.eps_bin)
end

ZBin(pos::PosVec, z_grid) = ZBin(pos[3], z_grid)
function ZBin(z, z_grid)
    # Don't need to be clever with the bins here since this is only called at initialisation.
    ind = BisectInterpInternal.BisectInsert(z, z_grid)
    return ind - 1
end

RBin(pos::XYZ, r_grid) = RBin(sqrt(pos[1]^2 + pos[2]^2), r_grid)
function RBin(r, r_grid)
    # Don't need to be clever with the bins here since this is only called at initialisation.
    ind = BisectInterpInternal.BisectInsert(r, r_grid)
    return ind - 1
end

# FindRegion(pos::PosVec, region_list) = FindRegion(pos[3], region_list)
# function FindRegion(z, region_list)
#     ind = findfirst(r -> r.left <= z <= r.right, region_list)
#     return region_list[ind]
# end

function ChangePType(part::PARTICLE, new_ptype::PARTTYPE, params::PARAMS)
    part.ptype_ind = new_ptype.ind
    # part.delay_time = new_ptype.delay_lookup()
    # This might not be the best idea in the future, but for now...
    part.vel *= 0.

    UpdateParticleLookups(part, params)
end


################################################################################
# ** Misc

@xport function TotDensity(params::PARAMS)
    ρ_list = getproperty.(params.gas_list, :ρ)
    if all(x -> x isa Number, ρ_list)
        return sum(ρ_list)
    else
        return :inhomogeneous
    end
end

export HomogeneousProp
HomogeneousProp(val) = :inhomogeneous
HomogeneousProp(val::Number) = val
HomogeneousProp(val::AbstractVector) = val
HomogeneousProp(val::Nothing) = nothing
# TODO Add in a special type that encodes both nothing and a unit

FindPtypeByName(params::PARAMS, name::String) = FindPtypeByName(params, Symbol(name))
@xport function FindPtypeByName(params::PARAMS, sym::Symbol)
    ind = findfirst(x -> GetSymbol(x) == sym, params.ptype_list)

    if ind == 0
        sym_list = GetSymbol.(params.ptype_list)
        error("Ptype $sym doesn't exist! List is $sym_list.")
    end

    return ind
end
