
################################################
# * Collision frequencies
#----------------------------------------------



const NULL_IND = 0

function EvalCollFreqs(eps, params::PARAMS, ptype_ind::Int, gas_ind::Int)
    cf_list = params.cf_matrix[gas_ind, ptype_ind]
    map(cf_list) do cf
        cf.func(eps)
    end
end

function GasCollFreq(eps, params::PARAMS, ptype_ind::Int, gas_ind::Int)
    cfs = EvalCollFreqs(eps, params, ptype_ind, gas_ind)
    isempty(cfs) && return 0.0*ufreq/uρ
    return sum(cfs)
end
function GasCollFreqWithDensity(pos::PosVec, eps, params::PARAMS, ptype_ind::Int, gas_ind::Int)
    density = GetVal(params.gas_list[gas_ind].ρ, pos)
    return density * GasCollFreq(eps, params, ptype_ind, gas_ind)
end

function TotalCollFreq(pos, eps, params::PARAMS, ptype_ind::Int)
    sum(GasCollFreqWithDensity.(Ref(pos), eps, params, ptype_ind, eachindex(params.gas_list)))
end

function ChooseCollFreq(part::PARTICLE, params::PARAMS)
    @unpack pos,ptype_ind = part
    eps = EpsFromVel(part)

    tcf = TotalCollFreq(pos, eps, params, ptype_ind)
    gases = GasCollFreqWithDensity.(Ref(pos), eps, params, ptype_ind, eachindex(params.gas_list))

    gas_ind = SelectWithNull(tcf, gases)
    gas_ind == NULL_IND && return gas_ind

    # Could just get this by dividing by density...
    gas_tot = GasCollFreq(eps, params, ptype_ind, gas_ind)
    gas_vals = EvalCollFreqs(eps, params, ptype_ind, gas_ind)

    return SelectWithNull(gas_tot, gas_vals)
end

function SelectWithNull(total, vals)
    R = rand() * total

    # Could use searchsortedlast here but it would require cumsum first
    for (ind,val) in enumerate(vals)
        R < val && return ind
        R -= val
    end

    return NULL_IND
end



########################################
# * Integrator things
#--------------------------------------

using OrdinaryDiffEq

function UnpackVector(u)
    pos = PosVec(u[SA[1,2,3]] * uL)
    vel = VelVec(u[SA[4,5,6]] * uVel)
    coll_C = u[7]

    (; pos, vel, coll_C)
end


function dfdt!(u,p,t)
    @unpack pos,vel,coll_C = UnpackVector(u)
    @unpack params,part = p

    eps = EpsFromVel(vel, part.cur_mass)
    acc = zero(XYZ{Q(uF)})
    if params.E_field !== nothing
        acc += GetVal(params.E_field, pos, t)
    end
    if params.B_field !== nothing
        B = GetVal(params.B_field, pos, t)
        acc += vel × B
    end
    acc *= p.part.cur_charge / part.cur_mass

    tcf = TotalCollFreq(pos, eps, params, part.ptype_ind)

    p.num_evals += 1

    du = SVector{7,Float64}(ustrip.(uVel, vel)...,
                            ustrip.(uAcc, acc)...,
                            -ustrip(ufreq,tcf))

    if IsTrue(params.has_cum_meas)
        @unpack weight,log2fac = part
        du_meas = dfdtMeas(params.meas_quants, pos,vel,eps,weight,log2fac)
        return [du ; du_meas]
    else
        return du
    end
end

@generated function dfdtMeas(quants, pos,vel,eps,weight,log2fac)
    block = []
    prep_steps = []

    for quant in quants.parameters
        label,unit,kind = quant.parameters

        sym = gensym()
        push!(prep_steps,
              :($sym = ustrip.(InstValue($(Val(label)), weight, log2fac, pos, vel, eps))))

        if kind == Float64
            push!(block, sym)
        elseif kind == XYZ
            push!(block, :($sym[1]))
            push!(block, :($sym[2]))
            push!(block, :($sym[3]))
        end
    end

    Expr(:block,
         prep_steps...,
         Expr(:call, :SVector, block...))
end

NumScalars(quants::Tuple) = NumScalars(typeof(quants))
NumScalars(quants::Type) = sum(NumScalars, tuple(quants.parameters...))
NumScalars(quant::Type{<:MEAS_QUANT{LABEL,UNIT,Float64}}) where {LABEL,UNIT} = 1
NumScalars(quant::Type{<:MEAS_QUANT{LABEL,UNIT,XYZ}}) where {LABEL,UNIT} = 3

@generated function UnpackMeas(quants, u)
    block = []
    
    ind = 8
    for quant in quants.parameters
        LABEL,UNIT,KIND = quant.parameters

        ftype = FullTypeCumul(KIND, UNIT)
        u = oneunit(oneunit(UNIT) * uN * uT)

        if KIND == Float64
            push!(block, :($ftype(u[$ind] * $u)))
            ind += 1
        elseif KIND == XYZ
            push!(block, :($ftype(((u[$ind], u[$ind+1], u[$ind+2]) .* $u)...)))
            ind += 3
        end
    end

    :(@inbounds tuple($(block...)))
end

function CollisionCond(u,t,int)
    -u[7]
end
function CollisionAffect!(int)
    @unpack part,params,extra_part_list = int.p

    @unpack pos,vel = UnpackVector(int.u)
    @pack! part = pos,vel
    part.time = int.t*uT

    cf_ind = ChooseCollFreq(part, params)

    if cf_ind == NULL_IND
        RecordEvent!(params, Event_null)
    else
        RecordEvent!(params, Event_coll ; part)

        # before_coll_eps = part.eps
        # before_coll_vel = part.vel
        # before_coll_weight = part.weight

        # Need to save this off, in case HandleCollFreq changes the ptype on us.
        this_cf = params.cf_list[part.ptype_ind][cf_ind]

        # FIXME: If HandleCollFreq wants to change the ptype, then it should do so by killing this particle and creating a new one.

        old_eps_bin = part.eps_bin

        HandleCollFreq(this_cf, params, part, extra_part_list)
        @unpack pos,vel = part

        # This is horrible so I should probably change it
        has_eps_bins = any(x -> IsTrue(x.with_eps), params.meas_bins)
        if has_eps_bins && old_eps_bin != part.eps_bin
            temp = part.eps_bin
            part.eps_bin = old_eps_bin
            RecordMeasurements!(int)
            part.eps_bin = temp
        end

        if log2(ustrip(uN,part.weight)) + part.log2fac < params.min_log2_weight
            part.weight = 0uN
        end

        if part.weight == 0uN
            terminate!(int)
        end

        #weight_prefac = before_coll_weight * rellog2fac
        # weight_prefac = before_coll_weight

        # if IsTrue(params.has_cum_meas)
        #     thisprops.cf[cf_ind].num              = weight_prefac
        #     thisprops.cf[cf_ind].numsqr           = weight_prefac * weight_prefac
        #     thisprops.cf[cf_ind].eps_cum          = weight_prefac * before_coll_eps
        #     thisprops.cf[cf_ind].vel_cum          = weight_prefac * before_coll_vel
        #     thisprops.cf[cf_ind].deps_cum         = weight_prefac * (part.eps - before_coll_eps)
        #     thisprops.cf[cf_ind].dvel_cum         = weight_prefac * (part.vel - before_coll_vel)
        #     thisprops.cf[cf_ind].dvelsqr_cum      = weight_prefac * (part.vel.^2 - before_coll_vel.^2)
        #     thisprops.cf[cf_ind].depsflux_cum     = weight_prefac * (part.vel*part.eps - before_coll_vel*before_coll_eps)
        #     thisprops.cf[cf_ind].depspressure_cum = weight_prefac * (part.vel.^2*part.eps - before_coll_vel.^2*before_coll_eps)
        #     thisprops.cf[cf_ind].time_cum         = weight_prefac * part.time
        #     thisprops.cf[cf_ind].numchange        = exp2(part.log2fac) * (part.weight - before_coll_weight)
        # end

        # TODO: measurements

        # @info "Collision" this_cf.colltype this_cf.name before_coll_eps part.eps (part.eps - before_coll_eps)
    end

    new_C = GenerateCollC()
    new_u_core = SVector{7,Float64}(ustrip.(uL, pos)..., ustrip.(uVel, vel)..., new_C)
    new_u_others = int.u[StaticArrays.SUnitRange(8,lastindex(int.u))]
    new_u = [new_u_core ; new_u_others]
    set_u!(int, new_u)
    # Forcing the integrator to forget about the event having already been triggered.
    int.event_last_time = 0
end
const CollisionCallback = ContinuousCallback(CollisionCond, CollisionAffect!, nothing, interp_points=0)

##################
function EpsBinCondition(out,u,t,int)
    @unpack vel = UnpackVector(u)
    @unpack params,part = int.p

    eps = EpsFromVel(vel, part.cur_mass)

    eps_min = params.eps_bin_grid[part.eps_bin]
    eps_max = params.eps_bin_grid[part.eps_bin+1]

    @inbounds out[1] = ustrip(uE, eps_min - eps)
    @inbounds out[2] = ustrip(uE, eps - eps_max)
end
function EpsBinAffect!(int, indx)
    RecordMeasurements!(int)

    @unpack part,params = int.p

    @unpack vel = UnpackVector(int.u)
    eps = EpsFromVel(vel, part.cur_mass)
    old_bin = part.eps_bin
    old_eps_min = params.eps_bin_grid[part.eps_bin]
    old_eps_max = params.eps_bin_grid[part.eps_bin+1]

    if indx == 1
        # Down
        part.eps_bin -= 1
        dir = :down
        # Need to update int so it understands which one fired
        int.vector_event_last_time = 2
    else
        # Up
        part.eps_bin += 1
        dir = :up
        int.vector_event_last_time = 1
    end

    eps_min = params.eps_bin_grid[part.eps_bin]
    eps_max = params.eps_bin_grid[part.eps_bin+1]

    RecordEvent!(params, Event_epsbin ; dir, part.vel,
                 part.eps_bin, eps, eps_min, eps_max,
                 old_bin, old_eps_min, old_eps_max,
                 int.event_last_time,
                 int.vector_event_last_time,
                 int.callback_cache)
end
const EpsBinCallback = VectorContinuousCallback(EpsBinCondition, EpsBinAffect!, nothing, 2)

function ZBinCondition(out,u,t,int)
    @unpack pos = UnpackVector(u)
    z = pos[3]

    @unpack params,part = int.p
    z_min = params.z_bin_grid[part.z_bin]
    z_max = params.z_bin_grid[part.z_bin+1]

    @inbounds out[1] = ustrip(uL, z_min - z)
    @inbounds out[2] = ustrip(uL, z - z_max)
end
function ZBinAffect!(int, indx)
    RecordMeasurements!(int)

    if indx == 1
        # Down
        int.p.part.z_bin -= 1
        # Need to update int so it understands which one fired
        int.vector_event_last_time = 2
    else
        # Up
        int.p.part.z_bin += 1
        int.vector_event_last_time = 1
    end
end
const ZBinCallback = VectorContinuousCallback(ZBinCondition, ZBinAffect!, nothing, 2)

function RBinCondition(out,u,t,int)
    @unpack pos = UnpackVector(u)
    rsqr = sum(pos.^2)

    @unpack params,part = int.p

    rsqr_min = params.r_bin_grid[part.r_bin]^2
    rsqr_max = params.r_bin_grid[part.r_bin+1]^2

    @inbounds out[1] = ustrip(uL^2, rsqr_min - rsqr)
    @inbounds out[2] = ustrip(uL^2, rsqr - rsqr_max)
end
function RBinAffect!(int, indx)
    RecordMeasurements!(int)

    if indx == 1
        # Down
        int.p.part.r_bin -= 1
        # Need to update int so it understands which one fired
        int.vector_event_last_time = 2
    else
        # Up
        int.p.part.r_bin += 1
        int.vector_event_last_time = 1
    end
end
const RBinCallback = VectorContinuousCallback(RBinCondition, RBinAffect!, nothing, 2)

function RecordMeasurements!(int)
    @unpack params,props_out,part = int.p

    IsTrue(params.has_cum_meas) || return
    
    UpdateAllMeasCumulative!(int, params, props_out, part)

    new_u_core = int.u[StaticArrays.SUnitRange(1,7)]
    new_u_others = zeros(SVector{NumScalars(params.meas_quants)})
    new_u = [new_u_core ; new_u_others]
    set_u!(int, new_u)
end

@AutoParm mutable struct INTEGRATOR
    params::AUTO
    part::PARTICLE
    extra_part_list::Vector{PARTICLE}
    props_out::PROPS_OUT
    num_evals::Int = 0
end

MyNorm(u::SVector,t) = sqrt(real(sum(abs2, u[SOneTo(7)])))
MyNorm(u,t) = norm(u)

function SetupIntegrator(params)
    # Note: on testing this separation of Setup and Reinit with Hardsphere
    # benchmark, the speedup was noticable but small (10%-ish).

    u0 = if IsTrue(params.has_cum_meas)
        @SVector zeros(7 + NumScalars(params.meas_quants))
    else
        @SVector zeros(7)
    end

    # Unfortunately need something to initialise the empty particle with
    part = PARTICLE()
    props_out = PROPS_OUT(params,1uN)

    int_data = INTEGRATOR(params=params,
                          part=part,
                          extra_part_list=[],
                          props_out=props_out)
    prob = ODEProblem(dfdt!, u0, (ustrip(uT,params.t_grid[begin]), ustrip(uT,params.t_grid[end])), int_data)

    bin_callbacks = []
    for (sym,callback) in [(:with_r, RBinCallback),
                           (:with_z, ZBinCallback),
                           (:with_eps, EpsBinCallback),
                           # (:with_costh, CosthBinCallback),
                           ]
        if any(bin -> IsTrue(bin.cumulative) && IsTrue(getproperty(bin, sym)), params.meas_bins)
            push!(bin_callbacks, callback)
        end
    end
    cb = CallbackSet(CollisionCallback, bin_callbacks..., params.user_callbacks...)

    int = OrdinaryDiffEq.init(prob,
                              # Need Tsit5 for magnetic fields... it's still not
                              # symplectic, but does much better than BS3.
                              # Probably could select BS3 if only a homogeneous
                              # electric field is present, but want to stay
                              # safe.
                              # BS3(),
                              Tsit5(),
                              save_everystep=false,
                              save_end=true,
                              callback=cb,
                              internalnorm=MyNorm,
                              # Need a choice for dt to avoid warnings from a bad dfdt returned from PARTICLE()
                              dt=1.0,
                              reltol=1e-6,
                              abstol=0)
end

function ReinitIntegrator!(int, part, props_out, time_final)
    u0_core = SVector{7,Float64}([ustrip.(uL, part.pos) ; ustrip.(uVel, part.vel) ; part.coll_C])
    u0_others = zeros(SVector{length(int.u)-7})
    u0 = [u0_core ; u0_others]

    int.p.part = part
    empty!(int.p.extra_part_list)
    int.p.props_out = props_out

    reinit!(int, u0,
            t0=ustrip(uT,part.time),
            tf=ustrip(uT,time_final))
    # TODO: Maybe submit a bug for this
    int.event_last_time = 0
end

function DoSolve(int, part)
    sol = solve!(int)
    @unpack pos,vel,coll_C = UnpackVector(int.u)
    @pack! part = pos,vel,coll_C
    part.time = last(sol.t)*uT

    RecordMeasurements!(int)

    dt = last(sol.t) - first(sol.t)
    return dt, int.p.extra_part_list
end

function GoToTime(int, params::PARAMS, part::PARTICLE, props_out::PROPS_OUT, time_final)
    # This is for my sanity. Not allowed to change log2w during a step.
    old_log2w = part.log2fac

    # int = SetupIntegrator(params, part, props_out, time_final)
    ReinitIntegrator!(int, part, props_out, time_final)

    dt,extra_part_list = DoSolve(int, part)

    CheckEnergyBin(params, part)
    CheckZBin(params, part)
    CheckCollC(params, part)

    @assert part.log2fac == old_log2w

    return extra_part_list
    

    # elseif event == Event_region
    #     if region_event == BinEvent_up
    #         newz = params.region_list[part.region_ind].right
    #     elseif region_event == BinEvent_down
    #         newz = params.region_list[part.region_ind].left
    #     else
    #         error("How?")
    #     end

    #     part.pos = XYZ(part.pos[1], part.pos[2], newz)

    #     # Find an appropriate region
    #     old_region = params.region_list[part.region_ind]
    #     new_region = FindNewRegion(part, params)

    #     # TODO:
    #     change_region = true

    #     # - Check for leaving events.
    #     HandleRegionEvent(old_region.leave_event, part, old_region, params, props_out)

    #     if new_region == nothing
    #         #@debug "Particle died so not entering new region: $(part.weight) $(part.log2fac)"
    #         @assert part.weight == 0.
    #     else
        
    #         # - Check for automatic bounce from difference in V.
    #         charge = params.ptype_list[part.ptype_ind].charge
    #         dV = new_region.V - old_region.V
    #         deps = -(dV * charge)

    #         # println("=========")
    #         # @show part.eps

    #         eps_z = 0.5 * m * part.vel[3]^2
    #         if deps < -eps_z
    #             # Gotta bounce
    #             change_region = false
    #             part.vel = XYZ(part.vel[1], part.vel[2], -part.vel[3])
    #         elseif deps == 0
    #             # Do nothing
    #         else
    #             new_eps_z = eps_z + deps
    #             new_vz = flipsign(sqrt(2 * new_eps_z / m), part.vel[3])

    #             part.vel = XYZ(part.vel[1], part.vel[2], new_vz)
    #             SetEps(part, EpsFromVel(part), params.eps_bin_grid)
    #         end

    #         # @show charge, deps
    #         # @show change_region, part.eps, new_region.ind

    #         if change_region
    #             part.region_ind = new_region.ind
    #             UpdateParticleLookups(part, params)

    #             HandleRegionEvent(new_region.enter_event, part, new_region, params, props_out)
    #         end
    #     end



    # elseif event == Event_delay
    #     ptype = params.ptype_list[part.ptype_ind]
    #     HandleDelayEvent(ptype.delay_event, part, params)
    # end

end



##############################
# * Other
#----------------------------
@inline function RecordEvent!(params, event_type ; extras...)
    if IsTrue(params.do_event_meas)
        params.event_counts[event_type] += 1
    end
    if IsTrue(params.show_events)
        @info "Event" event_type extras
    end
    nothing
end

function HasMeas(params, sym)
    any(params.meas_bins) do bin
        IsTrue(bin.cumulative) && IsTrue(getproperty(bin, sym))
    end
end

##############################
# * Checks
#----------------------------

function CheckEnergyBin(params::PARAMS, part::PARTICLE, args...)
    # TODO: In future this should also include other things that turn on eps bins.
    HasMeas(params, :with_eps) || return
    eps_min = params.eps_bin_grid[part.eps_bin]
    eps_max = params.eps_bin_grid[part.eps_bin+1]

    eps = EpsFromVel(part)
    if eps < eps_min*(1-1e-6)
        @warn "Particle energy too low" eps part.eps_bin eps_min eps_max args
        error("Stop")
    end
    if eps > eps_max*(1+1e-6)
        @warn "Particle energy too high" eps part.eps_bin eps_min eps_max args
        error("Stop")
    end
    nothing
end

function CheckZBin(params::PARAMS, part::PARTICLE, args...)
    HasMeas(params, :with_z) || return

    z_min = params.z_bin_grid[part.z_bin]
    z_max = params.z_bin_grid[part.z_bin+1]

    if isfinite(z_min) && z_min != 0.
        err = 1e8*eps(z_min)
    elseif isfinite(z_max) && z_max != 0.
        err = 1e8*eps(z_max)
    else
        err = 0.
    end

    if part.pos[3] < z_min - err
        @warn "Particle z too low" part.pos part.vel part.z_bin z_min z_max args
        error("Stop")
    end
    if part.pos[3] > z_max + err
        @warn "Particle z too high" part.pos part.vel part.z_bin z_min z_max args
        error("Stop")
    end
end

function CheckCollC(params::PARAMS, part::PARTICLE, args...)
    if part.coll_C < 0
        @warn "Collision C went negative" part.coll_C args
        error("Stop")
    end
    nothing
end





##########################################
# * Old needs checking
#----------------------------------------


@enum BinEvent BinEvent_up BinEvent_down BinEvent_neither


function FindNewRegion(part::PARTICLE, p::PARAMS)
    z = part.pos[3]
    for region in p.region_list
        # Skip the current region
        if region.ind == part.region_ind
            continue
        end

        if z ≈ region.left || z ≈ region.right
            return region
        end
    end

    return nothing
end
