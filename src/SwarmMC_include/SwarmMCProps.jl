
import Base: +,*
function Base.:+(x::PROPS_OUT, y::PROPS_OUT)
    meas = x.meas .+ y.meas

    N = x.num_particles + y.num_particles
    walltime = x.walltime + y.walltime

    PROPS_OUT(meas, walltime, N)
end

function CombineValWithLog2(val_A, log2w_A, val_B, log2w_B)
    if log2w_A == log2w_B
        return (val_A + val_B), log2w_A
    elseif log2w_A == typemin(Int)
        # FIXME: Do I need typemin checks in here?
        return val_B,log2w_B
    elseif log2w_B == typemin(Int)
        return val_A,log2w_A
    elseif log2w_A > log2w_B
        return (val_A + exp2(log2w_B - log2w_A)*val_B), log2w_A
    elseif log2w_B > log2w_A
        return (exp2(log2w_A - log2w_B)*val_A + val_B), log2w_B
    end
    error("Should never get here")
end

function Base.:+(x::T, y::T) where {T <: MEAS_LOG2}
    @argcheck x.label == y.label
    @argcheck x.ptype_ind == y.ptype_ind
    @argcheck x.cf_ind == y.cf_ind

    vals,log2w = DanUtilsInternal.UnpackBroadcast(CombineValWithLog2, x.vals, x.log2w, y.vals, y.log2w)

    T(x.label, x.mbin, x.ptype_ind, x.cf_ind, vals, log2w)
end

function DirToInt(dir::Symbol)
    if dir == :x
        return 1
    elseif dir == :y
        return 2
    elseif dir == :z
        return 3
    end

    error("Unknown direction $dir")
end

import Base: getindex

function getindex(props_out::PROPS_OUT, prop::Symbol, things... ; ptype=1, cf=0)
    # This is the "smart" version, which tries to disambiguate everything.

    vec_index = nothing

    if !isempty(things) && things[1] isa Integer
        vec_index = things[1]
        things = Base.tail(things)
    end

    meas_bin_syms = Symbol[]
    mbin = nothing
    raw_vals = false
    as_log2 = false

    ptype_assigned = false
    cf_assigned = false
    
    for thing in things
        if thing isa MEAS_BIN
            @assert mbin === nothing
            mbin = thing
        elseif thing ∈ [:ss, :cum]
            push!(meas_bin_syms, thing)
        elseif thing ∈ [:tbin, :rbin, :zbin, :epsbin, :costhbin]
            push!(meas_bin_syms, Symbol(string(thing)[begin:end-3]))
        elseif thing ∈ [:x, :y, :z]
            @assert vec_index == nothing
            vec_index = DirToInt(thing)
        elseif thing == :raw
            raw_vals = true
        elseif thing == :log2
            as_log2 = true
        elseif thing isa Integer
            if !ptype_assigned
                ptype = thing
                ptype_assigned = true
            elseif !cf_assigned
                cf = thing
                cf_assigned = true
            else
                @error "What's this integer mean?" things thing
                error("Stop")
            end
        else
            error("Don't know what thing $thing is")
        end
    end

    if !isempty(meas_bin_syms)
        @assert mbin === nothing
        mbin = MEAS_BIN(meas_bin_syms...)
    end

    GetProp(props_out, Val(prop), vec_index, mbin, ptype, cf ; raw_vals=raw_vals, as_log2=as_log2)
end

function FindMBin(props_out, prop, mbin, ptype_ind, cf_ind)
    ind = findfirst(props_out.meas) do m
        (m.label == prop &&
         (mbin === nothing || m.mbin == mbin) &&
         m.ptype_ind == ptype_ind &&
         m.cf_ind == cf_ind)
    end

    if ind === nothing
        # TODO: Check if prop is in there but not the right mbin
        # TODO: Check if ptype_ind is too large
        # TODO: Check if cf_ind is not there
        meases = getproperty.(props_out.meas, :mbin)
        @error "Can't find index" prop mbin ptype_ind cf_ind meases
        error("Bad prop list - need to improve this error message")
    end

    m = props_out.meas[ind]
end

function GetProp(props_out, prop::Val, vec_index, mbin, args... ; raw_vals=false, as_log2=false)
    m = FindMBin(props_out, prop, mbin, args...) 

    if vec_index !== nothing
        vals = getindex.(m.vals, vec_index)
    else
        vals = m.vals
    end

    log2w = m.log2w
        
    if !raw_vals && prop ∉ (Val(:denom), Val(:denom_dt), Val(:sqrdenom), Val(:invdenom))
        mdenom = FindMBin(props_out, Val(:denom), mbin, args...)
        vals = vals ./ mdenom.vals
        # Don't want to overwrite so copy array
        log2w -= mdenom.log2w
    end

    if as_log2
        out = log2.(ustrip.(vals)) .+ log2w
    else
        out = vals .* exp2.(log2w)
    end
    
    out
end

function GetProp(props_out, prop::Val{:walltime}, vec_index, mbin, rest... ; kwds...)
    @assert mbin === nothing
    return props_out.walltime
end

function GetProp(props_out, prop::Val{:num_init}, vec_index, mbin, rest... ; kwds...)
    @assert mbin === nothing
    return props_out.num_particles
end
function GetProp(props_out, prop::Val{:denom_log2}, vec_index, mbin, rest... ; kwds...)
    @assert mbin === nothing
    return props_out.num_particles
end
function GetProp(props_out, prop::Val{:eff_meas}, vec_index, rest... ; kwds...)
    m = FindMBin(props_out, Val(:denom), rest...)
    @assert !IsTrue(m.mbin.cumulative) "Can't work out effective number of measurements for a cumulative measurement"

    wtot = props_out[:denom, :log2, m.mbin]
    wsqr = props_out[:sqrdenom, :log2, m.mbin]

    exp2.(2*wtot .- wsqr)
end

@xport function NormalisedDists(params::PARAMS, props_out::PROPS_OUT, ptype::Int=1, meas_type=MEAS_TYPE(:ss,:eps,:cum))
    params = Finalise(params)

    meas_type_ind = findfirst(params.meas_types, meas_type)
    @assert meas_type_ind != 0

    base = props_out.props[meas_type_ind][2][ptype,:][:duration]

    # Throw away the final bin since it goes to infinty.
    base = base[1:end-1]
    E = EpsGrid(params)[1:end-1]
    dE = diff(params.eps_bin_grid)[1:end-1]

    num = base ./ dE
    norm = DanUtilsInternal.trapz(E, num)
    num /= norm
    num ./= sqrt.(E)
    
    #rest = dist[:,2:end] ./ base

    #return E, [num rest]
    return E,num
end
