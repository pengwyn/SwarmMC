module DanUtils

################################################################
# * Default imports that I expect
#--------------------------------------------------------------
using Reexport
using AndExport
# @reexport using AutoParameters
# @reexport using DelimitedFiles
using Printf
# @reexport using LinearAlgebra
# @reexport using Statistics
# @reexport using ArgCheck
# @reexport using EllipsisNotation
# # @reexport using AutoNamedTuples
# @reexport using IterTools
# @reexport using MsgWrap
# @reexport using Setfield
# @reexport using ProgressMeter
# @reexport using JLD2
# @reexport using UnPack

# using MsgWrap: FormatTime
# export FormatTime

# using Lazy: @>, @>>, @as
# export @>, @>>, @as

using MacroTools: @capture, @q, block, postwalk

# # Convenient iteration flatten
# # My version of flatten1, which doesn't remove blocks at the top level (this ruins structs)
# function my_flatten1(ex)
#   isexpr(ex, :block) || return ex
#   #ex′ = :(;)
#   ex′ = Expr(:block)
#   for x in ex.args
#     isexpr(x, :block) ? append!(ex′.args, x.args) : push!(ex′.args, x)
#   end
#   return ex′
# end
# iterflatten(ex) = postwalk(my_flatten1, block(ex)).args


##############################
# * TypeEnum
#----------------------------

abstract type TypeEnum end

import Base.==
==(x::T, y::Type{T2}) where {T<:TypeEnum,T2} = error("Can't compare a TypeEnum instance with a type")
==(y::Type{T2}, x::T) where {T<:TypeEnum,T2} = (x == y)

function TypeEnum(name::Symbol, symbols::Tuple{Vararg{Symbol}}, also_export=false)
    code = Expr(:block)

    push!(code.args, :(abstract type $(esc(name)) <: TypeEnum end))

    for symbol in symbols
        push!(code.args, :(struct $(esc(symbol)) <: $(esc(name)) end))
    end

    if also_export
        push!(code.args, :(export $(esc(name))))
        for symbol in symbols
            push!(code.args, :(export $(esc(symbol))))
        end
    end
    
    code
end
@xport macro TypeEnum(name, symbols...)
    TypeEnum(name, symbols)
end

@xport macro TypeEnumExport(name, symbols...)
    TypeEnum(name, symbols, true)
end


##############################
# * TypeBool
#----------------------------

import Base: convert, !, promote_rule
@TypeEnumExport TypeBool TypeTrue TypeFalse

@xport const BoolUnion = Union{TypeTrue, TypeFalse, Bool}

Base.convert(T::Type{TypeBool}, x::Bool) = x ? TypeTrue() : TypeFalse()
TypeBool(x::Bool) = convert(TypeBool, x)
Bool(x::TypeTrue) = true
Bool(x::TypeFalse) = false

Base.:!(x::TypeTrue) = TypeFalse()
Base.:!(x::TypeFalse) = TypeTrue()
#Base.promote_rule(x::Type{T}, y::Type{Bool}) where {T<:TypeBool} = Bool

export IsTrue
IsTrue(x::Val{true}) = true
IsTrue(x::Val{false}) = false
IsTrue(x::Bool) = x
IsTrue(x) = Bool(x)
    
#######################################################################
### The standard error. 
###
@xport StdErr(vec) = std(vec) / sqrt(length(vec))
# TODO: Should make a version of this that takes into account autocorrelation.


#######################################################################
### A generic copy of all fields 
###
export GenericCopy
@generated function GenericCopy(other::T) where T
    ex = Expr(:block)

    contents_expr = Expr(:tuple)
    for fname in fieldnames(other)
        push!(contents_expr.args, :(other.$fname))
    end
    
    push!(ex.args, :(contents = $(contents_expr)))
    push!(ex.args, :(obj = T(contents...)))

    ex
end


#######################################################################
### My favourite function as a macro
###
@xport macro printagain(duration::Real=5.)
    a = Ref(0.)

    :( PrintAgain($a, $duration) )    
end
@xport macro printagain(expr::Expr, duration::Real=5.)
    @q begin
        if @printagain $duration
            $(esc(expr))
        end
    end
end

@xport const print_again_override = Ref(false)
@xport function PrintAgain(time_ptr::Ref{Float64}, duration)
    if time() - time_ptr[] > duration || print_again_override[]
        time_ptr[] = time()
        true
    else
        false
    end
end

###############################################################################
### Averages neighbouring points. Useful with diff()
###
import Statistics.mean
@xport function RunningAvg(vec::AbstractVector)
    mean([ vec[1:end-1], vec[2:end] ])
end

@xport function DanMean(vec::AbstractVector, weights::AbstractVector)
    @assert length(vec) == length(weights)

    # Need to split this up as I want 0*NaN = 0 but 1*NaN = NaN
    numerator = sum(w==0 ? 0 : w*v for (v,w) = zip(vec,weights))

    denominator = sum(weights)

    numerator/denominator
end

##############################
# * Check turns
#----------------------------

mutable struct CHECK_TURNS
    tasknum::Int
end

function CheckMyTurn(tracker::CHECK_TURNS)
    #global CheckMyTurn_tasknum

    if ! haskey(ENV,"SLURM_ARRAY_TASK_ID")
        return true
    end

    maxtasks = parse(Int, ENV["SLURM_ARRAY_TASK_MAX"])
    thistask = parse(Int, ENV["SLURM_ARRAY_TASK_ID"]) - 1

    tracker.tasknum += 1

    if tracker.tasknum % maxtasks == thistask
        return true
    end

    return false
end

function CheckTurnsRemaining(tracker::CHECK_TURNS)
    #global CheckMyTurn_tasknum

    if ! haskey(ENV,"SLURM_ARRAY_TASK_ID")
        return false
    end

    maxtasks = parse(Int, ENV["SLURM_ARRAY_TASK_MAX"])

    return tracker.tasknum < maxtasks
end

@xport macro CheckTurns(expr)
    CheckTurns_macrofunc(expr)
end

function CheckTurns_macrofunc(expr)
    tracker = gensym()
    
    expr = deepcopy(expr)
    if expr.head != :block
        expr = Expr(:block, expr)
    end

    # Find lowest for
    inner = expr
    while true
        @assert inner.head == :block
        ind = findall(x -> x isa Expr && x.head == :for, inner.args)
        if length(ind) == 0
            # Do something to add in check
            pushfirst!(inner.args, :($CheckMyTurn($tracker) || continue))
            break
        end
        @assert length(ind) == 1
        ind = ind[1]
        inner = inner.args[ind]
        @assert inner.head == :for
        # Now go to the block inside the for
        inner = inner.args[2]
    end

    # Add in the bit at the end
    push!(expr.args, :($CheckTurnsRemaining($tracker) || break))

    # Also create the counter!
    while_loop = Expr(:while, true, expr)
    container = Expr(:block, :($tracker = $CHECK_TURNS(0)), while_loop)
    
    return esc(container)
end


##############################
# * Other
#----------------------------


################################################################################
# * Array stuff

@xport function Unpack(arr)
    """This is aimed at unpacking vectorised function calls that return tuples.
        And also to be faster than zips."""
    dims = size(arr)
    len = length(arr[1])
    @assert all(length.(arr) .== len)

    out_arr = map(1:len) do ind
        vals = getindex.(arr,ind)
        thistype = promote_type(typeof.(vals)...)

        Array{thistype}(undef, dims)
    end

    for i = eachindex(arr)
        for j = 1:len
            @inbounds out_arr[j][i] = arr[i][j]
        end
    end
    
    tuple(out_arr...)
end

@inline function UnpackBroadcast(func::Function, args...)
    BC = Base.Broadcast
    bc = BC.instantiate(BC.broadcasted(func, args...))
    return_type = BC.combine_eltypes(bc.f, bc.args)

    BroadcastTo(return_type, bc)
end
@inline function UnpackBroadcast(return_type::Type{<:Tuple}, func::Function, args...)
    BC = Base.Broadcast
    bc = BC.instantiate(BC.broadcasted(func, args...))

    BroadcastTo(return_type, bc)
end

BroadcastTo(return_type, bc) = copy(bc)
function BroadcastTo!(return_type, bc)
    @assert isconcretetype(return_type) "Return type $(return_type) is not concrete."
    copyto!(bc)
end

@inline function BroadcastTo(return_type::Type{<:Tuple}, bc)
    # This currently assumes that whatever the return type is, is correct. No
    # narrowing of the type is done.
    @assert isconcretetype(return_type)

    out_arrs = map(tuple(return_type.parameters...)) do param
        similar(bc, param)
    end

    BroadcastTo!(out_arrs, bc)
end

function BroadcastTo!(out_arrs::Tuple, bc)
    # Ignoring this next bit - does it matter?
    #bc′ = BC.preprocess

    @simd for I in eachindex(bc)
        @inbounds out_vals = bc[I]
        foreach(eachindex(out_arrs)) do out_ind
            @inbounds out_arrs[out_ind][I] = out_vals[out_ind]
        end
    end

    return out_arrs
end

#############################
# * Misc

@xport function PrettyTime(secs, sigfigs=2)
    time = secs*1e6
    
    factors = [1000, 1000, 60, 60, 24, Inf]
    cumfactors = cumprod(factors)
    names = ["μs", "ms", "secs", "mins", "hrs", "days"]
    ind = findfirst(time .< cumfactors)

    out = time
    if ind > 1
        out /= cumfactors[ind-1]
    end
    out = round(out, sigdigits=sigfigs)

    return string(out) * " " * names[ind]
end


##############################
# * Break into
#----------------------------
function SeparateBy(iter, separator_func::Function)
    el = eltype(iter)

    Channel(ctype=Vector{el}) do channel
        next = el[]

        for item in iter
            if separator_func(item)
                put!(channel, next)
                next = el[]
            else
                push!(next, item)
            end
        end

        put!(channel, next)
    end
end
function SeparateByHeader(iter, separator_func::Function)
    el = eltype(iter)

    Channel(ctype=Pair{Union{Nothing,el},Vector{el}}) do channel
        header = nothing
        next = el[]

        for item in iter
            if separator_func(item)
                if header != nothing || !isempty(next)
                    put!(channel, header => next)
                end
                next = el[]
                header = item
            else
                push!(next, item)
            end
        end

        put!(channel, header => next)
    end
end
function SeparateByFooter(iter, separator_func::Function)
    el = eltype(iter)

    Channel(ctype=Pair{Union{Nothing,el},Vector{el}}) do channel
        next = el[]

        for item in iter
            if separator_func(item)
                put!(channel, item => next)
                next = el[]
            else
                push!(next, item)
            end
        end

        if !isempty(next)
            put!(channel, nothing => next)
        end
    end
end


##############################
# * FilterTwo
#----------------------------

function FilterTwo(cond, iter)
    el = eltype(iter)
    true_list = el[]
    false_list = el[]

    for item in iter
        if cond(item)
            push!(true_list, item)
        else
            push!(false_list, item)
        end
    end

    return true_list,false_list
end




##########################################
# * Values with errors
#----------------------------------------
@xport struct VALERR{T,T2} <: Number where {T,T2}
    val::T
    sqr::T2 # Note T2 because Unitful will have different units.
    # Allowing for differentally weighted measurements here.
    count::Float64
end
(::Type{T})(val,sqr) where {T<:VALERR} = T(val,sqr,1.0)
(::Type{T})(val) where {T<:VALERR} = T(val,val^2)

@xport Value(x::VALERR) = x.val / x.count
(::Type{T})(x::VALERR) where {T<:Number} = convert(T, Value(x))
Value(x) = x

@xport Comb(x::VALERR, y::VALERR) = VALERR(x.val + y.val, x.sqr + y.sqr, x.count + y.count)
Comb(x::Number, y::Number) = (x + y)/2

@xport Variance(x::VALERR) = (x.sqr/x.count) - (x.val/x.count)^2
# @xport StdDev(x::VALERR) = (v=Variance(x) ; v > 0 ? sqrt(v) : (v > -1e-10*Float64(x) ? 0.0 : error("Negative variance! $v in $(Float64(x))")))
@xport function StdDev(x::VALERR)
    v = Variance(x)
    !isfinite(v) && return v
    v >= zero(v) && return sqrt(v)
    v > -1e-10*Value(x)^2 && return zero(v)
    @error "Negative variance!" v Value(x)
    return 0.
end
@xport StdErr(x::VALERR) = StdDev(x) / sqrt(x.count)

import Base
Base.promote_rule(x::Type{<:VALERR{T}}, y::Type{Ty}) where {T,Ty<:Number} = promote_type(T, y)
Base.:*(x::VALERR, y::Number) = VALERR(x.val*y, x.sqr*y^2, x.count)
Base.:*(x::Number, y::VALERR) = y * x
Base.:/(x::VALERR, y::Number) = VALERR(x.val/y, x.sqr/y^2, x.count)

# For convenience when manipulatingthings that have come from the same number of measurements
# No I can't calculate this!
# Base.:+(x::VALERR, y::VALERR) = VALERR(x.val + y.val, x.sqr/y^2, x.count)

# Base.show(io::IO,x::VALERR) = print(io, round(Value(x), sigdigits=3),"±",round(StdErr(x), sigdigits=3), "(n=$(x.count))")
function Base.show(io::IO,x::VALERR)
    val = Value(x)
    err = StdErr(x)

    val_digits = round(log10(abs(val)), RoundDown)
    err_digits = round(log10(err), RoundDown)
    digits = -min(val_digits, err_digits)
    digits = isfinite(digits) ? Int(digits) : 0
    # print(io, round(Value(x), digits=digits),"±",round(StdErr(x), sigdigits=2), "(n=$(x.count))")
    print(io, round(Value(x), digits=digits),"±",round(StdErr(x), sigdigits=2))
end

Base.parse(::Type{T}, s) where {T <: VALERR} = parse(VALERR{Float64,Float64}, s)
function Base.parse(::Type{T}, s) where {T2, T <: VALERR{T2,T2}}
    v = split(s, "±")
    T(parse.(T2,v)...)
end
    
using RecipesBase
@recipe f(::Type{<:VALERR}, valerr::VALERR) = Value(valerr)

# ** Specialised for Unitful

using Unitful
function Base.show(io::IO,x::VALERR{T,T2}) where {T<:Unitful.AbstractQuantity, T2}
    u = unit(T)
    show(io,VALERR(ustrip(u, x.val), ustrip(u^2, x.sqr), x.count))
    print(io, " ", u)
end
# Need this to disambiguate
Base.:*(x::VALERR, y::Unitful.AbstractQuantity) = invoke(*, Tuple{VALERR,Number}, x, y)
Base.:/(x::VALERR, y::Unitful.AbstractQuantity) = invoke(/, Tuple{VALERR,Number}, x, y)

Unitful.unit(x::VALERR) = unit(x.val)
Unitful.dimension(x::VALERR) = dimension(x.val)

Unitful.unit(x::Type{VALERR{T,T2}}) where {T,T2} = unit(T)
Unitful.dimension(x::Type{VALERR{T,T2}}) where {T,T2} = dimension(T)


Unitful.uconvert(u::Unitful.Units, x::VALERR) = VALERR(uconvert(u,x.val), uconvert(u^2,x.sqr), x.count)
Unitful.ustrip(u::Unitful.Units, x::VALERR) = VALERR(ustrip(u,x.val), ustrip(u^2,x.sqr), x.count)
# This is a copy paste from Unitful, because I can't have VALERR <: Quantity unfortunately.
Unitful.upreferred(x::VALERR) = uconvert(upreferred(unit(x)), x)

##############################
# * Claim filename
#----------------------------

ClaimNextFilename(prefix::AbstractString) = ClaimNextFilename(id -> DefaultClaimFileGen(id, prefix))
function ClaimNextFilename(filename_gen = id -> DefaultClaimFileGen(id, "File"))
    for i in 0:1000
        filename = filename_gen(i)

        try
            FS = Base.Filesystem
            file = FS.open(filename, FS.JL_O_CREAT | FS.JL_O_EXCL, 0o644)
            FS.close(file)
        catch exc
            (exc isa Base.IOError) && exc.code == Base.UV_EEXIST && continue

            rethrow()
        end

        return filename, i
    end

    error("Unable to find free file with prefix $prefix.")
    # return ""
end

DefaultClaimFileGen(id, prefix, suffix=".jld2") = prefix * "__" * gethostname() * @sprintf("%05d", id) * suffix


##############################
# * Clamp names
#----------------------------
@xport ClampMin = max
@xport ClampMax = min

####################################
# * Convenient map!
#----------------------------------
@xport Update!(func, vec) = map!(func, vec, vec)


##############################
# * LogRange
#----------------------------

export LogRange
                        
"""
    LogRange(start, stop, len ; as_val=false, prefactor=1.0)
    
A range with `len` log-spaced elements from `start` to `stop`. By default, will
assume `start` and `stop` are exponents, but if `as_val` is true then they
will be the values themselves.
    
If `prefactor` is not 1.0 then the entire range will be assumed to contain a
prefactor - useful for including, e.g. Unitful quantities, which don't play well
with `exp10` and `log10`.
"""
# Note: Can't make these AbstractRanges because some packages assume this means
# there is a constant stepsize.
struct LogRange{T,T_LIN <: LinRange{T}} <: AbstractVector{T}
    linrange::T_LIN
end
struct LogRangePrefactor{T,T_LIN <: LinRange} <: AbstractVector{T}
    linrange::T_LIN
    prefactor::T
end

function LogRange(start, stop, len ; as_val=false, prefactor=1.0)
    if as_val
        @assert start > zero(start) && stop > zero(stop)
        linrange = LinRange(log10(start), log10(stop), len)
    else
        linrange = LinRange(start, stop, len)
    end

    if prefactor != 1.0
        T = float(promote_type(typeof(prefactor * exp10(start)), typeof(prefactor * exp10(stop))))
        LogRangePrefactor{T,typeof(linrange)}(linrange, prefactor)
    else
        T = float(promote_type(typeof(start), typeof(stop)))
        LogRange{T,typeof(linrange)}(linrange)
    end
end

for func in [:length, :isempty, :firstindex, :lastindex, :size]
    @eval Base.$func(r::LogRange) = $func(r.linrange)
    @eval Base.$func(r::LogRangePrefactor) = $func(r.linrange)
end
Base.:(==)(r1::LogRange, r2::LogRange) = r1.linrange == r2.linrange
Base.:(==)(r1::LogRangePrefactor, r2::LogRangePrefactor) = (r1.linrange == r2.linrange) && (r1.prefactor == r2.prefactor)

Base.minimum(r::LogRange) = exp10(minimum(r.linrange))
Base.maximum(r::LogRange) = exp10(maximum(r.linrange))
Base.reverse(r::T) where {T <: LogRange} = T(reverse(r.linrange))

Base.reverse(r::T) where {T <: LogRangePrefactor} = T(reverse(r.linrange), r.prefactor)

Base.getindex(r::LogRange, i::Integer) = exp10(r.linrange[i])
Base.getindex(r::LogRangePrefactor, i::Integer) = r.prefactor * exp10(r.linrange[i])

Base.show(io::IO, r::LogRange) = print(io, "LogRange(", r.linrange, ")")
Base.show(io::IO, r::LogRangePrefactor) = print(io, "LogRange(", r.linrange, " with prefactor ", r.prefactor, ")")

Base.:*(r::LogRange{T,T_LIN}, x::Number) where {T,T_LIN} = LogRangePrefactor{typeof(one(T)*x),T_LIN}(r.linrange, x)
Base.:*(r::LogRangePrefactor{T,T_LIN}, x::Number) where {T,T_LIN} = LogRangePrefactor(r.linrange, x*r.prefactor)

Base.:*(x::Number, r::Union{LogRange,LogRangePrefactor}) = r*x

Base.:/(r::Union{LogRange,LogRangePrefactor}, x::Number) = r*inv(x)


# ** Unitful things
# Need to put this into a Requires block if I make this a separate module.

function LogRange(start::T, stop::T2, len::Integer ; as_val) where {T <: Unitful.AbstractQuantity, T2 <: Unitful.AbstractQuantity}
    @assert unit(start) == unit(stop)
    @assert as_val == true "When giving units, as_val needs to be true"
    u = unit(start)
    start,stop = promote(ustrip(u,start), ustrip(u,stop))
    LogRange(start, stop, len, as_val=true, prefactor=one(start)*u)
end
Base.:*(r::Union{LogRange,LogRangePrefactor}, x::Unitful.Units) = r*1x



end # module
