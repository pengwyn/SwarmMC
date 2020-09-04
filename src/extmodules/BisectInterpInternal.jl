module BisectInterpInternal

using ..DanUtilsInternal

export INTERP1D, INTERP2D, Interp
export Interp1D

struct INTERP1D{X,V}
    x::X

    val::V
    function INTERP1D{X,V}(x::X, val::V) where {X,V}
        @assert length(x) == length(val)
        @assert issorted(x)
        @assert hasmethod(+, Tuple{eltype(V),eltype(V)})
        new(x, val)
    end
end
INTERP1D(x::X,val::V) where {X<:AbstractVector,V<:AbstractVector} = INTERP1D{X,V}(x,val)

struct INTERP2D{X,Y,V}
    #x::AbstractVector{T}
    #y::AbstractVector{T}

    #mat::AbstractMatrix{T}

    x::Vector{X}
    y::Vector{Y}
    mat::Matrix{V}

    function INTERP2D{X,Y,V}(x, y, mat) where {X,Y,V}
        @assert size(mat) == (length(y), length(x))
        @assert all(isfinite.(mat))
        new(x, y, mat)
    end
end
INTERP2D(x::AbstractVector{X}, y::AbstractVector{Y}, mat::AbstractMatrix{V}) where {X,Y,V} = INTERP2D{X,Y,V}(x, y, mat)

Interp(int::INTERP1D, x, extrap=Val{:const}) = Interp1D(x, int.x, int.val, extrap)
(obj::INTERP1D)(x, extrap=Val{:const}) = Interp(obj, x, extrap)

function Interp(int::INTERP2D, x, y)
    xind = BisectInsert(x, int.x)
    xind = xind - 1
    xind = min(xind, length(int.x) - 1)
    xind = max(xind, 1)

    yind = BisectInsert(y, int.y)
    yind = yind - 1
    yind = min(yind, length(int.y) - 1)
    yind = max(yind, 1)

    frac_x = (x - int.x[xind]) / (int.x[xind+1] - int.x[xind])
    frac_y = (y - int.y[yind]) / (int.y[yind+1] - int.y[yind])

    return (1-frac_x) * ( (1-frac_y)*int.mat[yind,xind] + frac_y*int.mat[yind+1,xind] ) +
           frac_x * ( (1-frac_y)*int.mat[yind,xind+1] + frac_y*int.mat[yind+1,xind+1] )
end

struct FAKE_VECTOR{T} <: AbstractVector{T}
    func::Function
    len::Int
end
import Base: getindex, size
getindex(obj::FAKE_VECTOR, ind::Int) = obj.func(ind)
size(obj::FAKE_VECTOR) = (obj.len,)


function InterpRevY(int::INTERP2D, x, val)
    xind = BisectInsert(x, int.x)
    if xind == 1
         xind = 1
    elseif xind >= length(int.x)
        xind = length(int.x) - 1
    else
        xind = xind - 1
    end

    frac_x = (x - int.x[xind]) / (int.x[xind+1] - int.x[xind])

    func = yind -> (1-frac_x)*int.mat[yind,xind] + frac_x*int.mat[yind,xind+1]

    yind = BisectInsert(val, FAKE_VECTOR{eltype(int.mat)}(func, length(int.y)))

    if yind == 1
        yind = 1
    elseif yind >= length(int.y)
        yind = length(int.y) - 1
    else
        yind = yind - 1
    end

    frac_y = (val - func(yind)) / (func(yind+1) - func(yind))

    return (1-frac_y)*int.y[yind] + frac_y*int.y[yind+1]
end

ValTrueFalse = Union{Type{Val{false}}, Type{Val{true}}}

#function BisectInsert(x, vec, reversed=Val{false})
# TODO: Maybe replace this with the searchsorted functions?
function BisectInsert(x, vec, reversed::TypeBool=TypeFalse())
    # Assume that vec is sorted, find the ind which x would take after it is
    # inserted.

    N = length(vec)

    if !Bool(reversed)
        @assert vec[end] > vec[1]
    else
        @assert vec[end] < vec[1]
    end

    @inbounds begin
        if !Bool(reversed)
            if x < vec[1]
                return 1
            elseif x > vec[end]
                return N+1
            end
        else
            if x > vec[1]
                return 1
            elseif x < vec[end]
                return N+1
            end
        end

        left = 0
        right = N+1

        while right > left + 1

            point = (right+left)÷2

            if !Bool(reversed)
                if x >= vec[point]
                    left = point
                else
                    right = point
                end
            else
                if x <= vec[point]
                    left = point
                else
                    right = point
                end
            end
        end

        #return left+1
        return right

    end
end

@inline function Interp1D(x, xlist, ylist, ind::Integer)
    # This version uses knowledge of the index to go straight to the answer

    # No checks for out of bounds

    @inbounds begin
        frac = (x - xlist[ind]) / (xlist[ind+1] - xlist[ind])
        yinterp = ylist[ind] + frac*(ylist[ind+1] - ylist[ind])
    end

    return yinterp
end

function Interp1D(x, xlist, ylist, allow_extrap::Union{Bool,Symbol}, reversed)
    throw("You should be calling this with Val{true/false/:zero}!")
end

Opts_allow_extrap = Union{Type{Val{true}}, Type{Val{false}}, Type{Val{:zero}}, Type{Val{:const}}}

# Interp1D(xvec::AbstractVector, args...) = map(x -> Interp1D(x, args...), xvec)
function Interp1D(x, xlist, ylist, allow_extrap::Opts_allow_extrap=Val{true}, reversed::TypeBool=TypeFalse())

    ind = BisectInsert(x, xlist, reversed)

    if ind == 1 || ind == length(xlist)+1
        if allow_extrap == Val{false}
            throw("Not allowing extrapolation")
        elseif allow_extrap == Val{:zero}
            return 0.
        elseif allow_extrap == Val{:const}
            if ind == 1
                @inbounds return ylist[1]
            else # ind == length(xlist)+1
                @inbounds return ylist[end]
            end
        elseif ind == 1
            ind = 1
        else  #ind == length(xlist)+1
            ind = length(xlist)-1
        end
    else
        ind = ind - 1
    end

    return Interp1D(x, xlist, ylist, ind)
end



# function TestViewBisectNormal(x,y,N)
#     for i in 1:N
#         normalind = BisectInsert(5.5, x)
#     end
# end

# x = collect(1.:100000.)
# y = rand(100000)
# function TestViewBisectView(x,y,N)
#         viewx = BisectInsert(5.5, view(x, 2:10))
#     for i in 1:N
#         #viewind = BisectInsert(5.5, view(x, 2:10))
#         viewind = BisectInsert(5.5, viewx)
#     end
# end

# Module end
end
