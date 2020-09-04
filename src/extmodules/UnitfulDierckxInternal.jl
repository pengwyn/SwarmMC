module UnitfulDierckxInternal

using Unitful
import Dierckx


# TODO: This currently just turns Spline1D into a function producer. Really it should support "derivative" and all the other Spline1D functions.
# To do that, need to extend the Spline1D methods for different types.

function Dierckx.Spline1D(x::AbstractVector{T}, y::AbstractVector ; kwds...) where {T <: Unitful.AbstractQuantity}
    u = unit(T)
    spl = Dierckx.Spline1D(ustrip.(u,x), y ; kwds...)
    func = let spl=spl,u=u
        x -> spl(ustrip(u,x))
    end
end
function Dierckx.Spline1D(x::AbstractVector{Tx}, y::AbstractVector{Ty} ; kwds...) where {Tx <: Unitful.AbstractQuantity, Ty <: Unitful.AbstractQuantity}
    ux = unit(Tx)
    uy = unit(Ty)
    spl = Dierckx.Spline1D(ustrip.(ux,x), ustrip.(uy,y) ; kwds...)

    func = let spl=spl,ux=ux,uy=uy
        x -> spl(ustrip(ux,x))*uy
    end
end

end # module
