module ConstantsInternal

using Reexport
@reexport using Unitful, UnitfulAtomic

using AndExport

#################################################################
#################################################################
#################################################################

@xport begin
    const echarge = UnitfulAtomic.e_au
    const eV = Unitful.eV
    const Eₕ = UnitfulAtomic.Eh_au
    const Eh = Eₕ
    const amu = Unitful.u
    const electronmass = UnitfulAtomic.me_au
    const mₑ = electronmass
    const hconst = Unitful.h
    const hbar = UnitfulAtomic.ħ_au
    const ħ = hbar
    const a0 = UnitfulAtomic.a0_au
    const a₀ = a0
    const AvagadrosNumber = 6.02214179E23
    const clight = 1Unitful.c
    const kB = 1.38064852e-23 * u"J/K"
    const eps0 = 8.8541878176E-12 * u"F/m"
    const electron_r0 = 1/(4pi*eps0) * echarge^2 / (electronmass * clight^2)
    const Å = Unitful.Å
end

export amagat, Td, Hx

@unit    amagat    "amg"    Amagat      2.686_7805e25 * Unitful.m^-3    false
@unit    Td        "Td"     Townsend    1e-21 * u"V*m^2"                false
@unit    Hx        "Hx"     Huxley    1e-27 * u"T*m^3"                false

# const localunits = Unitful.basefactors
# function __init__()
#     Unitful.register(Constants)
#     merge!(Unitful.basefactors, localunits)
# end


# Allow for some more convenient commands
import Unitful, Base

Unitful.ustrip(u::Unitful.AbstractQuantity, val) = uconvert(NoUnits, val / u)
Base.round(x::Quantity ; kwds...) = round(unit(x), x ; kwds...)

@xport ispos(x) = x > zero(x)
@xport isneg(x) = x < zero(x)

# Trying a shorter display
string(x::Type{Q}) where {Q <: Quantity} = isconcretetype(x) ? string("Q{", unit(Q), "}") : invoke(string, Tuple{Any}, x)
Base.show_datatype(io::IO,x::Type{<:Quantity}) = isconcretetype(x) ? print(io, "Q{", Unitful.numtype(x), ",", unit(x), "}") : invoke(Base.show_datatype, Tuple{typeof(io),DataType}, io, x)

Base.typeinfo_implicit(::Type{T}) where {T<:Quantity} = isconcretetype(T)

end
