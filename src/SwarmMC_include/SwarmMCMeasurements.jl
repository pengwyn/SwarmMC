
####################################################
# * Cumulative measurements
#--------------------------------------------------

sym_list = [:vel, :pos, :energy]
push!(sym_list, [Symbol(sym1,"_",sym2) for sym1 in sym_list, sym2 in sym_list]...)
push!(sym_list, [Symbol(sym,"_pow",pow) for sym in [:vel, :pos, :energy],
                 pow in [2,3]]...)
push!(sym_list, :energy_vel_vel)

# function CumulValue end
function InstValue end
InstValue(quant::MEAS_QUANT{LABEL}, args...) where LABEL = InstValue(Val(LABEL), args...)
InstValue(prop::Val, part::PARTICLE) = InstValue(prop, part.weight, part.log2fac, part.pos, part.vel, EpsFromVel(part))
# Even though some of these look the same, they change with the log2facs
InstValue(::Val{:denom}, weight, log2fac, pos, vel, eps) = weight
InstValue(::Val{:sqrdenom}, weight, log2fac, pos, vel, eps) = (weight)^2
InstValue(::Val{:invdenom}, weight, log2fac, pos, vel, eps) = 1/weight
InstValue(::Val{:weight}, weight, log2fac, pos, vel, eps) = weight * weight * exp2(log2fac)
InstValue(::Val{:preweight}, weight, log2fac, pos, vel, eps) = weight * weight
InstValue(::Val{:log2weight}, weight, log2fac, pos, vel, eps) = weight * log2fac

Log2Fac(::Val, log2fac) = log2fac
Log2Fac(::Val{:sqrdenom}, log2fac) = log2fac*2
Log2Fac(::Val{:invdenom}, log2fac) = -log2fac

for sym in sym_list
    T = Val{sym}

    pieces = split(string(sym),"_")
    if length(pieces) == 1
        expr = sym
    else
        if startswith(string(last(pieces)), "pow")
            pow = parse(Int, string(last(pieces))[4:end])
            pieces = repeat(pieces[1:end-1], pow)
        end
        expr = :(broadcast(*, $(Symbol.(pieces)...)))
    end

    @eval InstValue(::$T, weight, log2fac, pos, vel, energy) = weight * ($expr)
end

@generated function MakeIndices(x::MEAS_BIN, t_bin, r_bin, z_bin, eps_bin, costh_bin)
    expr = Expr(:tuple)

    for sym in [:t, :r, :z, :eps, :costh]
        if fieldtype(x, Symbol(:with_,sym)) <: TypeTrue
            push!(expr.args, Symbol(sym, :_bin))
        end
    end

    expr
end


@inline function UpdateMeasInst!(meas, part::PARTICLE, steady_state, bins=(part.t_bin, part.r_bin, part.z_bin, part.eps_bin, 0))
    meas.cf_ind == 0 || return
    meas.ptype_ind == part.ptype_ind || return

    IsTrue(meas.mbin.steady_state) && !steady_state && return
    IsTrue(meas.mbin.cumulative) && return

    inds = MakeIndices(meas.mbin, bins...)

    any(inds .== 0) && return

    mass = part.cur_mass

    val = InstValue(meas.label, part)
    this_log2w = Log2Fac(meas.label, part.log2fac)

    meas.vals[inds...],meas.log2w[inds...] = CombineValWithLog2(meas.vals[inds...], meas.log2w[inds...], val, this_log2w)

    nothing
end


@inline function UpdateMeasCumul!(val, meas, part::PARTICLE, steady_state, bins)
    meas.cf_ind == 0 || return
    IsTrue(meas.mbin.steady_state) && !steady_state && return
    IsTrue(meas.mbin.cumulative) || return

    inds = MakeIndices(meas.mbin, bins...)
    any(inds .== 0) && return

    this_log2w = Log2Fac(meas.label, part.log2fac)

    meas.vals[inds...],meas.log2w[inds...] = CombineValWithLog2(meas.vals[inds...], meas.log2w[inds...], val, this_log2w)

    nothing
end

function UpdateAllMeasCumulative!(int, params, props_out, part, bins=(part.t_bin, part.r_bin, part.z_bin, part.eps_bin, 0))
    steady_state = part.time > params.t_grid[end]*params.steady_state_timefrac

    vals = UnpackMeas(params.meas_quants, int.u)

    for q_ind in eachindex(params.meas_quants)
        val = vals[q_ind]
        map(props_out.meas[q_ind,:,part.ptype_ind]) do meas
            UpdateMeasCumul!(val, meas, part, steady_state, bins)
        end
    end
end


function NDims(mtype::MEAS_BIN)
    dims = sum(MakeIndices(mtype,
                           1,
                           1,
                           1,
                           1,
                           1))
end

function Dims(params::PARAMS, mtype::MEAS_BIN)
    dims = MakeIndices(mtype,
                       length(params.t_grid) - 1,
                       length(params.r_bin_grid) - 1,
                       length(params.z_bin_grid) - 1,
                       length(params.eps_bin_grid) - 1,
                       length(params.costh_bin_grid) - 1)
end


function PROPS_OUT(params::PARAMS, N=0uN)
    props = [MEAS_LOG2(params, quant, mbin, ptype.ind, 0)
             for quant in params.meas_quants,
             mbin in params.meas_bins,
             ptype in params.ptype_list]
    
    PROPS_OUT(props, 0u"s", N)
end
