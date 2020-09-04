using OffsetArrays
using .NonNegLeastSquares

DeltaAngleDistFromMTAndViscosity(energy, sigma_raw ; kwds...) = DeltaAngleDistFromMTAndViscosity(energy, cols(sigma_raw)... ; kwds...)
@xport function DeltaAngleDistFromMTAndViscosity(energy, σ_tot, σ_mt ; kwds...)
    σ_l = [σ_tot (σ_tot - σ_mt)]
    MatchWithDeltas(energy, σ_l ; kwds...)
end
@xport function DeltaAngleDistFromMTAndViscosity(energy, σ_tot, σ_mt, σ_vis ; kwds...)
    σ_l = [σ_tot (σ_tot - σ_mt) (σ_mt - 3/2*σ_vis)]
    MatchWithDeltas(energy, σ_l ; kwds...)
end


using LegendrePolynomials
function MatchWithDeltas(eps_list, sigma_l ; num_grid=size(sigma_l,2), normalise=true)
    @assert size(sigma_l,1) == length(eps_list)

    if normalise
        # Sigma 0 is a normalisation factor
        valid_inds = sigma_l[:,1] .> 0uσ

        # Need to create a new array because of units
        temp = zeros(Float64, size(sigma_l)...)
        temp[valid_inds,:] .= sigma_l[valid_inds,:] ./ sigma_l[valid_inds, 1]
        sigma_l = temp
    end

    sigma_l::Matrix{Float64}
    
    maxl = size(sigma_l,2) - 1
    @assert maxl > 0

    x_locs = LinRange(-1, 1, num_grid)

    coeffsmat = hcat([Pl(x, lmax=maxl) for x in x_locs]...)

    Rvals = zeros(length(eps_list), length(x_locs))

    for epsind in 1:length(eps_list)
        out = nonneg_lsq(coeffsmat, sigma_l[epsind,:])

        check = coeffsmat * out - sigma_l[epsind,:]

        out = vec(out)

        if maximum(abs.(check)) > 1e-5
            println("Changing up num_grid to $num_grid")
            @show out
            @show check
            if num_grid > maxl*10
                @error "Too many grid points required!" epsind maximum(abs.(check))
                error("Too many grid points required!")
            end
            return MatchWithDeltas(eps_list, sigma_l, num_grid=num_grid+1, normalise=false)
        end

        cumout = cumsum(out)

        Rvals[epsind,:] = cumout
    end

    return DELTA_ANGLE_DIST(eps_list, x_locs, Rvals)
end

function AngleFromDeltaRvals(eps, R, eps_list, x_locs, Rvals)
    epsind = searchsortedlast(eps_list, eps)
    if 1 <= epsind <= length(eps_list)-1
        alpha = (eps - eps_list[epsind]) / (eps_list[epsind+1] - eps_list[epsind])

        indlow = searchsortedfirst(Rvals[epsind,:], R)
        indhigh = searchsortedfirst(Rvals[epsind+1,:], R)

        indlow = min(indlow, length(x_locs))
        indhigh = min(indhigh, length(x_locs))

        return (1-alpha)*x_locs[indlow] + alpha*x_locs[indhigh]
    elseif epsind == 0
        indhigh = searchsortedfirst(Rvals[1,:], R)
        indhigh = min(indhigh, length(x_locs))

        return x_locs[indhigh]
    else
        indlow = searchsortedfirst(Rvals[end,:], R)
        indlow = min(indlow, length(x_locs))

        return x_locs[indlow]
    end
end
