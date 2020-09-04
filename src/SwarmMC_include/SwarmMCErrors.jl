#################################################################
# * Error analysis

#################################################################
# ** Exp fitting

function TestFitSuitability(times, vals, tol=1e-3 ; show_plots=:failure, raise_errors=false)
    all(vals .== 0) && return true
    
    # TODO: Need to first check to see if a slope of zero can fit the data.
    # If so, then need to treat this as if we are in the steady-state from the beginning.
    
    params = ExpFit(times,vals)
    c,A,lam = params

    dx = A*exp(-lam*times[end])

    ratio = abs(dx / c)

    if c == 0 || ratio < tol
        pass = true
    else
        pass = false
    end

    if show_plots === true || (show_plots == :failure && !pass)
        if !isdefined(:Plots)
            error("Need to have Plots included for this function to work with show_plots!")
        else
            eval(Expr(:import, :Plots))
        end
        fitted = ExpModel(times, params)
        resids = vals - fitted

        Plots.plot(times, vals, show=true)
        Plots.plot!(times, fitted)

        #Plots.plot(times,resids, reuse=false, show=true)
    end

    if !pass && raise_errors
        error("Fit is not suitable! Ratio is $(round(ratio,sigdigits=3)) > $tol")
    end

    return pass
end

        
#################################################################

function GetValErrAndTest(times, vals, ind, weights=ones(size(vals)))
    # FIXME: Should do the suitability check here.
    norm = maximum(weights[ind:end])
    
    tot = sum(vals[ind:end] .* weights[ind:end]) / norm
    sqrs = sum(vals[ind:end].^2 .* weights[ind:end]) / norm
    count = sum(weights[ind:end]) / norm

    return VALERR(tot,sqrs,count)
end

function AutocorrelationFactorWithFunc(times, m)
    # TODO: Should try bootstrap

    # Assumes times is evenly spaced
    dt = times[2] - times[1]
    @assert all(diff(times) .≈ dt)
    
    # Figuring out tau from fft
    function onedt(m,dt)
        tspan = length(m)-dt
        m1 = m[1:tspan]
        m2 = m[end-tspan+1:end]; #@show length(m1), length(m2)
        out = 1/tspan * sum((m1 - mean(m1)) .* (m2 - mean(m2)))
    end

    chi = [onedt(m,dt) for dt in 0:length(m)-1]

    # Going to only take the first half, since the lengths get too small otherwise
    #intchi = trapz(dt*(0:length(chi)-1), chi)
    intchi = DanUtilsInternal.trapz(dt*(0:length(chi)÷2-1), chi[1:length(chi)÷2])
    tauapprox = intchi / chi[1]

    return tauapprox, chi
end
AutocorrelationFactorWithFunc(m) = AutocorrelationFactorWithFunc(1:length(m), m)

function AutocorrelationFactor(args...)

    return AutocorrelationFactorWithFunc(args...)[1]
end

function AutocorrelationFromBinning(m, maxpow=10)

    if !isdefined(:Plots)
        error("Need to have Plots included for this function to work with show_plots!")
    else
        eval(Expr(:import, :Plots))
    end

    themean = mean(m)

    std_list = [std(m)]

    Plots.plot(1:length(m), m)
    
    for pow = 2:maxpow
        @show pow
        n = 2^pow
        if 2n >= length(m)
            break
        end
        
        binned = [mean(m[(i-1)*n+1:i*n]) for i in 1:(length(m)÷n)]
        push!(std_list, std(binned) * sqrt(n))

        Plots.plot!((1:length(binned)) * n - n/2, binned)
    end

    std_list
end
