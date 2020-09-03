using Plots

using SwarmMC

function LogData(x,y)
    ind = (x .> 0) & (y .> 0)
    return x[ind], y[ind]
end

Temps = [0., 77., 293.]

global f, data

if true

    plotvars = [:eps, :flux_W, :flux_DL, :cf_1_rate, :cf_1_dvz, :cf_1_deps, :cf_2_rate, :cf_2_dvz, :cf_2_deps]
    plotlist = Dict(var => plot(title = string(var), reuse=false) for var in plotvars)

    for Temp in Temps, gridsize in [0]#[0, 10^4]
        if gridsize == 0
            setname = "SuperModelSingle_T$(Temp)"
            restrict = Regex(setname)
            exclude = r"gridsize="
        else
            setname = "SuperModelSingle:T=$(Temp).*gridsize=$(gridsize)\$"
            restrict = Regex(setname)
            exclude = r"^$"
        end

        println("Doing set: ", setname)

        data = GetAll(:ETd, plotvars..., restrict=restrict, exclude=exclude, sortby=:ETd, quiet=true, include_prefix=false, include_sortby=false, n0_factors=true)

        # Got to convert the array of Etd to a single value
        data[1,:] = map(x->x[1], data[1,:])
        # Get only the values not the errors.
        data = Matrix{Float64}(data)

        if Temp > 0
            datasuper = GetAll(:ETd, :cf_3_rate, :cf_3_dvz, :cf_3_deps, restrict=Regex(setname), sortby=:ETd, quiet=true, include_prefix=false, include_sortby=false, n0_factors=true)
            datasuper = Matrix{Float64}(datasuper)
            data = [data ; datasuper[2:end,:]]
        else
            data = [data ; zeros(3, size(data,2))]
        end

        for (i,var) in enumerate(plotvars)
            #plot!(plotlist[var], data[1,:], data[1+i,:], label="el/inel=$(elasticT)K/$(inelasticT)K")
            #plot!(plotlist[var], abs.(data[1,:]), abs.(data[1+i,:]), label="el/inel=$(elasticT)K/$(inelasticT)K", xscale=:log10, yscale=:log10)
            x,y = LogData(abs.(data[1,:]), abs.(data[1+i,:]))
            if !isempty(x)
                plot!(plotlist[var], x, y, label="el/$(Temp)K", xscale=:log10, yscale=:log10)
            end
        end

        # Create the vel and energy transfer rates.
        el_dvzdt = data[5,:] .* data[6,:]
        el_depsdt = data[5,:] .* data[7,:]
        inel_dvzdt = data[8,:] .* data[9,:]
        inel_depsdt = data[8,:] .* data[10,:]
        super_dvzdt = data[11,:] .* data[12,:]
        super_depsdt = data[11,:] .* data[13,:]

        data = [ data[1:7,:] ; el_dvzdt' ; el_depsdt' ; data[8:10,:] ; inel_dvzdt' ; inel_depsdt' ; data[11:13,:] ; super_dvzdt' ; super_depsdt']
        
        delim = ",\t"
        
        f = open("SuperModel_data_T$(Temp).csv", "w")
        write(f, "Data for supermodel with T = $(Temp)\n")
        write(f, join(["E (Td)", "eps (eV)", "W (m/s)", "n0*DL (1/ms)", "elrate/n0 (m^3/s)", "el_dvz (m/s per coll)", "el_deps (eV per coll)", "el_dvzdt (m/s * m^3/s)", "el_depsdt (eV * m^3/s)", "inelrate/n0 (m^3/s)", "inel_dvz (m/s per coll)", "inel_deps (eV per coll)", "inel_dvzdt (m/s * m^3/s)", "inel_depsdt (eV * m^3/s)", "superrate/n0 (m^3/s)", "super_dvz (m/s per coll)", "super_deps (eV per coll)", "super_dvzdt (m/s * m^3/s)", "super_depsdt (eV * m^3/s)"], delim) * "\n")

        # Stop weird things with print_shortest
        #data[1,:] = signif.(data[1,:], 5)
        data = round.(data, sigdigits=5)

        writedlm(f, data', delim)
        close(f)


        timeplot = plot(title="Time Plots Check", reuse=false)

        for name in keys(PrefixSets(Regex("SuperModelSingle_T$(Temp)")))
            params,props = ReadAll(name, quiet=true)
            plot!(timeplot, TGrid(params), props[:eps], xscale=:log10, yscale=:log10) 
        end
        savefig(timeplot, "SuperModel_T$(Temp)_timeplots.png")
    end

    for (var,p) in plotlist
        try
            savefig(p, "SuperModel_figs_" * String(var) * ".png")
            display(p)
        catch exc
            print_with_color(:red, "Error in displaying plot for var=$var. $exc", bold=true)
        end
    end
end

if false

    for elasticT in Temps, inelasticT in Temps
        setname = "SuperModel_elT$(elasticT)_inelT$(inelasticT)"
        println("Doing set: ", setname)

        filenames = Regex("$(setname).*_E[0.]*9[0.]*(_|\$)")
        plot(reuse=false)
        SwarmMC.TimePlotAllSets(:eps, restrict=filenames)

        plot!(title=setname)
    end
end
