
println("Be aware - this file doesn't run anything! It needs to be included")

using SwarmMC
function CollectMyResults()
    for F in [0., 0.5, 1.0]
        local params,props
        try
            params,props = ReadAll("LucasSaeleeTransient:F=$(F)")
        catch exc
            if exc isa ErrorException
                println("Ignoring error reading for F=$F.")
                continue
            else
                rethrow()
            end
        end
            

        T = params.time_unit
        L = params.len_unit

        n0 = TotDensity(params) / L^3
        dt = DT(params)

        data = hcat(params.t_grid.grid[2:end]*T*n0,
                props[:eps],
                props[(:vel,:z)]*L/T,
                props[:count]/props.num_particles,
                props[(:pos,:z)]*L*n0,
                props[(:posvel,:z)]*L^2/T*n0,
                props[(:possqr,:z)]*L^2 * n0^2,

                props[:time, (:cum,:t)] * T * n0,
                props[:eps, (:cum,:t)],
                props[(:vel,:z), (:cum,:t)]*L/T,
                props[:duration, (:cum,:t)]/props.num_particles./dt,
                props[(:pos,:z), (:cum,:t)]*L*n0,
                props[(:posvel,:z), (:cum,:t)]*L^2/T*n0,
                props[(:possqr,:z), (:cum,:t)]*L^2 * n0^2,

                props[:cf_1_rate, (:cum,:t)]/T/n0,
                props[:cf_1_deps, (:cum,:t)],
                props[(:cf_1_dvel,:z), (:cum,:t)],

                props[:cf_2_rate, (:cum,:t)]/T/n0,
                props[:cf_2_deps, (:cum,:t)],
                props[(:cf_2_dvel,:z), (:cum,:t)],

                props[:cf_3_rate, (:cum,:t)]/T/n0,
                props[:cf_3_deps, (:cum,:t)],
                props[(:cf_3_dvel,:z), (:cum,:t)],

                props[:field_deps, (:cum,:t)]/T/n0,
                props[(:field_dvel,:z), (:cum,:t)]/T/n0,
                )

        headers = ["time*n0",
                   "energy",
                   "vz",
                   "count/initial",
                   "z*n0",
                   "z*vz*n0",
                   "z^2*n0^2",

                   "meas time * n0 (cum)",
                   "energy (cum)",
                   "vz (cum)",
                   "count/initial (cum)",
                   "z*n0 (cum)",
                   "z*vz*n0 (cum)",
                   "z^2*n0^2 (cum)",

                   "elastic rate / n0",
                   "elastic deps",
                   "elastic dvz",

                   "inelastic rate / n0",
                   "inelastic deps",
                   "inelastic dvz",

                   "ionisation rate / n0",
                   "ionisation deps",
                   "ionisation dvz",

                   "deps rate from field / n0",
                   "dvz rate from field / n0"]

        data = signif.(data, 4)

        open("LS_transient_results_F$(F).txt","w") do file
            println(file, "Everything in m, s, eV: " * join(headers, ", "))
            writecsv(file, data)
        end
    end
end


function CompareWithDale(F, var=:eps)

    if F == 0
        Fstr = "0"
    elseif F == 0.5
        Fstr = "0_5"
    elseif F == 1.0
        Fstr = "1"
    end

    Dale1 = "DaleLSTransient/F=$(Fstr)-tn0=10^12.csv"
    Dale2 = "DaleLSTransient/F=$(Fstr)-tn0=10^14.csv"

    Dale_data1 = readcsv(Dale1, skipstart=1)
    Dale_data2 = readcsv(Dale2, skipstart=1)

    Dale_data = [Dale_data1 ; Dale_data2]

    my_data = readcsv("LS_transient_results_F$(F).txt", skipstart=1)

    time = my_data[:,1]
    if var == :eps
        my_y = my_data[:,2]
        Dale_y = Dale_data[:,1]
    elseif var == :count
        my_y = my_data[:,4]
        Dale_y = Dale_data[:,6]
    elseif var == :vz
        my_y = my_data[:,3]
        Dale_y = Dale_data[:,2]
    elseif var == :z
        my_y = my_data[:,5]
        Dale_y = Dale_data[:,3] / 10.
    elseif var == :zsqr
        my_y = my_data[:,7]
        Dale_y = Dale_data[:,4] / 100.
    end
    plot(time, my_y, "o", label="Danny")
    plot(time, Dale_y, "x", label="Dale")
end

function AllDalePlots(F)
    for var in [:eps, :count, :vz, :z, :zsqr]
        figure()
        CompareWithDale(F, var)
        if var in [:z, :zsqr]
            title("Comparison of $(string(var)) (factors of 10 added for density)")
        else
            title("Comparison of $(string(var))")
        end
        xscale("log")
        yscale("log")
        legend()
        savefig("DaleLSTransient/compare_F$(F)_$(string(var)).png")
        
    end
end
