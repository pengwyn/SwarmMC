using SwarmMC


for exctype in ["JustEl", "WithExc", "All"], ETd in ["1.0", "5.0", "10.0"]
    local props, params
    try
        props,params = SwarmMC.ReadAll("ArForAri$(exctype)_E$ETd")
    catch exc
        warn("Couldn't load files for set ArForAri$(exctype)_E$ETd")
        continue
    end
        

    #quants = SwarmMC.Quants(params, props)
    time = props[:time] * params.time_unit
    eps = props[:eps]
    velz = props[(:vel,:z)] * params.len_unit / params.time_unit
    DL = (props[(:posvel,:z)] - props[(:pos,:z)].*props[(:vel,:z)]) * params.len_unit^2 / params.time_unit
    DT = (props[(:posvel,:x)] + props[(:posvel,:y)])/2. * params.len_unit^2 / params.time_unit

    open("Argon_$(exctype)_E$(ETd).csv", "w") do file
        write(file, "Argon with cross section set $(exctype) at E/n0 = $(ETd)Td with n0 = 10^-24 m^-3. time (s), eps (eV), flux_W (m/s), flux_DL (m^2/s), flux_DT (m^2/s)\n")
        writecsv(file, [time eps velz DL DT])
    end

    E,dists = NormalisedDists(params, props)

    open("Argon_$(exctype)_E$(ETd)_dists.csv", "w") do file
        write(file, "Argon with cross section set $(exctype) at E/n0 = $(ETd)Td with n0 = 10^-24 m^-3. Energy (eps), distribution f(energy) normalised such that int(sqrt(en)*f(en)) = 1, <v_z>(en), <v_z^2>(en)... (only up to v_z^1 at the moment)\n")

        writecsv(file, [E dists])
    end

end
