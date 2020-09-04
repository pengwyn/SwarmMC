
using .MAT

@xport function ReadMaddyCSSet(filename)
    m = matopen(filename)
    data = read(m, "Nu")

    mass = data["m0"] / electronmass

    cs_list = map(zip(data["energy"], data["xsect"], data["type"], data["name"], data["energyThresh"], data["nD_excit"])) do (en,cs,typ,name,threshold,degen)
        manual_degen = 0.
        if typ == 1
            colltype = CFS_ELASTIC()
        elseif typ == 2
            if startswith(name, "Rotation")
                colltype = CFS_INELASTIC_MANUAL_DEGEN()
                manual_degen = degen
            else
                colltype = CFS_INELASTIC()
            end
        elseif typ == 4
            colltype = CFS_IONISATION()
        end
        
        CSData(name, threshold, colltype, en[:], cs[:], manual_degen)
    end

    cs_list = vec(cs_list)

    return cs_list, mass
end
