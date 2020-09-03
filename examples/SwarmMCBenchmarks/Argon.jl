module Argon

using DanUtils, Constants
@reexport using SwarmMC

const Ar_mass = 72820.74916940005amu

function LoadBiagiCS(params, gas, ptype, just_elastic=false)
	arr = readdlm("Ar_el_extract_RWfilewithBiagiV10.txt", Float64)
	elastic = CreateCollFreq(params, gas, ptype, "Ar_elastic_Biagi", CFS_ELASTIC(), GENCF_INTERP(arr[:,1]*eV, arr[:,2]*u"Å^2"))
	arr = readdlm("Ar_exc1_extract_RWfilewithBiagiV10.txt", Float64)
	exc1 = CreateCollFreq(params, gas, ptype, "Ar_exc1_Biagi", CFS_INELASTIC(), GENCF_INTERP(arr[:,1]*eV, arr[:,2]*u"Å^2"), threshold=11.55eV)
	arr = readdlm("Ar_exc2_extract_RWfilewithBiagiV10.txt", Float64)
	exc2 = CreateCollFreq(params, gas, ptype, "Ar_exc2_Biagi", CFS_INELASTIC(), GENCF_INTERP(arr[:,1]*eV, arr[:,2]*u"Å^2"), threshold=13.0eV)
	arr = readdlm("Ar_exc3_extract_RWfilewithBiagiV10.txt", Float64)
	exc3 = CreateCollFreq(params, gas, ptype, "Ar_exc3_Biagi", CFS_INELASTIC(), GENCF_INTERP(arr[:,1]*eV, arr[:,2]*u"Å^2"), threshold=14.0eV)
	arr = readdlm("Ar_ion_extract_RWfilewithBiagiV10.txt", Float64)
	ion = CreateCollFreq(params, gas, ptype, "Ar_ion_Biagi", CFS_IONISATION(), GENCF_INTERP(arr[:,1]*eV, arr[:,2]*u"Å^2"), threshold=15.7eV, ion_sharing_style=ISS_AFE(), new_ptype=ptype.ind)

    if just_elastic == true
        return [elastic]
    elseif just_elastic == :allexc
        return [elastic, exc1, exc2, exc3]
    elseif just_elastic == false
        return [elastic, exc1, exc2, exc3, ion]
    else
        error("Unknown value for just_elastic: $just_elastic")
    end
end

# function LoadMcEachranCS(anisostyle=:stepfunc)
#     # I think this is broken

# 	gas = GAS("argon", Ar_mass, 1e-6u"Å^-3", Temp)
#     ptype = PARTTYPE()
# 	#arr = readdlm("ARMD6.CS", Float64, skipstart=3)
# 	arr = readdlm("ARMD6.CS.unix", Float64, skipstart=3)
# 	eps = arr[:,1]
# 	elastic = arr[:,4]
# 	mt = arr[:,6]

# 	if anisostyle == :stepfunc
# 		cs = COLLFREQ(CFS_ELASTIC(), "Ar_elastic_Bob", gas, ptype, ("Interp", eps, elastic), 0., true)
# 		# TODO: Assign mt here.
# 	elseif anisostyle == :mtisotropic
# 		cs = COLLFREQ(CFS_ELASTIC(),"Ar_elastic_Bob", gas, ptype, ("Interp", eps, mt), 0., true)
# 	else
# 		error("Unknown anisostyle: $anisostyle.")
# 	end

# 	push!(ptype.cf_list, cs)

#     ptype
# end
	
@xport function SetupParams(ETd=1Td ; cs="Biagi", just_elastic=false, tmtr=90u"K")
    if ETd < 0.1Td
        maxtime = 1e12
    elseif ETd < 1.0Td
        maxtime = 1e11
    else
        maxtime = 1e10
    end
    maxtime *= SwarmMC.uT

    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_ISOTROPIC(1e-3eV)))

    p.eps_grid = 1eV * LogRange(-5, 5, 1001)
    p.t_grid = LinRange(zero(maxtime), maxtime, 101)

    ρ = 1e-6u"Å^-3"
    p.E_field = ETd * ρ * [0,0,1]

    gas = GAS(; name=:argon, m0=Ar_mass, ρ, tmtr)
    electron = PARTTYPE(:electron)

	push!(p.gas_list, gas)
	push!(p.ptype_list, electron)

    if cs == "Biagi"
        if just_elastic == false
            prefix = "ArBiagi"
        elseif just_elastic == true
            prefix = "ArBiagiJustEl"
        elseif just_elastic == :allexc
            prefix = "ArBiagiAllExc"
        else
            error("Unknown just_elastic: $just_elastic")
        end

        p.gen_cf = (args...) -> LoadBiagiCS(args..., just_elastic)
    else
        prefix = "ArBobMT"
        p.gen_cf = (args...) -> LoadMcEachranCS(args..., :mtisotropic)
    end

    p.save_name = MakeSaveName("runs_Argon/$prefix", ETd=ETd, T=tmtr, separate_dir=false)

    p
end

end
