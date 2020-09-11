var documenterSearchIndex = {"docs":
[{"location":"#SwarmMC.jl","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"Documentation for SwarmMC.jl","category":"page"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"Simple API list in order of use.","category":"page"},{"location":"#Setup-of-simulation","page":"SwarmMC.jl","title":"Setup of simulation","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"SwarmMC.PARAMS\nSwarmMC.Finalise\nSwarmMC.PARTTYPE\nSwarmMC.GAS\nSwarmMC.CreateCollFreq\nSwarmMC.CF_STYLE_ENUM\nSwarmMC.ISS_ENUM\nSwarmMC.GNS_ENUM\nSwarmMC.MEAS_BIN\nSwarmMC.MEAS_QUANT\nSwarmMC.MakeSaveName","category":"page"},{"location":"#SwarmMC.PARAMS","page":"SwarmMC.jl","title":"SwarmMC.PARAMS","text":"PARAMS(; kwds...)\n\nThe main parameters object. The parameters can be set as keywords or by mutating the structure afterwards. When finished setting the parameters, call Finalise(...) on the object which will set the structure type parameters and fill in the necessary interal-use-only fields.\n\nKeywords that must be set (no defaults):\n\ninit_style - how particles are generated\nt_grid - the time grid used for the simulation and for measurements\nz_grid, r_grid - the grids used for measurements\ngas_list - a list of GAS objects\nptype_list - a list of PARTTYPE objects\ngen_cf - a function that returns the list of collision frequencies for each gas and particle type combination. Function signature (params, gas, ptype) -> ...\n\nKeywords that have defaults:\n\neps_grid - the energy grid used for measurements\nmeas_bins - the bins for measurements\nmeas_quants - the quantities measured\nsteady_state_timefrac - defines the time at which steady-state measurements are taken from.\nE_field - the electric field. Can be nothing, an XYZ vector, or a function   of signature (pos,time) -> ....\nB_field - the magnetic field. Can be nothing, an XYZ vector, or a function   of signature (pos,time) -> ....\nmin_log2_weight - the point at which if a particle's weight is less than 2^min_log2_weight\nweight_reduction - the percentage of the weight that is lost in a loss collision.\nionisation_as_weight - whether to either create a secondary particle in ionisation, or to double the weight. \ngenerate_next_step_style - one of GNS_ENUM types for refreshing the   particle list for each timestep.\nonly_regen_for_too_many - when true, don't ever generate extra particles   when using GNS_REGEN_ALL.\nmax_substepsize - the maximum timestep. Use this when the time grid is too   coarse to handle changing numbers of particles.\nrunaway_threshold - the multiplicative factor to judge when a runaway of particles is occuring, which will stop the simulation.\nstatic_structure_factor - the static structure factor used for elastic collisions\nsave_name - the name for saving, see MakeSaveName\ndo_event_meas - if true, store information in event_counts\nshow_events - log whenever an event occurs.\nhandle_temperature_bias - set this to false to force sampling gas velocities directly (and incorrectly) from a Maxwellian distribtion.\nuser_callbacks - a list of additional callbacks to provide to the ODE solver. These callbacks will directly interface with DifferentialEquations.jl. For reference, the p in the integrator object will be of type INTEGRATOR. \n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.Finalise","page":"SwarmMC.jl","title":"SwarmMC.Finalise","text":"Finalise(params)\n\nTurn a PARAMS object into the most narrow set of types and also fill in all internal-use parts of the object.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.PARTTYPE","page":"SwarmMC.jl","title":"SwarmMC.PARTTYPE","text":"PARTTYPE(;name=:particle, mass=1mₑ, charge=1q)\n\nThe definition of the type of particle, giving a mass and charge.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.GAS","page":"SwarmMC.jl","title":"SwarmMC.GAS","text":"GAS(; name::Symbol=:gas, m0, ρ, tmtr)\n\nThe definition of a gas. The temperature tmtr defaults to zero, which itself must be specified as nothing.\n\nBoth ρ and tmtr can accept a number or a function which itself takes one argument: pos.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.CreateCollFreq","page":"SwarmMC.jl","title":"SwarmMC.CreateCollFreq","text":"CreateCollFreq(params, gas, ptype,                name,                colltype::CF_STYLE_ENUM,                 gencf::GENCF_PARAMS                ;                threshold=0.0uE,                ion_sharing_style=ISS_AFE(),                ion_sharing_ratio=0.5,                is_cross_section=true)\n\nGenerate a COLLFREQ object, should be used for creating the list to be given to gen_cf of the params object.\n\nThe collision type should be chosen from one of those in CF_STYLE_ENUM. The details of the form of the collision frequency should be provided with gencf. Options are:\n\nGENCF_MANUAL - an arbitrary functional form of the cross section\nGENCF_INTERP - interpolate a given list of energy and cross section\nGENCFINTERPLOG - as with GENCFINTERP but interpolate on the log values\n\nSome convenience structs:\n\nGENCF_HARDSPHERE\nGENCF_MAXWELL\nGENCF_ATTACH - for the Ness-Robson attachment benchmark\n\nBy default, these are interpreted as a cross section. If a collision frequency is directly required, then set is_cross_section to false.\n\nFor inelastic cross sections, threshold is important. For ionisation cross sections, ion_sharing_style and ion_sharing_ratio is important.\n\nDCS values can be provided using angle_dist_cum.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.CF_STYLE_ENUM","page":"SwarmMC.jl","title":"SwarmMC.CF_STYLE_ENUM","text":"CF_STYLE_ENUM\n\nThe type of collision frequency process.\n\nThese types are hopefully obvious:\n\nCFS_ELASTIC\nCFS_INELASTIC\nCFS_LOSS\nCFS_IONISATION\n\nThe types that might need some more explanation are:\n\nCFS_NULL - for introducing manually configured null collisions.\nCFSINELASTICMANUAL_DEGEN - an inelastic cross section with manually set ground/excitation populations, as opposed to using the gas temperature.\nCFS_CHANGETYPE - changes the type of particle.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.ISS_ENUM","page":"SwarmMC.jl","title":"SwarmMC.ISS_ENUM","text":"ISS_ENUM\n\nOptions for ionisation energy sharing\n\nISS_FRACTION - specify an exact constant fraction\nISS_AFE - all fractions equiprobable\nISS_FUNC - the fraction is chosen at the time of collision.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.GNS_ENUM","page":"SwarmMC.jl","title":"SwarmMC.GNS_ENUM","text":"GNS_ENUM\n\nThe style for generating or removing particles for the next time bin.\n\nGNS_NOTHING - don't ever generate or remove new particles\nGNSUPDATELOG2FAC - as with GNS_NOTHING but update the log2fac of the particles.\nGNS_DOUBLE - when there are too few particles, double those that are left.\nGNSREGENALL - when the effective number of particles gets too high or low,                 regenerate by sampling the same weights of particles.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.MEAS_BIN","page":"SwarmMC.jl","title":"SwarmMC.MEAS_BIN","text":"MEAS_BIN([:cum][,:ss][,:t][,:r][,:z][,:eps][,:costh])\n\nDefine the times, positions or energies in which measurements are to be divided into. It also selects whether measurements should be cumulative or instantaneous.\n\nThe constructor takes symbols, e.g. MEAS_BIN(:cum,:t) represents cumulative measurements in time bins and MEAS_BIN(:ss,:eps,:r) are instantaneous measurements in energy and radial bins, only for the steady-state times, defined by steady_state_frac in PARAMS.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.MEAS_QUANT","page":"SwarmMC.jl","title":"SwarmMC.MEAS_QUANT","text":"MEAS_QUANT(label::Symbol, unit, kind::Type=Float64)\n\nReferences the quantity that is measured and will be placed into MEAS_BIN. A quantity defined this way should then extend InstMeas using label to implement the measure itself.\n\nThe storage of this quantity is defined by unit and kind. kind can be Float64 or XYZ.\n\n\n\n\n\n","category":"type"},{"location":"#SwarmMC.MakeSaveName","page":"SwarmMC.jl","title":"SwarmMC.MakeSaveName","text":"MakeSaveName(prefix ; separate_dir=true, kwds...)\n\nMake up a filename prefix to save as, which is: \"prefix:kwd1=val1:kwd2=val2:...\"\n\nIf separate_dir is true then prefix all of the filenames with \"runs_prefix/\"\n\n\n\n\n\n","category":"function"},{"location":"#Running-of-simulation","page":"SwarmMC.jl","title":"Running of simulation","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"SwarmMC.Save\nSwarmMC.LoopMaxTime\nSwarmMC.BunchedPropagate\nSwarmMC.PROPS_OUT","category":"page"},{"location":"#SwarmMC.Save","page":"SwarmMC.jl","title":"SwarmMC.Save","text":"Save(params, props ; mmap=false)\n\nFind a new filename and save the params and props into it using JLD2.\n\nChoose whether to save with IOStream or mmap from the keyword.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.LoopMaxTime","page":"SwarmMC.jl","title":"SwarmMC.LoopMaxTime","text":"LoopMaxTime(params, N ; walltime=:auto)\n\nRuns simulations repeatedly until walltime (a Unitful quantity) is reached.\n\nIf walltime is :auto then try to lookup the environment variable SWARMMC_WALLTIME and interpret it as a time in minutes.\n\nN and params are passed through to BunchedPropagate.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.BunchedPropagate","page":"SwarmMC.jl","title":"SwarmMC.BunchedPropagate","text":"BunchedPropagate(params, N ; return_part_list=false)\n\nRun a simulation of N particles using params. Returns a PROPS_OUT object.\n\nIf you pass true to return_part_list then the initial and final particle lists will also be returned.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.PROPS_OUT","page":"SwarmMC.jl","title":"SwarmMC.PROPS_OUT","text":"PROPS_OUT\n\nStructure to store measurements taken during a simulation. This struct supports convenient querying of different quantities for different bins.\n\nThe user should not need to manually create or populate this, but will need to query it to obtain results.\n\n\n\n\n\n","category":"type"},{"location":"#Loading/Combining-results","page":"SwarmMC.jl","title":"Loading/Combining results","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"SwarmMC.ReadAll\nSwarmMC.GetAll\nSwarmMC.CollectUp\nSwarmMC.NameElements\nSwarmMC.PrefixSets\n","category":"page"},{"location":"#SwarmMC.ReadAll","page":"SwarmMC.jl","title":"SwarmMC.ReadAll","text":"ReadAll(prefix)\n\nLoad and combine all params and props from files belonging to a particular prefix, i.e. the part of the filenames before the \"__*.jld2\" part.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.GetAll","page":"SwarmMC.jl","title":"SwarmMC.GetAll","text":"GetAll(names... ; kwds...)\n\nRead all sets of results in the current directory, extracting quantities given by names. Each name can be either:\n\na symbol, in which case a key from Quants is used.\na string, in which case a parameter from the filename prefix is used. See NameElements\na function, f(params,props,prefix,quants) which produces a value.\n\nA name may also be a Pair, name => value, where value is one of the three choices above.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.CollectUp","page":"SwarmMC.jl","title":"SwarmMC.CollectUp","text":"CollectUp()\n\nCombine files from the same prefix into one file.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.NameElements","page":"SwarmMC.jl","title":"SwarmMC.NameElements","text":"NameElements(str) NameElements(str, index::String)\n\nSplit and identify the different parts of a prefix.\n\nIn the second form, return a specific part given by index.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.PrefixSets","page":"SwarmMC.jl","title":"SwarmMC.PrefixSets","text":"PrefixSets()\n\nReturn all the different unique sets of prefixes in the current directory.\n\n\n\n\n\n","category":"function"},{"location":"#Analysing-results","page":"SwarmMC.jl","title":"Analysing results","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"SwarmMC.Quants\nSwarmMC.EpsFromVel\n\nSwarmMC.ΔR\nSwarmMC.RGrid\n\nSwarmMC.TGrid\nSwarmMC.TInstGrid\nSwarmMC.ΔT\n\nSwarmMC.EpsGrid\nSwarmMC.ΔEps\n\nSwarmMC.ZGrid\nSwarmMC.ΔZ\n\nSwarmMC.CosthGrid\nSwarmMC.ΔCosth","category":"page"},{"location":"#SwarmMC.Quants","page":"SwarmMC.jl","title":"SwarmMC.Quants","text":"Quants(params, props ; mbin, SI=true, n0_factors=false, include_errors=true)\n\nGenerate a dictionary of commonly used transport quantities from time averages in props. The type of measurement bin can be selected in mbin, although it must be either an instantaneous or cumulative time-only set of measurements.\n\nn0_factors indicate that the appropriate density-scaling factors (e.g. n0*DT) will be given instead of the raw values.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.EpsFromVel","page":"SwarmMC.jl","title":"SwarmMC.EpsFromVel","text":"EpsFromVel(vel, mass) EpsFromVel(vel, ::PARTTYPE) EpsFromVel(::PARTICLE)\n\nCalculate the energy of a particle from the velocity and mass.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.ΔR","page":"SwarmMC.jl","title":"SwarmMC.ΔR","text":"ΔR(params)\n\nReturn the size of the radial measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.RGrid","page":"SwarmMC.jl","title":"SwarmMC.RGrid","text":"RGrid(params)\n\nReturn the midpoint of the radial measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.TGrid","page":"SwarmMC.jl","title":"SwarmMC.TGrid","text":"TGrid(params)\n\nReturn the midpoint of all time bins, which is the value relevant for cumulative measurements.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.TInstGrid","page":"SwarmMC.jl","title":"SwarmMC.TInstGrid","text":"TInstGrid(params)\n\nReturn the points at which instantaneous measurements correspond to.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.ΔT","page":"SwarmMC.jl","title":"SwarmMC.ΔT","text":"ΔT(params)\n\nReturn the size of the cumulative time measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.EpsGrid","page":"SwarmMC.jl","title":"SwarmMC.EpsGrid","text":"EpsGrid(params)\n\nReturn the midpoint of the energy measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.ΔEps","page":"SwarmMC.jl","title":"SwarmMC.ΔEps","text":"ΔEps(params)\n\nReturn the size of the energy measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.ZGrid","page":"SwarmMC.jl","title":"SwarmMC.ZGrid","text":"ZGrid(params)\n\nReturn the midpoint of the spatial z measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.ΔZ","page":"SwarmMC.jl","title":"SwarmMC.ΔZ","text":"ΔZ(params)\n\nReturn the size of the spatial z measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.CosthGrid","page":"SwarmMC.jl","title":"SwarmMC.CosthGrid","text":"CosthGrid(params)\n\nReturn the midpoint of the costheta measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.ΔCosth","page":"SwarmMC.jl","title":"SwarmMC.ΔCosth","text":"ΔCosth(params)\n\nReturn the size of the costheta measurement bins.\n\n\n\n\n\n","category":"function"},{"location":"#Advanced","page":"SwarmMC.jl","title":"Advanced","text":"","category":"section"},{"location":"#Definining-new-measurement-types","page":"SwarmMC.jl","title":"Definining new measurement types","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"SwarmMC.Log2Fac\nSwarmMC.InstValue","category":"page"},{"location":"#SwarmMC.Log2Fac","page":"SwarmMC.jl","title":"SwarmMC.Log2Fac","text":"Log2Fac(::Val{sym}, log2fac)\n\nWhat log2fac should be assigned to this quantity. Only relevant for measurements that are considering the weights.\n\n\n\n\n\n","category":"function"},{"location":"#SwarmMC.InstValue","page":"SwarmMC.jl","title":"SwarmMC.InstValue","text":"InstValue(::Val{sym}, weight, log2fac, pos, vel, eps)\n\nThe value of the quantity (sym) being measured, given the current particle properties are pos, vel, eps and the particle's weight is given by weight + 2^log2fac.\n\n\n\n\n\n","category":"function"},{"location":"#User-callbacks","page":"SwarmMC.jl","title":"User callbacks","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"SwarmMC.INTEGRATOR","category":"page"},{"location":"#SwarmMC.INTEGRATOR","page":"SwarmMC.jl","title":"SwarmMC.INTEGRATOR","text":"INTEGRATOR\n\nThe parameter object passed to DifferentialEquations.jl to provide information about the simulation. The user should not need to care about this, unless they are trying to implement custom callbacks.\n\nA custom callback will have access to this structure and hence the following properties:\n\nparams - the PARAMS used for the simulation\npart - the current state of the particle\nextra_part_list - any additional particles to pass back to the simulation (e.g. created during ionisation).\nprops_out - the measurements object\n\n\n\n\n\n","category":"type"},{"location":"#TODO","page":"SwarmMC.jl","title":"TODO","text":"","category":"section"},{"location":"","page":"SwarmMC.jl","title":"SwarmMC.jl","text":"[ ] Custom measurements\n[ ] CF + nonconservative measurements","category":"page"}]
}