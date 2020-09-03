
@xport function Finalise(p::PARAMS ; convert=true)
    p.finalised && return p

    if convert
        # Need this to prevent a non-typed version being changed.
        p = deepcopy(p)
    end

    println("Doing PARAMS finalisation on worker $(myid()).")

    if length(p.t_grid) == 0
        error("Currently need to have some times set to run simulation.")
    end

    # Explicitly putting in the infinite limits - this helps with later checks
    p.eps_bin_grid = [0.0*uE ; p.eps_grid ; Inf*uE]
    p.z_bin_grid = [-Inf*uL ; p.z_grid ; Inf*uL]
    p.r_bin_grid = [0.0*uL ; p.r_grid ; Inf*uL]

    # fakegas = GAS("fakegas", 0., 1.)
    # if p.adapt_noncons_style == ANS_FAKE_NONCONS()
    #     # Need to have n0 = 1 for the cf stuff to work out.
    #     # This will only work with one particle type now.
    #     thecf = COLLFREQ(CFS_FAKE_NONCONS(), "FakeNoncons", fakegas, p.ptype_list[1], tuple("FakeNoncons"))
    #     p.fake_noncons_cf = thecf
    #     push!(p.ptype_list[1].cf_list, thecf)
    # end

    # Make sure all gases are in their lists
    @assert all(p.ptype_list[ind].ind == ind for ind in eachindex(p.ptype_list))
    @assert all(p.gas_list[ind].ind == ind for ind in eachindex(p.gas_list))
    # @assert all(p.region_list[ind].ind == ind for ind in eachindex(p.region_list))

    p.has_cum_meas = TypeBool(any(x -> IsTrue(x.cumulative), p.meas_bins))

    # for region in p.region_list, gas in p.gas_list
    #     dens = p.chars[:dens, GetSymbol(region), GetSymbol(gas)]
    #     tmtr = p.chars[:temp, GetSymbol(region), GetSymbol(gas)]
    #     w = wFromT(gas, tmtr)

    #     p.density[region.ind, gas.ind] = dens
    #     p.temperature[region.ind, gas.ind] = tmtr
    #     p.gas_w[region.ind, gas.ind] = w
    # end

    p.cf_matrix = Matrix{Vector{COLLFREQ}}(undef, length(p.gas_list), length(p.ptype_list))
    for gas in p.gas_list, ptype in p.ptype_list

        cflist = p.gen_cf(p,gas,ptype)

        # Add on superelastics
        super_list = COLLFREQ[]
        for cf in cflist
            if cf.colltype isa CFS_INELASTIC_ABSTRACT && gas.tmtr !== nothing
                kT = gas.tmtr * kB
                if cf.colltype == CFS_INELASTIC()
                    # This is a crude approximation! Valid only when the
                    # exictations are high enough in energy to not significantly
                    # deplete the ground state.
                    boltzfac = exp(-cf.threshold / kT) / (1 + exp(-cf.threshold / kT))

                    func = let boltzfac=boltzfac, threshold=cf.threshold, oldfunc=cf.com_func
                        eps -> boltzfac * sqrt(1 + threshold/eps) * oldfunc(eps + threshold)
                    end
                elseif cf.colltype == CFS_INELASTIC_MANUAL_DEGEN()
                    # This has to keep the original for the partition function,
                    # but corrects the factor on top.
                    boltzfac = exp(-cf.threshold / kT)

                    func = let boltzfac=boltzfac, threshold=cf.threshold, func=cf.func
                        eps -> boltzfac * sqrt(1 + threshold/eps) * func(eps + threshold)
                    end
                end

                super = COLLFREQ(CFS_SUPERELASTIC(), cf.name*"_Super", cf.ptype_ind, cf.gas_ind, -cf.threshold, cf.manual_degen, cf.angle_dist_cum, cf.ionisation_sharing_ratio, cf.ionisation_sharing_style, cf.new_ptype_ind, func, cf.com_func)

                push!(super_list, super)
            end
        end

        append!(cflist, super_list)

        p.cf_matrix[gas.ind, ptype.ind] = cflist

        for cf in cflist
            if cf.colltype == CFS_IONISATION()
                @assert 1 <= cf.new_ptype_ind <= length(p.ptype_list) || cf.new_ptype_ind == -1
            end
        end
    end
    p.cf_list = reduce(vcat, p.cf_matrix, dims=1, init=COLLFREQ[]) |> x -> dropdims(x, dims=1)

    # @msgwrap "Generating total coll freqs" begin
    #     p.tot_cf_matrix = Matrix{Vector{Q(ufreq)}}(undef, length(p.region_list), length(p.ptype_list))
    #     for region in p.region_list, ptype in p.ptype_list
    #         p.tot_cf_matrix[region.ind, ptype.ind] = GenerateCFTotMax(p, region, ptype)
    #     end
    # end


    # for region in p.region_list
    #     if region.leave_event isa REGION_EVENT_MEASURE
    #         #region.leave_event.meas = CreateArray(p, region.leave_event.meas_type)
    #         name = "region_" * GetName(region) * "_leave"
    #         p.meas_types = (p.meas_types..., MEAS_TYPE(Symbol(name), GetSymbols(region.leave_event.meas_type)...))
    #         region.leave_event.meas_ind = length(p.meas_types)
    #     end
    #     if region.enter_event isa REGION_EVENT_MEASURE
    #         #region.enter_event.meas = CreateArray(p, region.enter_event.meas_type)
    #         name = "region_" * GetName(region) * "_enter"
    #         # p.meas_types = (p.meas_types..., MEAS_TYPE(:region, GetSymbols(region.enter_event.meas_type)...))
    #         p.meas_types = (p.meas_types..., MEAS_TYPE(Symbol(name), GetSymbols(region.enter_event.meas_type)...))
    #         region.enter_event.meas_ind = length(p.meas_types)
    #     end
    # end

    p.finalised = true

    if convert
        return ConvertToTyped(p)
    else
        return p
    end
end

function ConvertToTyped(obj)
    args = [getfield(obj,i) for i = 1:nfields(obj)]    

    T = typeof(obj).name.wrapper
    T(args...)
end

# function GenerateCFTotMax(p::PARAMS, region::REGION, ptype::PARTTYPE)
#     cf_tot_max = zeros(Q(ufreq), length(p.eps_grid) - 1)

#     for ind in 1:length(p.eps_grid)-1
#         eps_list = linspace(p.eps_grid[ind], p.eps_grid[ind+1], 101)
#         cf_tot_max[ind] = 1.1 * maximum(TotalCollFreq.(eps_list, Ref(p), region.ind, ptype.ind))
#     end

#     # We want to easily look up with the eps_bin, so add on the extremities
#     #cf_tot_max = [0.0 ; cf_tot_max ; 0.0]
#     cf_tot_max = [cf_tot_max[1] ; cf_tot_max ; cf_tot_max[end]]

#     return cf_tot_max
# end
