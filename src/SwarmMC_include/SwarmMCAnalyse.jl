
using Match

import DataStructures: OrderedDict
@xport function Quants(params::PARAMS, props::PROPS_OUT ; mbin::MEAS_BIN=first(props.meas).mbin, SI=true, n0_factors=false, include_cf=false, steady_state=true, n0_scale=1., include_errors=true)
    params = Finalise(params)

    if IsTrue(mbin.cumulative)
        times = TGrid(params)
        weights = DT(params)
    else
        times = TInstGrid(params)
        weights = ones(length(times))
    end

    tavg = RunningAvg(times)
    diffweights = RunningAvg(weights)

    if steady_state
        ind = findfirst(t -> t > params.steady_state_timefrac * times[end], times)
    else
        ind = 1
    end

    DoItWithErrors(vals, w=weights) = GetValErrAndTest(times, vals, ind, w)
    DoItWithoutErrors(vals, w=weights) = DanMean(vals[ind:end], w[ind:end])
    DoIt = include_errors ? DoItWithErrors : DoItWithoutErrors

    DoIt2WithErrors(vals, w=diffweights) = GetValErrAndTest(tavg, vals, ind, w)
    DoIt2WithoutErrors(vals, w=diffweights) = DanMean(vals[ind:end], w[ind:end])
    DoIt2 = include_errors ? DoIt2WithErrors : DoIt2WithoutErrors

    quants = OrderedDict()

    quants[:E] = HomogeneousProp(params.E_field)
    quants[:ETd] = ETdFromE(params)
    quants[:B] = HomogeneousProp(params.B_field)
    quants[:BHx] = BHxFromB(params)
    quants[:n0] = TotDensity(params)

    if n0_factors && !(quants[:n0] isa Number)
        error("Can't use n0 factors with inhomogeneous n0s")
    else
        n0 = quants[:n0] .* n0_scale
    end
    if quants[:ETd] isa Number
        quants[:ETd] /= n0_scale
    end

    for (label,n0fac,prefac,val) in [
                 ["walltime", 0, 1, props[:walltime]],
                 ["num_init", 0, 1, props[:num_init]],
                 [:weight, 0, 1, :weight],
                 [:log2weight, 0, 1, :log2weight],

                 [:energy, 0, 1, :energy],
                 [:flux_W, 0, 1, (:vel,:z)],

                 [:flux_DTx, 1, 1, (:pos_vel,:x)],
                 [:flux_DTy, 1, 1, (:pos_vel,:y)],
                 ["flux_DT", 0, 1, () -> Comb(quants[:flux_DTx], quants[:flux_DTy])],
                 [:flux_DL, 1, 1, props -> props[:pos_vel, :z, mbin] - props[:pos,:z,mbin].*props[:vel,:z,mbin]],

                 [:epsflux, 0, 1, (:energy_vel, :z)],
                 [:epspressure, 0, 1, (:energy_vel_vel, :z)],
                 [:velsqrx, 0, 1, (:vel_pow2, :x)],
                 [:velsqry, 0, 1, (:vel_pow2, :y)],
                 [:velsqrz, 0, 1, (:vel_pow2, :z)]]

        try
            if label isa Symbol
                if val isa Function
                    quants[label] = prefac*DoIt(val(props))
                elseif val isa Tuple
                    quants[label] = prefac*DoIt(props[val..., mbin])
                else
                    quants[label] = prefac*DoIt(props[val, mbin])
                end
            else
                # Note: be careful of double application of n0 factors.
                label = Symbol(label)
                if val isa Function
                    quants[label] = prefac*val()
                else
                    quants[label] = prefac*val
                end
            end

            if n0_factors
                quants[label] *= n0[]^n0fac
            end
        catch exc
            @info "Ignoring exc for" label exc
        end
    end

    for (label,n0fac,prefac,val) in [
        [:countR, -1, 1, props -> log(2) * (props[:denom, :log2, mbin] - log2.(ustrip.(weights)))],
        ["countR_alt", -1, 1, () -> begin
         w = props[:denom] ./ weights
         DoIt2(diff(w) ./ RunningAvg(w) ./ diff(times))
         end],
        [:bulk_W, 0, 1, (:pos, :z)],
        [:bulk_Wx, 0, 1, (:pos, :x)],
        [:bulk_Wy, 0, 1, (:pos, :y)],
        [:bulk_DTx, 1, 1/2, (:pos_pow2, :x)],
        [:bulk_DTy, 1, 1/2, (:pos_pow2, :y)],
        ["bulk_DT", 0, 1, () -> Comb(quants[:bulk_DTx], quants[:bulk_DTy])],
        [:bulk_DL, 1, 1/2, props -> props[:pos_pow2, :z, mbin] - props[:pos, :z, mbin].^2],
    ]
        try
            if label isa Symbol
                if val isa Function
                    vals = val(props)
                elseif val isa Tuple
                    vals = props[val..., mbin]
                else
                    vals = props[val, mbin]
                end
                quants[label] = prefac*DoIt2(diff(vals) ./ diff(times))
            else
                # Note: be careful of double application of n0 factors.
                label = Symbol(label)
                quants[label] = prefac*val()
            end

            if n0_factors
                quants[label] *= n0[]^n0fac
            end
        catch exc
            @info "Ignoring exc for" label exc
        end
    end

                 # [:alt_bulk_W, 1, alt_bulk_W],
                 # [:alt_bulk_DT, 1, alt_bulk_DT],
                 # [:alt_bulk_DL, 1, alt_bulk_DL],
                 # [:alpha_Townsend_linapprox => alpha_Townsend_linapprox],
                 # [:alpha_Townsend_quadapprox => alpha_Townsend_quadapprox]
    
    # epsz = 0.5*params.ptype_list[1].mass * velsqrz

    # S0 = DoIt(props[:noncons, mbin])

    # S1Tx = DoIt(props[:noncons_pos, :x, mbin])
    # S1Ty = DoIt(props[:noncons_pos, :y, mbin])
    # S1L = DoIt(props[:noncons_pos, :z, mbin] - props[:noncons, mbin].*props[:pos, :z, mbin])

    # S2Tx = 0.5 * DoIt(props[:noncons_possqr, :x, mbin] - props[:noncons, mbin].*props[:pos_pow2, :x, mbin])
    # S2Ty = 0.5 * DoIt(props[:noncons_possqr, :y, mbin] - props[:noncons, mbin].*props[:pos_pow2, :y, mbin])
    # S2L = 0.5 * DoIt(props[:noncons_possqr, :z, mbin] - props[:noncons, mbin].*props[:pos_pow2, :z, mbin]
    #                     - 2*props[:pos, :z, mbin].*props[:noncons_pos, :z, mbin] + 2*props[:noncons, mbin].*props[:pos, :z, mbin].^2)

    # flux_DT = mean([flux_DTx, flux_DTy])
    # S1T = mean([S1Tx, S1Ty])
    # S2T = mean([S2Tx, S2Ty])
    # flux_DT = 1/2 * DoIt(props[:pos_vel, :x, mbin] + props[:pos_vel, :y, mbin])
    # S1T = 1/2 * DoIt(props[:noncons_pos, :x, mbin] + props[:noncons_pos, :y, mbin])
    # S2T = 1/2 * 0.5 * DoIt(  props[:noncons_possqr, :x, mbin] - props[:noncons, mbin].*props[:pos_pow2, :x, mbin]
    #                        + props[:noncons_possqr, :y, mbin] - props[:noncons, mbin].*props[:pos_pow2, :y, mbin])


    # meas_eff_alt = GetValErrAndTest(times, props[:meas_eff_alt, mbin], ind)
    # include_errors || (meas_eff_alt = float(meas_eff_alt))
    
    if SI
        for (key,val) in quants
            val isa Number || continue
            quants[key] = @match key begin
                :ETd => uconvert.(u"Td", val)
                :E => uconvert.(u"V/cm", val)
                :num_init => uconvert.(uN, val)
                :weight => uconvert.(uN, val)
                _ => begin
                    val = map(val) do val
                        if any(x -> x isa Unitful.Dimension{:Mass}, typeof(dimension(val)).parameters[1])
                            val = upreferred.(val/1eV) * 1eV
                        else
                            val = upreferred.(val)
                        end
                    end
                end
            end
        end
    end
    
    # alt_bulk_W = Value(flux_W) + Value(S1L)
    # alt_bulk_DT = Value(flux_DT) + Value(S2T)
    # alt_bulk_DL = Value(flux_DL) + Value(S2L)

	# alpha_Townsend_linapprox = Value(S0) / alt_bulk_W

	# if 4*Value(bulk_DL)*Value(S0) > Value(bulk_W)^2
	# 	alpha_Townsend_quadapprox = 0.
	# else
	# 	alpha_Townsend_quadapprox = Value(bulk_W) / (2*Value(bulk_DL)) * (1 - sqrt(1 - 4*Value(bulk_DL)*Value(S0) / Value(bulk_W)^2))
    # end



    # if include_cf
    #     cf_dict_list = Dict[]
    #     # Cf stuff
    #     # Need to handle multiple ptypes
    #     ptype_ind = 1

    #     for cf_ind = 1:params.num_cf[ptype_ind]
    #         rate = DoIt(props[Symbol("cf_$(cf_ind)_rate"), mbin])
    #         ratevec = props[Symbol("cf_$(cf_ind)_rate"), mbin]
    #         deps = GetValErrAndTest(times, props[Symbol("cf_$(cf_ind)_deps"), mbin], ind, ratevec)
    #         dvz = GetValErrAndTest(times, props[Symbol("cf_$(cf_ind)_dvel"), :z, mbin], ind, ratevec)
    #         dvelsqrz = GetValErrAndTest(times, props[Symbol("cf_$(cf_ind)_dvelsqr"), :z, mbin], ind, ratevec)
    #         depsfluxz = GetValErrAndTest(times, props[Symbol("cf_$(cf_ind)_depsflux"), :z, mbin], ind, ratevec)
    #         depspressurez = GetValErrAndTest(times, props[Symbol("cf_$(cf_ind)_depspressure"), :z, mbin], ind, ratevec)

    #         if scale
    #             rate *= 1. / T
    #             dvz *= L / T
    #             dvelsqrz *= (L / T)^2
    #             depsfluxz *= L / T
    #             depspressurez *= (L / T)^2
    #         end

    #         if n0_factors
    #             rate /= n0common
    #         end

    #         cf_dict = Dict(:rate => rate,
    #                        :deps => deps,
    #                        :dvz => dvz,
    #                        :dvelsqrz => dvelsqrz,
    #                        :depsfluxz => depsfluxz,
    #                        :depspressurez => depspressurez)

    #         push!(cf_dict_list, cf_dict)
    #     end

    #     quants[:cf] = cf_dict_list
    # end

    quants
end

using DataFrames
@xport function GetAll(names... ; restrict=r"", exclude=r"^$", sortby=:prefix, n0_factors=false, quiet=false, include_prefix=false, include_sortby=false, quants_kwds=Any[:SI=>true, :n0_factors=>true, :include_errors=>false], prefix_list=:auto, as_df=true)

    if include_sortby === nothing
        if sortby == :prefix
            include_sortby = false
        elseif sortby ∈ names
            include_sortby = false
        else
            include_sortby = true
        end
    end
    include_sortby::Bool
    
    if prefix_list == :auto
        prefix_dict = PrefixSets()
        prefix_list = keys(prefix_dict) |> collect |> sort
    end

    num_rows = length(names)
    include_prefix && (num_rows += 1)
    include_sortby && (num_rows += 1)
        
    mat = Matrix{Any}(undef, num_rows, 0)
	sortvals = []

    for prefix in prefix_list
        if !occursin(restrict, prefix) || occursin(exclude, prefix)
            if !quiet
                println("Skipping $prefix.")
            end
            continue
        end

        local props, params
        try 
            @show prefix
            params,props = ReadAll(prefix, quiet=true)
        catch exc
            if isa(exc, NoValidFiles)
                println("No valid files found for $prefix. Skipping")
                continue
            end
            rethrow()
        end

        quants = try
            Quants(params, props, n0_factors=n0_factors, include_cf=true ; quants_kwds...)
        catch exc
            @error "Can't process quants" prefix exc
            continue
        end

        function GetName(name::Symbol)
            namestr = string(name)
            if startswith(namestr, "cf_")
                _,cfind,cfprop = split(namestr, '_')
                cfind = parse(Int, cfind)

                return quants[:cf][cfind][Symbol(cfprop)]
            else
                val = quants[name]
                if isa(val, AbstractArray)
                    @assert all(val .== val[1])
                    val = val[1]
                end

                return val
            end
        end
        function GetName(func::Function)
            val = func(params, props, prefix, quants)
        end
        function GetName(str::AbstractString)
            val = NameElements(prefix, str)
        end
        function GetName(pair::Pair)
            val = GetName(pair.second)
        end
        
        col = Any[GetName(name) for name in names]
        if include_prefix
            col = [prefix ; col]
        end

        if !(sortby isa Vector)
            sortby = [sortby]
        end

        sortby_val = map(sortby) do sortby
            if isa(sortby, Function)
                return sortby(prefix, quants)
            elseif isa(sortby, AbstractString)
                return NameElements(prefix, sortby)
            elseif sortby == :prefix
                return prefix
            else
                sortby::Symbol
                temp = GetName(sortby)
                if isa(temp, Vector)
                    temp = tuple(sortby_val...)
                end
                return temp
		    end
        end

        push!(sortvals, sortby_val)
        if include_sortby
            col = [string(sortby_val) ; col]
        end

        mat = [mat col]
    end

	ind = sortperm(sortvals)
	mat = mat[:,ind]

    if as_df
        names = map(names) do name
            if name isa Pair
                name.first
            else
                Symbol(string(name))
            end
        end |> collect

        if include_prefix
            names = [:Prefix ; names]
        end
        if include_sortby
            names = [:SortBy ; names]
        end

        # Inefficient - can't be bothered fixing for now.
        as_cols = collect([z...] for z in rows(mat))
        df = DataFrame(as_cols, names) 
    else
        return mat
    end
end

@xport function ShowAll(names... ; outfile=stdout, with_errors=true, signif_digits=4, include_prefix=true, include_sortby=nothing, sortby=:prefix, kwds...)

    if include_sortby == nothing 
        include_sortby = !(sortby == :prefix)
    end

    mat = GetAll(names...; include_prefix=include_prefix, include_sortby=include_sortby, sortby=sortby, kwds...)

    if length(mat) == 0
        println(outfile, "No entries to show")
        return
    end

	# Convert VALERRS to floats to avoid printing errors later
	if ! with_errors
		for i = 1:size(mat, 1)
			if isa(mat[i,1], VALERR)
				mat[i,:] = Value.(mat[i,:])
			end
		end
	end

    header_rows = 0
    include_prefix && (header_rows += 1)
    include_sortby && (header_rows += 1)

    if signif_digits != nothing
        mat[1+header_rows:end,:] .= round.(float.(mat[1+header_rows:end,:]), sigdigits=signif_digits)
    end
    mat = string.(mat)

    maxsizes = mapslices(x->maximum(length(z) for z in x), mat, dims=2)

    # Forcing arrays
    names = collect(names)

	printnames = string.(names)

    if include_prefix
        printnames = ["Prefix", printnames...]
    end
    if include_sortby
        printnames = ["SortBy", printnames...]
    end

    # Also include the header
    header_maxsizes = length.(printnames)
    maxsizes = max.(vec(maxsizes), vec(header_maxsizes))

	printnames = lpad.(printnames, maxsizes)

    println(outfile, join(printnames, ", "))

    for col in 1:size(mat, 2)
		println(outfile, join([lpad(mat[i,col], maxsizes[i]) for i in 1:size(mat,1)], ", "))
    end
end

@xport function NameElementForGetAll(name)
    (params::PARAMS, props::PROPS_OUT, prefix, quants) -> NameElements(prefix, name)
end

function TimePlotAllSets(variable, meas_type=MEAS_TYPE(:cum,:t) ; restrict::Regex=r"", scale=:log, deriv=false, plotkwds...)
    prefix_dict = PrefixSets()

    if !isdefined(Main, :Plots)
        error("Need to have Plots included for this function to work!")
    end
    Plots = Main.Plots

    p = Plots.plot(xscale=scale, yscale=scale ; plotkwds...)

    didone = false
    for prefix in sort(collect(keys(prefix_dict)))
        if !occursin(restrict, prefix)
            println("Skipping $prefix.")
            continue
        end

        didone = true
        
        params,props = ReadAll(prefix)

        if variable isa Function
            variable(params, props)
        elseif variable == :norm_duration
            val = props[:duration, meas_type] ./ DT(params)
        elseif variable == :flux_DL
            val = props[:int_posvel, :z, meas_type] - props[:int_pos, :z, meas_type].*props[:int_vel, :z, meas_type]
        else
            val = props[variable, meas_type]
        end

        time = TGrid(params)

        if deriv
            val = diff(val) ./ diff(time)
            time = mean([time[2:end], time[1:end-1]])
        end
        
        if scale == :log
            val = abs.(val)
        end
        

        Plots.plot!(time, val, label=prefix)
    end

    if !didone
        error("There were no prefixes that matched restriction $restrict.")
    end

    Base.invokelatest(display, (p,))

    p
end

@xport function Dists(params::PARAMS, props::PROPS_OUT, meas_type=(:eps,:cum,:ss), part_ind=1 ; conv_units=false)
    eps = EpsGrid(params)[1:end-1]
    deps = DEps(params)[1:end-1]

    duration = props[:duration, meas_type, part_ind][1:end-1]

    mass = params.ptype_list[part_ind].mass
    vz = props[(:vel,:z), meas_type, part_ind][1:end-1]
    vz2 = props[(:velsqr,:z), meas_type, part_ind][1:end-1]
    vz3 = props[(:vel3,:z), meas_type, part_ind][1:end-1]
    vz4 = props[(:vel4,:z), meas_type, part_ind][1:end-1]

    vz_pow = [duration vz vz2 vz3 vz4]
    vz_pow[:,2:end] .*= duration
    vz_pow[duration .== 0, :] .= 0.
    
    f = Dists(eps, deps, vz_pow, mass)

    if conv_units
        f .*= (params.time_unit / params.len_unit)^3
        eps .*= ustrip(u"J", 1eV)
    end

    return eps,f
end

# using NumericalIntegration
# using OffsetArrays
# # While I'm fixing the other stuff this will be usign the wrong module
# using LegendrePolys
# @xport function Dists(eps, deps, vz_pow, mass)
#     v = sqrt.(2/mass * eps)

#     max_l = size(vz_pow,2) - 1
#     l = OffsetArray(0:max_l, 0:max_l)'

#     # This is just a power at the moment, not really l.
#     cos_pow = OffsetArray(copy(vz_pow), 1:length(eps), 0:max_l)
#     cos_pow ./= (v.^l)

#     f = similar(cos_pow)

#     for l′ in l
#         poly = LegendrePoly(l′)
#         thisf = sum(poly[l2]*cos_pow[:,l2] for l2 in 0:l′)
#         thisf *= (2*l′ + 1) / 4pi

#         f[:,l′] = thisf
#     end

#     int_factor = sqrt.(2*eps/mass^3)
#     f ./= deps
#     f ./= int_factor

#     norm = NumericalIntegration.integrate(eps, f[:,0] .* int_factor * 4pi)
#     f /= norm

#     return f
# end

using DataStructures
@xport function CFCounts(params::PARAMS, props::PROPS_OUT, part_ind=1, meas_type=MEAS_TYPE(:cum,:t))

    params = Finalise(params)

    cfs = params.cf_list[part_ind]

    out = OrderedDict()
    for (cf_ind,cf) in enumerate(cfs)
        name = cf.name
        count = sum(props[Symbol("cf_$(cf_ind)_num"), meas_type])

        out[name] = count
    end

    return out
end

############################################################
# * For the database benchmarks
#----------------------------------------------------------
using Dates, CSV
@xport function SaveBenchmark(filename, desc, df::AbstractDataFrame ; header_filter=x->nothing, header_lookup=Dict(), round_ETd=:auto, error_scaling=:eps, extra_params=[], date=today(), round_sigdigits=6)
    extra_params = String.(extra_params)
    @assert all(extra_params .∈ Ref(names(df)))

    # May modify df so copy it to be sure
    df = deepcopy(df)

    # Check if replacing ETd by a simpler form is allowed.
    if :ETd ∈ names(df)
        if round_ETd == true || (round_ETd == :auto && all(round.(df[!,:ETd], sigdigits=2) .≈ df[!,:ETd]))
            println("Rounding ETd")
            Update!(x -> round.(x,sigdigits=2), df[!,:ETd])
        end
    end

    if error_scaling != false && error_scaling ∈ names(df) && :numinit ∈ names(df)
        vals = df[!, error_scaling]
        relerr = StdErr.(vals) ./ Value.(vals)
        result = relerr .* sqrt.(df[!, :numinit])
        df[!, Symbol(String(error_scaling) * "_scaling")] = result
    end
            

    # Assuming scaled is on
    units_lookup = Dict("T" => ("T", u"K"),
                        "ETd" => ("E/N", u"Td"),
                        "energy" => ("En", u"eV"),
                        "flux_W" => ("Flux W", u"m/s"),
                        "bulk_W" => ("Bulk W", u"m/s"),
                        "flux_DT" => ("Flux n0*DT", u"1/(m*s)"),
                        "bulk_DT" => ("Bulk n0*DT", u"1/(m*s)"),
                        "flux_DL" => ("Flux n0*DL", u"1/(m*s)"),
                        "bulk_DL" => ("Bulk n0*DL", u"1/(m*s)"),
                        "countR" => ("NonCons Rate / n0", u"m^3/s"),
                        "num_init" => ("Initial particles", uN))
    header = map(names(df)) do name
        name = string(name)
        
        out,u,ustr = if name ∈ keys(header_lookup)
            header_lookup[name]
        elseif name ∈ keys(units_lookup)
            out,u = units_lookup[name]
            # Annoying that unit() doesn't accept Units...
            ustr = (u isa Unitful.Units) ? u : unit(u)
            out,u,string(ustr)
        else
            (string(name), NoUnits, "")
        end

        coltype = eltype(df[!,name])
        if u == nothing || !isconcretetype(coltype)
            # Do nothing!
        elseif coltype == Symbol
            @argcheck dimension(u) == NoDims
        else
            @argcheck dimension(coltype) == dimension(u)
            df[!,name] = ustrip.(u,df[!,name])
        end

        if round_sigdigits !== nothing
            # VALERRs have their own rounding figured out
            if coltype <: Number && !(coltype <: VALERR)
                df[!,name] = round.(df[!,name], sigdigits=round_sigdigits)
            end
        end

        temp = header_filter(name)
        out = (temp === nothing) ? "$out ($ustr)" : temp

        if name ∈ extra_params
            out = "@ $out"
        end

        return out
    end 

    # Temporary buffer to align columns - let writedlm do the heavy lifting
    io = PipeBuffer()
    println(io, join(header, ","))
    CSV.write(IOContext(io, :compact=>true), df, append=true)

    mat = readdlm(io, ',', String)
    colsize = maximum(length.(mat), dims=1)
    # Odd way to do this?
    Update!(s -> s * ",", @view mat[:,1:end-1])
    mat = rpad.(mat, colsize .+ 1)
    


    open(filename, "w") do file
        println(file, desc)
        git_commit = read(`git rev-parse HEAD`) |> String |> strip
        git_branch = read(`git symbolic-ref --short HEAD`) |> String |> strip
        println(file, "Danny Cocks, Date: $date, $git_branch branch of code, git commit: $git_commit")
        println(file, "============")
        writedlm(file, mat)
    end

    return mat
end


using Glob
@xport function SaveBenchmark(filename, desc, args... ; getall_kwds=[], kwds...)
    df = GetAll(args... ; as_df=true, quants_kwds=Any[:SI=>true, :n0_factors=>true, :include_errors=>true], getall_kwds...)

    @show eltype(df[!,:energy])
    filenames = glob("*__*.jld2")
    date = maximum(@. Date(unix2datetime(mtime(filenames))))
    SaveBenchmark(filename, desc, df ; date=date, kwds...)
end
