
export ReadCSFile, ReadIlijaCSFile, ReadLXCatCSFile, ReadLXCatSwarmSet

# This is an old version - I don't konw what it's used for now.
function ReadCSFile(filename)
    fullset = []
    
    open(filename) do file
        while !eof(file)
            name = readline(file) |> strip

            threshold = parse(Float64, readline(file))

            en = readline(file) |> strip
            en = [parse(Float64, z) for z in split(en, ',')]

            cs = readline(file) |> strip
            cs = [parse(Float64, z) for z in split(cs, ',')]

            push!(fullset, (name, threshold, en, cs))
        end
    end

    fullset
end


@AutoParm mutable struct CSData
    name::String
    threshold::Q(uE)
    colltype
    energy::Vector{Q(uE)}
    cs::Vector{Q(uσ)}
    m0

    manual_degen::Float64 = 0.
end

function ReadIlijaCSFile(filename)

    lines = readlines(filename)
    lines = strip.(lines)
    filter!(x -> length(x) > 0, lines)

    num_cs = parse(Int, popfirst!(lines))

    global_mass = 0.0

    cs_list = CSData[]
    
    for cs_ind = 1:num_cs
        @assert popfirst!(lines) == "No_$(cs_ind)"
        name = popfirst!(lines)

        # Skip the header
        popfirst!(lines)

        thresholdstr, ndatastr, scalestr, cstypestr, massstr = popfirst!(lines) |> split
        threshold = parse(Float64, thresholdstr)
        ndata = parse(Int, ndatastr)
        scale = parse(Float64, scalestr)
        cstype = parse(Int, cstypestr)
        mass = parse(Float64, massstr)

        if cs_ind != 1
            @assert mass == global_mass
        end
        global_mass = mass

        @assert occursin(r"^-+$", popfirst!(lines))

        csdatastr = join(lines[1:ndata], "\n")
        lines = lines[ndata+1:end]

        csdata = readdlm(IOBuffer(csdatastr))
        energy = csdata[:,1]
        cs = csdata[:,2]


        colltype = "Unknown"
        if cstype == 10
            colltype = "ElasticWTBoth"
        elseif cstype == 11
            colltype = "ElasticWTMomentumOnly"
        elseif cstype == 12
            colltype = "ElasticWTEnergyOnly"
        elseif cstype == 0
            colltype = CFS_INELASTIC()
        elseif cstype == 1
            colltype = CFS_IONISATION()
        end

        push!(cs_list, CSData(name, threshold, colltype, energy, cs))
    end

    @assert length(lines) == 0

    cs_list, global_mass
end

function SkipDatabase(lines)
    @assert occursin(r"DATABASE:", lines[1])
    while !occursin(r"^x+$", popfirst!(lines)) end
end

function SkipHeader(lines)
    @assert occursin(r"^\*{6,}", lines[1])
    while !occursin(r"^\*+\s*\w+\s*\*+$", popfirst!(lines)) end
end

function ReadLXCatCSFile(filename)

    lines = readlines(filename)
    lines = strip.(lines)

    cs_list = CSData[]
    
    filter!(x -> !isempty(x), lines)

    # Go up to the first database definition
    while !occursin(r"^x+$", popfirst!(lines)) end

    while length(lines) > 0
        SkipDatabase(lines)

        while true
            # Skip past cs header
            SkipHeader(lines)
            # In the future, I should make this read the name of t

            while true
                colltypestr = popfirst!(lines)
                name = popfirst!(lines)

                threshold = 0.
                m0 = nothing
                if colltypestr == "ELASTIC"
                    colltype = CFS_ELASTIC()
                    temp = parse(Float64, popfirst!(lines))
                    m0 = 1/temp * electronmass
                elseif colltypestr == "EXCITATION"
                    colltype = CFS_INELASTIC()
                    threshold = parse(Float64, popfirst!(lines))
                elseif colltypestr == "IONIZATION"
                    colltype = CFS_IONISATION()
                    threshold = parse(Float64, popfirst!(lines))
                elseif colltypestr == "ATTACHMENT"
                    colltype = CFS_LOSS()
                else
                    error("Unknown colltype $colltypestr.")
                end

                while true
                    line = popfirst!(lines)
                    if occursin(r"^-+$", line)
                        break
                    end

                    # Making sure I've not gone off the deep end.
                    prefix = split(line,':')[1]
                    @assert prefix in ["SPECIES", "PROCESS", "PARAM.", "UPDATED", "COLUMNS", "COMMENT"] "Unknown prefix $prefix: line was $line, in CS $name."
                end

                csdatastr = IOBuffer()

                while true
                    line = popfirst!(lines)
                    if occursin(r"^-+$", line)
                        break
                    end

                    write(csdatastr, line * "\n")
                end

                seekstart(csdatastr)
                
                csdata = readdlm(csdatastr)
                energy = csdata[:,1]
                cs = csdata[:,2]

                push!(cs_list, CSData(name=name, threshold=threshold*eV, colltype=colltype, energy=energy*eV, cs=cs*u"m^2", m0=m0))

                if occursin(r"^[\*x]{6,}", lines[1])
                    break
                end
            end

            if occursin(r"^x*$", lines[1])
                popfirst!(lines)
                break
            end
        end
    end

    cs_list
end

mutable struct SwarmData
    name::String
    process::String
    species::String
    ETd::Vector{Float64}
    val::Vector{Float64}
end

function ReadLXCatSwarmSet(filename)
    lines = readlines(filename)
    lines = strip.(lines)
    filter!(x -> !isempty(x), lines)

    data_list= []

    # Go up to the first database definition
    while !occursin(r"^x+$", popfirst!(lines)) end

    while length(lines) > 0
        SkipDatabase(lines)
        
        while true
            SkipHeader(lines)
            while occursin(r"^SPECIES:", lines[1])

                species = "Unassigned"
                process = "Unassigned"
                name = "Unassigned"
                while true
                    line = popfirst!(lines)
                    if occursin(r"^-+$", line)
                        break
                    end

                    # Making sure I've not gone off the deep end.
                    m = match(r"([^:]*):\s*(.*)\s*$", line)
                    @assert m != nothing "Invalid line $line found"
                    prefix,thestr = m.captures

                    @assert prefix in ["SPECIES", "PROCESS", "UPDATED", "COLUMNS", "COMMENT", "PARAM."] "Unknown prefix $prefix: line was $line, in swarm set $name."

                    if prefix == "SPECIES"
                        species = thestr
                    elseif prefix == "PROCESS"
                        process = thestr
                    elseif prefix == "COMMENT"
                        name = thestr
                    end
                end

                datastr = IOBuffer()

                while true
                    line = popfirst!(lines)
                    if occursin(r"^-+$", line)
                        break
                    end

                    write(datastr, line * "\n")
                end

                seekstart(datastr)

                data = readdlm(datastr)
                energy = data[:,1]
                val = data[:,2]

                push!(data_list, SwarmData(name, process, species, energy, val))
            end

            if occursin(r"^x*$", lines[1])
                popfirst!(lines)
                break
            end
        end
    end

    data_list
end

@xport function SetupCollFreqsFromLXCat(params, gas, ptype, lxcat_filename, temperature=0., ionise_kwds=[])
    set = ReadLXCatCSFile(lxcat_filename)

    cfset = []
    for cs in set
        eps = cs.energy * ustrip(u"J",1eV) / params.eps_unit
        csval = cs.cs/params.len_unit^2
        threshold = cs.threshold * ustrip(u"J",1eV) / params.eps_unit

        cf = CreateCollFreq(params, gas, ptype, cs.name, cs.colltype, GENCF_INTERP(eps, csval), threshold=threshold, temperature=temperature ; ionise_kwds...)
        push!(cfset, cf)
    end

    return cfset
end

using LegendrePolynomials

@xport function GenAngleDist(eps_list::Vector{Float64}, sigma_l::Matrix{Float64})
    @assert length(eps_list) == size(sigma_l, 2)

    # Sigma 0 is a normalisation factor (know integrate P_0 and get 2 so need extra factor of 2)
    valid_inds = sigma_l[1,:] .> 0

    sigma_l[:,valid_inds] ./= 2*sigma_l[1:1,valid_inds]

    maxl = size(sigma_l,1) - 1
    legendre_poly = LegendrePoly.(0:maxl+1)

    # Checking whether there will be any negative distribution values.
    for eps_ind in 1:size(sigma_l,2)
        s_l = sigma_l[:,eps_ind]

        s_l[1] == 0 && continue

        tot = Poly(0.)

        if false
            # This is the derivative route
            costh = Poly([0., 1.])
            for l in 1:maxl
                lind = l+1
                tot += s_l[lind] * l*(2l+1) / 4pi * (costh * legendre_poly[lind] - legendre_poly[lind-1])
            end

            valid = true
            for r in roots(tot)
                abs(r) >= 1 && continue

                val = polyval(tot, r)
                if val < 0.
                    valid = false
                    break
                end
            end

            if !valid
                error("A negative angle distribution!")
            end
        else
            # This is the function val
            for l in 0:maxl
                lind = l+1
                tot += s_l[lind] * (2l+1) / 4pi * legendre_poly[lind]
            end

            valid = true
            for r in roots(tot)
                imag(r) != 0 && continue
                abs(r) >= 1 && continue

                # Technically this could be just touching but it is annoying to work that out.
                valid = false
                break
            end

            if !valid
                @show eps_ind
                @show eps_list[eps_ind]
                @show sigma_l[:,eps_ind]
                @show roots(tot)
                error("A negative angle distribution!")
            end
        end
    end

    # Otherwise go with the flow and generate a coeffs lookup table
    legendre_poly_int = polyint.(legendre_poly)
    legendre_poly_int .-= polyval.(legendre_poly_int, -1.)

    newsize = size(sigma_l)
    newsize = (newsize[1]+1, newsize[2])
    coeffs = zeros(newsize)
    for epsind = 1:length(eps_list)
        thispoly = Poly(0.)
        for l in 0:maxl
            lind = l+1
            thispoly += (2l+1) / 4pi * sigma_l[lind, epsind] * legendre_poly_int[lind]
        end

        coeffs[:,epsind] = coeffs(thispoly)
    end

    coeff_interps = [INTERP1D(eps_list, coeffs[ind,:]) for ind in 1:size(coeffs,1)]

    return (eps,R) -> AngleDistPoly(eps, R, coeff_interps)
end

function PolyFromSigmaL(sigma_l)
    maxl = length(sigma_l) - 1

    tot = Poly(0.)
    for l in 0:maxl
        lind = l+1
        tot += sigma_l[lind] * (2l+1) / 4pi * LegendrePoly(l)
    end

    return tot
end

function AngleDistPoly(eps, R, coeff_interps)
    coeffs = Interp.(coeff_interps, eps)

    coeffs[1] -= R

    r = roots(Poly(coeffs))
    r = filter(isreal, r)

    chosen = findfirst(abs.(r) <= 1)
    @assert chosen == findlast(abs.(r) <= 1)

    return chosen
end


function ReadMyCSSet(filename)
    file = open(filename)
    sections = DanUtilsInternal.SeparateBy(readlines(file), x->match(r"^===+$", strip(x)) != nothing)

    cs_sets = map(Iterators.drop(sections, 1)) do lines
        lines = strip.(lines)
        filter!(!isempty, lines)

        keyword_lines,data_lines = DanUtilsInternal.FilterTwo(x->occursin("=", x), lines)

        keyword_pairs = map(keyword_lines) do line
            parts = split(line,"=")
            @argcheck length(parts) == 2
            Pair(lowercase(parts[1]), parts[2])
        end |> Dict

        name = pop!(keyword_pairs, "name")
        coll_type = pop!(keyword_pairs, "coll_type", "unknown")
        threshold = pop!(keyword_pairs, "threshold", "0") |> x->parse(Float64, x)

        if !isempty(keyword_pairs)
            @warn "Unknown keyword pairs left for cross section $name" collect(keys(keyword_pairs))
        end

        data = readdlm(IOBuffer(join(data_lines,"\n")))
        @argcheck size(data,2) == 2

        CSData(name, threshold, coll_type, data[:,1], data[:,2], 0.)
    end

    return cs_sets
end
