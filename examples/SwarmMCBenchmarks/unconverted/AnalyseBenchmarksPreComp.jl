__precompile__()

module AnalyseBenchmarksPreComp

using SwarmMC

function PrintSection(title)
    println()
    println()
    println("===============")
    println(title)
    println("===============")
end

# This makes me angry
#import Base.transpose
#transpose(x::AbstractString) = x

const numsigfig = 5

function PrintCompare(data, extdata)
    if extdata == nothing
        extdata = Matrix(0,0)
    else
        #extdata = readcsv(extfile)
    end

    allreldiffs = zeros(size(data,1) - 1)
    
    for i = 1:size(data,2)
        item = data[:,i]

        name = item[1]
        vals = item[2:end]

        output = Array{Array{Any,1},1}()
        
        # First print my results
        errs = getfield.(vals, :err)
        vals = Float64.(vals)

        myrelerr = abs.(errs ./ vals)

        #println(name, ", ", join(string.(vals), ", "))
        push!(output, [name ; signif.(vals, numsigfig)])

        # Now print all of the loaded stuff
        diff = ones(length(vals)) * Inf
        for extline = 1:size(extdata,1)
            if extdata[extline, 1] == name
                #println(join(extdata[extline, 2:end], ", "))
                push!(output, extdata[extline, 2:end])

                #thisvals = tryparse.(Float64,extdata[extline, 2:end])
                try
                    thisvals = Vector{Float64}(extdata[extline, 3:end])
                    thisdiff = abs.(vals - thisvals)
                    diff = min.(diff, thisdiff)
                end
            end
        end
        reldiff = abs.(diff ./ vals) * 100
        #println("MinDiff(%), ", join(string.(signif.(reldiff, numsigfig)), ", "))
        #println("MyRelErr(%), ", join(string.(signif.(myrelerr*100, numsigfig)), ", "))
        push!(output, ["MinDiff(%)" ; round.(reldiff, sigdigits=3)])
        push!(output, ["MyRelErr(%)" ; round.(myrelerr*100, sigdigits=3)])

        #display(output)
        temp = hcat(output...)
        Base.print_matrix(STDOUT, permutedims(temp, [2,1]))
        #Base.showarray(STDOUT, permutedims(temp, [2,1]), false, header=false)
        println()
        println()

        reldiff[isinf.(reldiff) .| isnan.(reldiff)] = 0.

        allreldiffs = max.(allreldiffs, reldiff)
    end

    allreldiffs /= 100.
    println("LargestRelDiffs(not %), ", join(string.(signif.(allreldiffs, 2)), ", "))
end
    
DataFuncs = Dict( "ReidRamp" => (quants) -> GetAll(quants..., restrict=r"ReidRamp_", n0_factors=true, quiet=true),
                  "ReidRampB" => (quants) -> GetAll(quants..., restrict=r"ReidRampB_", n0_factors=true, quiet=true),
                  "LucasSaelee" => (quants) -> GetAll(quants..., restrict=r"LucasSaelee_", n0_factors=true, quiet=true),
                  "LucasSaeleeB" => (quants) -> GetAll(quants..., restrict=r"LucasSaeleeB_", n0_factors=true, quiet=true),
                  "NessRobsonSharing" => (quants) -> GetAll(quants..., restrict=r"NRS", n0_factors=true, quiet=true),
                  "PercusYevick" => (quants) -> GetAll(quants..., restrict=r"PercusYevick", n0_factors=true, quiet=true),
                  "ReidAniso" => (quants) -> GetAll(quants..., restrict=r"ReidAniso", n0_factors=true, quiet=true),
                  # "NessRobsonLoss" => (quants) -> GetAll(quants..., restrict=r"NRL", n0_factors=true, quiet=true),
                  "Maxwell" => (quants) -> GetAll(quants..., restrict=r"Maxwellmodel", n0_factors=true, quiet=true),
                  "Hardsphere" => (quants) -> GetAll(quants..., restrict=r"Hardspheremodel", n0_factors=true, quiet=true),
                  "MassRatio" => (quants) -> GetAll(quants..., restrict=r"MassRatio", n0_factors=true, quiet=true),
                  "BoyleSpacePaper" => (quants) -> GetAll(quants..., restrict=r"BoyleSpacePaper", n0_factors=true, quiet=true),
                  "MagneticField" => (quants) -> GetAll(quants..., restrict=r"MagneticField", n0_factors=true, quiet=true),
                 )

function NRLData(quants)
    data = GetAll(quants..., restrict=r"NRL.*SUBSTEPS.*DOUBLE", n0_factors=true, quiet=true)

    valid_inds = filter(1:size(data,2)) do col_ind
        !contains(data[1,col_ind], "T=")
    end

    data = data[:,valid_inds]

    data[1,:] = map(data[1,:]) do name
        temp = split(name, ":")
        temp = temp[1:3]
        temp = join(temp, ":")
    end

    data
end
DataFuncs["NessRobsonLoss"] = NRLData

function DoOne(x)
    local quants, extdata
    open("extdata_$(x).csv") do file
        line = readline(file)
        quants = split(line, ",") .|> strip
        @assert all(x -> startswith(x,":"), quants)
        quants = [Symbol(z[2:end]) for z in quants]

        extdata = readcsv(file)
    end

    PrintSection(x)
    local data
    try
        data = DataFuncs[x](quants)
    catch exc
        print_with_color(:red, "Exception in '$x': $exc")
        #return
        rethrow()
    end
    #PrintCompare(data, "extdata_$(x).csv")
    PrintCompare(data, extdata)
end


function Main()
    if length(ARGS) > 0
        for arg in ARGS
            @assert arg in keys(DataFuncs)
        end
        for arg in ARGS
            DoOne(arg)
        end
    else
        for name in keys(DataFuncs)
            DoOne(name)
        end
    end
end
    
end
