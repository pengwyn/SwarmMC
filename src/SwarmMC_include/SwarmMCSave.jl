
@xport function MakeSaveName(prefix ; separate_dir=true, kwds...)
    items = map(collect(kwds)) do pair
        name,val = pair
        suffix = ""
        if val isa Unitful.AbstractQuantity
            # suffix = replace(string(unit(val)), " " => "")
            val = ustrip(val)
        end
        
        if val isa Float64
            val = round(val, sigdigits=6)
        elseif val isa AbstractString
            val = replace(val, ":" => "_")
            val = replace(val, "=" => "")
            val = replace(val, "/" => "_")
            val = replace(val, "__" => "_")
        end
        # Singleton types
        if isdefined(typeof(val), :instance)
            val = typeof(val)
        end
        item = string(name) * "=" * string(val) * suffix
    end

    items = [prefix ; items]

    thestr = join(items, ":")

    if separate_dir
        thestr = "runs_" * prefix * "/" * thestr
    end

    return thestr
end

function SavePrefix(p::PARAMS)
    return p.save_name
end

@xport function Save(in_params::PARAMS, props::PROPS_OUT, prefix::String=SavePrefix(in_params) ; mmap::Bool=false)
    in_params = Finalise(in_params)

    pdict = ConvertToDict(in_params)

    path = dirname(prefix)
    if path != "" && !ispath(path)
        mkpath(path)
    end

    filename = DanUtils.ClaimNextFilename(prefix)

    tempfilename = filename * ".TMPSAVE"
    if mmap
        @msgwrap "Saving to $filename" JLD2.@save tempfilename props pdict
    else
        # FIXME: Change to the new syntax for IOStream
        @msgwrap "Saving (IOStream) to $filename" begin
            jldopen(tempfilename, true, true, true, IOStream) do file
                # write(file, "props_save", props_save)
                write(file, "props", props)
                write(file, "pdict", pdict)
            end
        end
    end

    #mv(tempfilename, filename, remove_destination=true)
    # Going copy instead, due to mmap/NFS/JLD2 problems
    cp(tempfilename, filename, force=true)
    rm(tempfilename)
end

using Glob
@xport function PrefixSets(restrict=r"", quiet=false)
    filelist = glob("*__*.jld2")
    sort!(filelist)

    prefixset = Set{String}()
    for filename in filelist
        ind = findlast("__", filename)
        push!(prefixset, filename[1:ind.start-1])
    end

    filter!(x -> occursin(restrict, x), prefixset)

    thedict = Dict{String,Vector{String}}()

    for prefix in prefixset
        prefixfilelist = filter(filelist) do filename
            startswith(filename, prefix * "__")
        end

        prefixfilelist = prefixfilelist

        thedict[prefix] = prefixfilelist
    end

    thedict
end

@xport function CollectUp(args... ; recurse=true, kwds...)
    if recurse
        for (dirname,) in walkdir(".")
			basename(dirname) == "bad" && continue
            @show dirname
            cd(dirname) do
                CollectUpOneDir(args... ; kwds...)
            end
        end
    else
        CollectUpOneDir(args... ; kwds...)
    end
end
@xport function CollectUpOneDir(quiet::Union{Bool,String}=false, force::Bool=false, BeforePrefix::Union{Nothing,Function}=nothing ; move_bad_out=false, summarise=true, mmap=true)
    prefix_dict = PrefixSets(r"", quiet)

    num_sets_collected = 0
    num_files_collected = 0
    num_sets_ignored = 0
    num_errors = 0
    
    for prefix in sort(collect(keys(prefix_dict)))
        filelist = prefix_dict[prefix]
        if ! (BeforePrefix isa Nothing)
            BeforePrefix(prefix)
        end

        if length(filelist) == 0
            @error "Shouldn't get here anymore"
            quiet == false && println("No files finished saving.")
            num_errors += 1
            continue
        end
        
        if !force && length(filelist) == 1
            quiet == false && println("Not collecting $prefix: only one file.")
            if filelist[] != DanUtils.DefaultClaimFileGen(0, prefix)
                quiet == false && println("However, renaming.")

                filename = DanUtils.ClaimNextFilename(prefix)
                mv(filelist[1], filename, force=true)
            end

            num_sets_ignored += 1
            continue
        end

        !(quiet == "extra") && println("Collecting up $prefix ($(length(filelist)) files)")
        num_sets_collected += 1

        local props, params
        try
            params,props = ReadAll(prefix, print_dots=true, move_bad_out=move_bad_out, keep_on_error=false)
            num_files_collected += length(filelist)
        catch exc
            @error "Unable to collect up files. Skipping set (nothing modified)." prefix stacktrace(catch_backtrace())
            continue
        end

        # Not deleting all the files just yet
        for filename in filelist
            # Files could have been moved out if they were bad
            isfile(filename) && mv(filename, filename * ".stored")
        end

        try
            Save(params, props, prefix, mmap=mmap)
        catch exc
            @error "Error when saving. Skipping set (nothing modified)." prefix stacktrace(catch_backtrace())
            continue
        end

        # Now delete!
        for filename in filelist
            # Files could have been moved out if they were bad
            del_filename = filename * ".stored"
            isfile(del_filename) && rm(del_filename)
        end
    end

    if summarise
        println("Collected $num_sets_collected sets with $num_files_collected files ignoring $num_sets_ignored sets. There were $num_errors errors.")
    end
end

struct NoValidFiles <: Exception
    prefix::String
end

ReadAll(params::PARAMS; kwds...) = ReadAll(SavePrefix(params); kwds...)
@xport function ReadAll(prefix::String; quiet=false, print_dots=false, move_bad_out=false, keep_on_error=true)
    filelist = glob(prefix * "__*.jld2")

    quiet == false && println("Found $(length(filelist)) files.")

    if length(filelist) == 0
        throw(NoValidFiles(prefix))
    end

    bad_dir = "bad"

    # Unfortunately can't do this with possible errors.
    #props = mapreduce(x->JLD2.load(x, "props"), +, filelist)
    props = nothing
    params = nothing
    i = 1
    for filename in sort(filelist)
        if quiet == false
            if print_dots
                i%10 == 0 ? print(i) : print(".")
                i += 1
            else
                println("Loading from $filename")
            end
        end

        if filesize(filename) == 0
            println("File $filename saving wasn't completed. Skipping.")
            if move_bad_out
                println("Moving $filename to '$(bad_dir)'.")
                ispath(bad_dir) || mkdir(bad_dir)
                mv(filename, joinpath(bad_dir, filename), force=true)
            end
            continue
        end

        try
            if params == nothing
                pdict = load(filename, "pdict")
                params = ParseFromDict(pdict, keep_on_error)
            end

            thisprops = load(filename, "props")

            if props == nothing
                props = thisprops
            else
                try
                    props += thisprops
                catch exc
                    @error "Can't combine props together." filename
                    rethrow()
                end
            end
        catch exc
            if isa(exc, JLD2.InvalidDataException) || isa(exc, JLD2.UnsupportedVersionException) || isa(exc, FileIO.UnknownFormat) || isa(exc, MethodError) || isa(exc, EOFError)
                if move_bad_out
                    println("Got exception $exc - moving $filename to '$(bad_dir)'.")
                    if !ispath(bad_dir)
                        mkdir(bad_dir)
                    else
                        @assert isdir(bad_dir)
                    end

                    mv(filename, bad_dir * "/" * filename, force=true)
                else
                    showerror(stdout, exc)
                    @error "Ignoring exception - skipping file $filename."
                end
                continue
            end
            @error "Problem in file $filename"
            rethrow()
        end
    end

    if params == nothing
        throw(NoValidFiles(prefix))
    end

    params, props
end



using ObjectSaving
import ObjectSaving: ShouldConvertToDict, ConvertToDict

# Need to have the dimensions accessed from ObjectSaving in here.
using Unitful: 𝐓, 𝐋, 𝐌


############################################
# * SwarmMC conversions
#------------------------------------------

ConvertToDict(p::PARAMS) = ConvertToDict(p, [:tot_cf_matrix])

ShouldConvertToDict(::PARAMS) = true

ShouldConvertToDict(::GAS) = true
# ShouldConvertToDict(::REGION) = true
ShouldConvertToDict(::PARTTYPE) = true
ShouldConvertToDict(::COLLFREQ) = true
ShouldConvertToDict(::MEAS_BIN) = true

ShouldConvertToDict(::PARTICLE_INIT_STYLE) = true
ShouldConvertToDict(::PARTICLE_INIT_VEL_STYLE) = true
ShouldConvertToDict(::PARTICLE_INIT_SPATIAL_STYLE) = true
ShouldConvertToDict(::PARTICLE_INIT_TIME_STYLE) = true

ShouldConvertToDict(::ISS_FUNC) = true

function CreateObjectFromDict(::Type{PARAMS}, dict::Dict)
    # Throw away any keys that aren't in PARAMS.
    for key in keys(dict)
        if key ∉ fieldnames(PARAMS)
            @warn "Removing key $key for creation of PARAMS."
            delete!(dict, key)
        end
    end

    p = PARAMS(;dict...)
    return p
end


##############################
# * NameElements
#----------------------------

using DataStructures
@xport function NameElements(str::AbstractString ; doparse=true)
    if occursin(":", str)
        splitchar = ':'
    else
        splitchar = '_'
    end
    
    elements = split(str,splitchar)

    name = elements[1]
    others = map(elements[2:end]) do item
        if occursin("=", item)
            parts = split(item, "=")
            @assert length(parts) == 2
            first = parts[1]
            second = parts[2]
        else
            m = match(r"^(\pL*)(.*)", item)
            if isempty(m.captures[1])
                first = second = m.captures[2]
            elseif isempty(m.captures[2])
                first = second = m.captures[1]
            else
                first = m.captures[1]
                second = m.captures[2]
            end
        end

        if doparse
            local val
            try
                val = Meta.parse(second)
                if isa(val, Symbol)
                    val = string(val)
                end
            catch
                val = second
            end
        else
            val = second
        end

        first => val
    end

    OrderedDict("name" => name, others...)
end
function NameElements(str::AbstractString, index::AbstractString, force_T::Type=Any)
    val = NameElements(str, doparse=false)[index]
    if force_T == Any
        try
            Meta.parse(val)
        catch
            # Redo as a string
            string(val)
        end
    elseif force_T <: AbstractString
        force_T(val)
    else
        parse(force_T, val)
    end
end
function NameElements(str::AbstractString, indices::Vector{T} where T<:AbstractString, force_T::Type=Any)
    out = map(x -> NameElements(str, x, force_T), indices)
    tuple(out...)
end

