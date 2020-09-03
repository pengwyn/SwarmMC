using SwarmMC

const numsigfig = 5

dataexists = false
try
    _ = data
    if data !== false
        dataexists = true
    end
catch exc
end

quants = ["p", "a", "ans", "gns", "N", "red", "split", :countR, :eps, :bulk_W, :bulk_DT, :bulk_DL, :flux_W, :flux_DT, :flux_DL, :walltime, :numinit, :meas_eff_alt]
if !dataexists
    data = GetAll(quants... ; restrict=r"NRL", n0_factors=true)

    valid_inds = filter(1:size(data,2)) do col_ind
        !contains(data[1,col_ind], "T=")
    end

    data = data[2:end,valid_inds]
end



# # Loading Greg data
# gregdata = readcsv("extdata_NessRobsonLoss.csv", skipstart=1)
# gregdata = gregdata[gregdata[:,1] .!= "", :]

# gregdata = mapslices(gregdata, 2) do set
#     name = set[1]
#     m = match(r"NRL_p([-.\d]*)_a([-e.\d]*)", name)
#     p = parse(Float64, m.captures[1])
#     a = parse(Float64, m.captures[2])

#     out = [p;a;set[3];set[3];set[4:end]]
# end



paset = Set([data[:,i][[1,2]] for i in 1:size(data,2)])
    
compare_params = ["SUBSTEPS", "DOUBLE", 1000, 1.0, false]

function rows(A)
    [A[i,:] for i = 1:size(A, 1)]
end

for pa in sort(collect(paset), by=x->tuple(x...))
    println("Doing pa = $pa.")

    #thisdata = [data[3:end,i] for i in 1:size(data,2) if data[1:2,i] == pa]
    thisdata = [data[:,i] for i in 1:size(data,2) if data[1:2,i] == pa]

    thisdata = matrix(thisdata)
    #thisdata[6:end,:] = Float64.(thisdata[6:end,:])
    thisdata[8:end,:] = Float64.(thisdata[8:end,:])

    # gregind = findfirst(x -> x[1:2] == pa, rows(gregdata))
    # @show gregind
    # correct = gregdata[gregind, 3:end]
    # correct = [correct ; 1 ; 1 ; 1]

    # errs = (thisdata[6:end, :] .- correct) ./ correct
    # errs = signif.(errs, 2)
    # errs = [thisdata[1:5, :] ; errs]

    
    if true
        thisdata = permutedims(thisdata, [2,1])
        # errs = permutedims(errs, [2,1])

        #thisdata = ["ANS" "GNS" "RED" "SPLITFAKE" "N" string.(hcat(quants...)) ; thisdata ; ["Greg" "Greg" "Greg" "Greg" "Greg" correct']]
        #thisdata = [string.(hcat(quants[3:end]...)) ; thisdata]
        thisdata = [string.(hcat(quants...)) ; thisdata]
        # errs = ["ANS" "GNS" "RED" "SPLITFAKE" "N" string.(hcat(quants...)) ; errs]

        io = IOContext(STDOUT, limit=false)
        Base.showarray(io, thisdata, false)
        println()
        # display(errs)
        #Base.showarray(STDOUT,errs)
    end
end
    
