using DanUtils
using SwarmMC, StaticArrays, Dierckx, BisectInterp
using DataFrames, CSV


df = cd("runs_SuperModel") do
    # Got to get all and then filter on those where the temperatures are equal.

    df = GetAll("elT", "inelT", "inelrecoil", :ETd, :eps, :flux_W, :bulk_W, sortby=["elT", "inelrecoil", :ETd], as_df=true, quants_kwds=[:include_errors => false])

    df[!,:ETd] .*= 2
    Update!(x -> round(x,sigdigits=2), df[!,:ETd])
    df = df[df.elT .== df.inelT,:]
    df = df[!,setdiff(names(df),[:inelT])]

    temp = names(df)
    temp[1] = :T
    names!(df, temp)

    CSV.write("supermodel_reruns.csv", df)

    # I'm so lazy
    df = GetAll("elT", "inelT", "inelrecoil", :ETd, :eps, :flux_W, :bulk_W, sortby=["elT", "inelrecoil", :ETd], as_df=true)

    df[!,:ETd] .*= 2
    Update!(x -> round(x,sigdigits=2), df[!,:ETd])
    df = df[df.elT .== df.inelT,:]
    df = df[!,setdiff(names(df),[:inelT])]

    temp = names(df)
    temp[1] = :T
    names!(df, temp)

    CSV.write("supermodel_reruns_witherrors.csv", df)

    return df
end
