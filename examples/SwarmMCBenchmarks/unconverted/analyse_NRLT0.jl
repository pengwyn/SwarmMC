using SwarmMC

GetT(params, props, prefix, quants) = contains(prefix, "T=") ? NameElements(prefix, "T") : 293.

#ShowAll("p","a","T", :countR, :eps, :bulk_W, :bulk_DT, :bulk_DL, :flux_W, :flux_DT, :flux_DL, :numinit,
ShowAll("p","a", GetT, :countR, :eps, :bulk_W, :bulk_DT, :bulk_DL, :flux_W, :flux_DT, :flux_DL, :numinit,
		#restrict=r"NRL.*T=0.",
		restrict=r"NRL.*ans=SUBSTEPS.*gns=REGEN_ALL.*red=0.5.*N=1000",
		quiet=true, include_prefix=false, include_sortby=false, n0_factors=true,
		sortby=(x,y)->(NameElements(x,"p"),NameElements(x,"a")))



asdf


#prefixes = filter((x,y)->startswith(x, "NRL") && ismatch(r"T=0.",x), PrefixSets()) |> keys |> collect |> sort
prefixes = PrefixSets(r"^NRL.*T=0\.") |> keys |> collect |> sort

using PyPlot
figure()
for prefix in prefixes
    params,props = ReadAll(prefix)
    plot(TGrid(params), props[:eps])
end
