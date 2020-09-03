using SwarmMC

for T = [0., 293.]
    params,props = ReadAll("Maxwellmodel:T=$T")

    eps,f = Dists(params,props)

    open("Maxwellmodel_T=$(T)_dists.csv","w") do file
        println(file,"""Distributions for Maxwell model.
                        f0 is normalised to be ∫√(2ϵ/m³) f₀ dϵ = 1.
                        Energy is in eV, and distributions in (mₑ / eV)^(3/2)
                        eps, f0, f1""")
        writecsv(file, [eps f])
    end


    buf = IOBuffer()
    show(IOContext(buf, :limit=>true, :displaysize=>(typemax(Int),typemax(Int))), MIME"text/plain"(), Quants(params,props,include_errors=false,n0_factors=true))

    text = String(take!(buf))
    text = split(text, '\n')[2:end]
    unshift!(text, """Quantites for Maxwell model.
                    Factors of n0 have been included where appropriate.
                    SI units apart from energy in eV.""")
    text = join(text,"\n")

    write("Maxwellmodel_T$(T)_quants.txt", text)
end
