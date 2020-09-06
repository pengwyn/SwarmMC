# Does the comparison to the crossed E,B results

using SwarmMC
using DataFrames

function MaybePlot(df, args...)
    # Only plot those ones which have all the crossing angles
    if(all(∈(df[:Bθ]), [0,30,60,90]))
        Plot(df, args...)
    end
end

function Plot(df, p, q)
    F=unique(df[:F])[]
    B=unique(df[:BHx])[]

    cols = Dict(100 => :red,
                200 => :blue,
                1000 => :green)
    lines = Dict(0.0 => :solid,
                 0.5 => :dash,
                 1.0 => (:solid,2))
    
    # func = (q == :countR) ? semilogy! : plot!
    func = plot!
    func(p,
          df.Bθ, ustrip.(Value.(df[q])),
          yerror=ustrip.(StdErr.(df[q])),
          label="F=$F, B=$B Hx",
          color=cols[B],
          line=lines[F])

end

df = cd("runs_LucasSaelee") do
    GetAll("F", "BHx", "Bθ", "ion_weight", :energy, :countR, :bulk_W, :bulk_Wx, :bulk_Wy, :flux_W)
end
groups = DataFrames.groupby(df, [:F,:BHx,:ion_weight])

function PlotQuant(q)
    p = plot(title=q)
    foreach(g -> MaybePlot(g, p, q), groups)
    p
end

plot(PlotQuant.([:energy,:bulk_W,:bulk_Wx,:bulk_Wy,:countR])...)
