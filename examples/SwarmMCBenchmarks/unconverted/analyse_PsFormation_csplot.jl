using SwarmMC
using PsFormation

params = SetupParams(init_eps=1e4,
                     eps_spread=:gauss_tenth,
                     ion_A=0.3,
                     ion_λ=30.,
                     Q=1.,
                     exc_A=0.5,
                     exc_threshold=10.,
                     exc_λ=12.)


params = Finalise(params)

# Doing the plotting manually to make it look nice
# using Plots
using PyPlot


en = EpsGrid(params)

function Filt(x,y)
    # inds = (x .> 0) .& (y .> 1e-10)
    inds = (x .> 0)
    x[inds],max.(y[inds], 1e-10)
end

# fac = params.density[1,1] * sqrt.(2*en/1.)
fac = sqrt.(2*en/1.)

PyPlot.plt[:rc]("font", size=18, family="serif")
PyPlot.plt[:rc]("mathtext", fontset="cm")
PyPlot.plt[:rc]("lines", linewidth=2.5)
PyPlot.plt[:rc]("axes", labelsize=18)

figure(figsize=(8,5))
loglog(Filt(en, fac.\params.cf_list[1][1].func.(en))..., label="Elastic")
loglog(Filt(en, fac.\params.cf_list[1][2].func.(en))..., label="Ps")
loglog(Filt(en, fac.\params.cf_list[1][3].func.(en))..., label="Ionisation")
loglog(Filt(en, fac.\params.cf_list[1][4].func.(en))..., label="Excitation")
xlim(1e0, 1e4)
ylim(1e-3,1.1)
xlabel(raw"Energy (eV)")
ylabel(raw"Cross section ($\AA^2$)")
grid()
tight_layout()

savefig("paper_fig_cs.eps")
savefig("paper_fig_cs.pdf")
