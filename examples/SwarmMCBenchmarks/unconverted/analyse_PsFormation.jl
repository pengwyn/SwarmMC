
ENV["MPLBACKEND"] = "Agg"
using PyPlot

# islesstest(x::Real, y::Real) = isless(x,y)
# islesstest(x::T, y::T) where {T<:AbstractArray} = all(islesstest.(x,y))
# islesstest(x::Symbol, y::Real) = false
# islesstest(x::Real, y::Symbol) = true
import Base: isless
isless(x::Symbol, y::Real) = false
isless(x::Real, y::Symbol) = true

using JLD2, FileIO, DataFrames
@load "analyse_PsFormation.jld2"

df = DataFrame(permutedims(data), labels)
# # An issue I need to fix
# df[:ion_lam] = Float64.(df[:ion_lam])

# Plotting opposite
df[:PsFrac] = 1 .- df[:PsFrac]

coarse_temperature = [0.]
coarse_init = [1e3, 1e4, 1e5]
coarse_spread = [:gauss_twiceth, :gauss_tenth]
coarse_Q = [0.5, 0.8, 0.9, 0.99, 1.0]

meas_params = [:PsFrac,:NumEl,:NumExc,:NumIon,:AvgPsEnergy,:AvgThermalEnergy]

# overall_params = [:spread => :gauss_tenth,
overall_params = [:spread => :gauss_twiceth,
                  :init => 1e4]
# Set A

setA_params = [:ion_lam => 30.,
               :ion_A => 0.3,
               :exc_A => 0.,
               :exc_thresh => 0.,
               :exc_lam => 1.,
               :T => 0.]

setB_params = [:ion_lam => 30.,
               :ion_A => 0.3,
               :exc_A => 0.5,
               :exc_thresh => 10.,
               :exc_lam => 12.,
               :T => 0.]

setC_params = [:ion_lam => 30.,
               :ion_A => 0.3,
               :exc_A => 1.0,
               :exc_thresh => 2.,
               :exc_lam => 10.]

for (setname,set_params,plotvars) in [
                                      ("SetA", setA_params, [:ion_A,:ion_lam,:Q]),
                                      ("SetB", setB_params, [:exc_A,:exc_lam,:exc_thresh]),
                                      ("SetC", setC_params, [:exc_A,:exc_lam,:T])
    ],
    plotvar in plotvars

    thisdf = df
    for (param,val) in [set_params ; overall_params]
        param == plotvar && continue
        # Temporary issue
        # if plotvar == :ion_lam && param == :ion_A
        #     val = 1.0
        # end
        # if plotvar ∈ [:exc_thresh,:exc_lam] && param == :exc_A
        #     val = 1.0
        # end
        # if plotvar == :exc_lam && param == :exc_thresh
        #     val = 15.0
        # end

        thisdf = thisdf[thisdf[param] .== val, names(thisdf) .!= param]
    end

    others = names(thisdf)
    others = filter!(x->x ∉ [plotvar ; meas_params], others)

    if isempty(others)
        groups = [thisdf]
    else
        groups = groupby(thisdf, others)
        groups = collect(groups)
        # Sort by the first row in each gruop
        sort!(groups, by=x->vec(convert(Array, x[1,others])))
    end

    @info "Doing $(length(groups)) sets" plotvar
    figure(figsize=(6,4)); axfrac = PyPlot.axes()
    figure(figsize=(6,4)); axnum = PyPlot.axes()
    figure(figsize=(6,4)); axenergy = PyPlot.axes()
    for set in groups
        size(set,1) > 3 || continue

        set = set[sortperm(set[plotvar]), :]

        t = join(["$z=$(set[1,z])" for z in others], ", ")

        # figure()
        axfrac[:plot](set[plotvar], set[:PsFrac], marker="o", label=t)

        l, = axnum[:plot](set[plotvar], set[:NumEl], marker="o", label=t)
        axnum[:plot](set[plotvar], set[:NumIon], linestyle="--", marker="v", color=l[:_color])

        axenergy[:plot](set[plotvar], set[:AvgPsEnergy], marker="o", label=t)
        l, = axenergy[:plot](set[plotvar], set[:AvgThermalEnergy], linestyle="--", marker="v", color=l[:_color])

        #display(set)
    end
    for ax in [axfrac,axnum,axenergy]
        ax[:legend]()
        ax[:set_xlabel](plotvar)
        ax[:grid]()
    end

    # axfrac[:set_ylabel]("Ps formation %")
    axfrac[:set_ylabel]("Thermal positron %")
    axfrac[:figure][:savefig]("Psform_$(setname)_$(plotvar)_PsFrac.png")

    axnum[:set_title]("Number of elastic (solid) vs ionising (dashed) collisions")
    axnum[:set_yscale]("log")
    axnum[:figure][:savefig]("Psform_$(setname)_$(plotvar)_Num.png")

    axenergy[:set_title]("Average energy when forming Ps (solid) or 'thermalising' (dashed)")
    axenergy[:set_ylabel]("Energy (eV)")
    axenergy[:figure][:savefig]("Psform_$(setname)_$(plotvar)_AvgEnergy.png")
end


# Figure plots

plot_size = (3,4)

# * Set A
set_params = [overall_params ;
              setA_params ;
              :Q => 1.]


####
# ** Ion A

thisdf = FilterDF(df, set_params, :ion_A)

figure(figsize=plot_size)
plot(thisdf[:ion_A], thisdf[:PsFrac], marker="o")

xlabel(raw"$A_\mathrm{ion}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetA_ionA.png")


####
# ** Ion Lam

thisdf = FilterDF(df, set_params, :ion_lam)

figure(figsize=plot_size)
plot(thisdf[:ion_lam], thisdf[:PsFrac], marker="o")

xlabel(raw"$\lambda_\mathrm{ion}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetA_ionlam.png")

####
# ** Q

set_params = [overall_params ;
              setA_params]

# thisdf = FilterDF(df, set_params, :Q)
# This is annoying
thisdf = df[isa.(df[:ion_lam], Float64), :]
thisdf = FilterDF(thisdf, set_params, :Q)

H_val = thisdf[thisdf[:Q] .== :hydrogen_dist, :PsFrac]

# Get rid of the non-numeric Q
thisdf = thisdf[isa.(thisdf[:Q],Real), :]
# But convert to non-numeric so that we can see all the fine detail
#thisdf[:Q] = string.(thisdf[:Q])

# Restrict Q range
thisdf = thisdf[thisdf[:Q] .>= 0.5, :]

# This version used to plot as a histogram type thing
# thisdf = thisdf[thisdf[:Q] .!= :hydrogen_mean, :]
# thisdf[thisdf[:Q] .== :hydrogen_dist,:Q] = :H

figure(figsize=plot_size)
plot(thisdf[:Q], thisdf[:PsFrac], marker="o")#, linestyle="")
plot(0.55, H_val, "x", markersize=10, markeredgewidth=5)

xlabel(raw"$\tilde{Q}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetA_Q.png")

# * Set B
set_params = [overall_params ;
              setB_params ;
              :Q => 1.]
hyd_set_params = [overall_params ;
              setB_params ;
              :Q => :hydrogen_dist]

hyd = "r--"

####
# ** Exc A

figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :exc_A)
plot(thisdf[:exc_A], thisdf[:PsFrac], marker="o")

# Hydrogen
thisdf = FilterDF(df, hyd_set_params, :exc_A)
plot(thisdf[:exc_A], thisdf[:PsFrac], hyd, marker="o")

xlabel(raw"$A_\mathrm{exc}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetB_excA.png")


####
# ** Exc threshold

figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :exc_thresh)
plot(thisdf[:exc_thresh], thisdf[:PsFrac], marker="o")

thisdf = FilterDF(df, hyd_set_params, :exc_thresh)
plot(thisdf[:exc_thresh], thisdf[:PsFrac], hyd, marker="o")

xlabel(raw"$\epsilon_\mathrm{exc}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetB_excthreshold.png")

####
# ** Exc Lam

figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :exc_lam)
plot(thisdf[:exc_lam], thisdf[:PsFrac], marker="o")

thisdf = FilterDF(df, hyd_set_params, :exc_lam)
plot(thisdf[:exc_lam], thisdf[:PsFrac], hyd, marker="o")

xlabel(raw"$\lambda_\mathrm{exc}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetB_exclam.png")


# * Set C
set_params = [overall_params ;
              setC_params ;
              :Q => 1.]
hyd_set_params = [overall_params ;
              setC_params ;
              :Q => :hydrogen_dist]

hyd = "r--"

####
# ** Exc A

figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :exc_A)
plot(thisdf[:exc_A], thisdf[:PsFrac], marker="o")

# Hydrogen
thisdf = FilterDF(df, hyd_set_params, :exc_A)
plot(thisdf[:exc_A], thisdf[:PsFrac], hyd, marker="o")

xlabel(raw"$A_\mathrm{exc}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetC_excA.png")


####
# ** Exc threshold

figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :exc_thresh)
plot(thisdf[:exc_thresh], thisdf[:PsFrac], marker="o")

thisdf = FilterDF(df, hyd_set_params, :exc_thresh)
plot(thisdf[:exc_thresh], thisdf[:PsFrac], hyd, marker="o")

xlabel(raw"$\epsilon_\mathrm{exc}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetC_excthreshold.png")

####
# ** Exc Lam

figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :exc_lam)
plot(thisdf[:exc_lam], thisdf[:PsFrac], marker="o")

thisdf = FilterDF(df, hyd_set_params, :exc_lam)
plot(thisdf[:exc_lam], thisdf[:PsFrac], hyd, marker="o")

xlabel(raw"$\lambda_\mathrm{exc}$")
ylabel("Thermal positron %")
grid()

tight_layout()
savefig("Psform_paper_figs_SetC_exclam.png")

# * Set A Times
set_params = [overall_params ;
              setA_params ;
              :Q => 1.]
hyd_set_params = [overall_params ;
              setA_params ;
              :Q => :hydrogen_dist]


####
# ** Ion A


figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :ion_A)
plot(thisdf[:ion_A], thisdf[:TimeToThermal], marker="o")

# Hydrogen
thisdf = FilterDF(df, hyd_set_params, :ion_A)
plot(thisdf[:ion_A], thisdf[:TimeToThermal], hyd, marker="o")

xlabel(raw"$A_\mathrm{ion}$")
ylabel("Thermalisation time")
grid()

tight_layout()
savefig("Psform_paper_figs_SetAtime_ionA.png")

####
# ** Ion Lam


figure(figsize=plot_size)

thisdf = FilterDF(df, set_params, :ion_lam)
plot(thisdf[:ion_lam], thisdf[:TimeToThermal], marker="o")

# Hydrogen
thisdf = FilterDF(df, hyd_set_params, :ion_lam)
plot(thisdf[:ion_lam], thisdf[:TimeToThermal], hyd, marker="o")


xlabel(raw"$\lambda_\mathrm{ion}$")
ylabel("Thermalisation time")
grid()

tight_layout()
savefig("Psform_paper_figs_SetAtime_ionlam.png")

####
# ** Q

set_params = [overall_params ;
              setA_params]

# thisdf = FilterDF(df, set_params, :Q)
# This is annoying
thisdf = df[isa.(df[:ion_lam], Float64), :]
thisdf = FilterDF(thisdf, set_params, :Q)

H_val = thisdf[thisdf[:Q] .== :hydrogen_dist, :TimeToThermal]

# Get rid of the non-numeric Q
thisdf = thisdf[isa.(thisdf[:Q],Real), :]
# But convert to non-numeric so that we can see all the fine detail
#thisdf[:Q] = string.(thisdf[:Q])

# Restrict Q range
thisdf = thisdf[thisdf[:Q] .>= 0.5, :]

# This version used to plot as a histogram type thing
# thisdf = thisdf[thisdf[:Q] .!= :hydrogen_mean, :]
# thisdf[thisdf[:Q] .== :hydrogen_dist,:Q] = :H

figure(figsize=plot_size)
plot(thisdf[:Q], thisdf[:TimeToThermal], marker="o")#, linestyle="")
plot(0.85, H_val, "x", markersize=10, markeredgewidth=5)

xlabel(raw"$\tilde{Q}$")
ylabel("Thermalisation time")
grid()

tight_layout()
savefig("Psform_paper_figs_SetAtime_Q.png")



# * Set A distributions

# TODO: Need to save off the particular distributions in grabdata.


# * The end

# Copy everything to the paper folder
using Glob
for filename in glob("Psform_paper_*.png")
    cp(filename, expanduser("~/Dropbox/Physics/Papers/PsFormation/paper_figs/$filename"), force=true)
end
