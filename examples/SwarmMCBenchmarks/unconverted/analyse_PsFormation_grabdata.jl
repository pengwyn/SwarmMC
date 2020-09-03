
using SwarmMC, StaticArrays, PsFormation, JLD2
# using DataFrames


# Also save off a couple of distributions
try
    params,props = cd("runs_PsFormation") do
        ReadAll("PsFormation:T=0.0:recoil=false:Q=1.0:ion_A=0.3:ion_lam=30.0:exc_A=0.0:exc_thresh=0.0:exc_lam=1.0:init=10000.0:spread=gauss_tenth")
    end
    dist = props[:cf_5_num, (:eps,:cum)]
    dist ./= DEps(params)
    energy = EpsGrid(params)

    @save "analyse_PsFormation_setA_distQ1.jld2" energy dist

    params,props = cd("runs_PsFormation") do 
        ReadAll("PsFormation:T=0.0:recoil=false:Q=hydrogen_dist:ion_A=0.3:ion_lam=30.0:exc_A=0.0:exc_thresh=0.0:exc_lam=1.0:init=10000.0:spread=gauss_tenth")
    end
    dist = props[:cf_5_num, (:eps,:cum)]
    dist ./= DEps(params)
    energy = EpsGrid(params)

    @save "analyse_PsFormation_setA_disthyd.jld2" energy dist
catch exc
    @error "There was a problem trying to save the distributions - maybe the files don't exist?" exc
end




# Original stuff

function PsFrac(params, props, args...)
    ps = props[:cf_2_num, :cum]
    an = props[:cf_5_num, :cum]

    frac = ps / (ps + an)
end

NumEl(params, props, args...) = props[:cf_1_num, :cum] / props.num_particles
NumIon(params, props, args...) = props[:cf_3_num, :cum] / props.num_particles
NumExc(params, props, args...) = props[:cf_4_num, :cum] / props.num_particles
function AvgPsEnergy(params, props, args...)
    sum((props[:cf_2_num, (:eps,:cum)] .* EpsGrid(params))[2:end-1]) / props[:cf_2_num, :cum]
end
function AvgThermalEnergy(params, props, args...)
    sum((props[:cf_5_num, (:eps,:cum)] .* EpsGrid(params))[2:end-1]) / props[:cf_5_num, :cum]
end
TimeToThermal(params, props, args...) = props[:cf_5_time, :cum]
TimeToPs(params, props, args...) = props[:cf_2_time, :cum]

FinalNumLeft(params, props, args...) = props[:duration, (:cum,:t)][end] / DT(params)[end]

filename_elements = ["T", "Q", "ion_A", "ion_lam", "exc_A", "exc_thresh", "exc_lam", "init", "spread", "recoil"]


data = cd("runs_PsFormation") do
    GetAll(PsFrac, NumEl, NumIon, NumExc, AvgPsEnergy, AvgThermalEnergy, TimeToThermal, TimeToPs, :numinit, FinalNumLeft, NameElementForGetAll.(filename_elements)..., include_prefix=false, quants_kwds=[:include_errors=>false])
end

labels = [:PsFrac ; :NumEl ; :NumIon ; :NumExc ; :AvgPsEnergy ; :AvgThermalEnergy ; :TimeToThermal ; :TimeToPs ; :numinit ; :FinalNumLeft ; Symbol.(filename_elements)]
#df = DataFrame(permutedims(data), labels)

using JLD2,FileIO
@save "analyse_PsFormation.jld2" data labels


