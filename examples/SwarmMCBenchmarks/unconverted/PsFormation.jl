
module PsFormation

using Reexport ; @reexport using SwarmMC
using DanUtils, Constants
using Kwonly

const ion_threshold = 13.6
const ps_threshold = ion_threshold - 6.8
#const exc_threshold = 8.

function SurgeExp(eps, threshold, p, λ)
    x = eps - threshold
    A = (ℯ/λ)^p
    return A * x^p * exp(-p*x / λ)
end
function SurgePower(eps, threshold, p, λ)
    x = eps - threshold
    A = λ^p * (p+1)^(p+1)
    return A * x / (x + p*λ)^(p+1)
end


Elastic(eps) = 0.4 * SurgePower(eps, -20, 2, 20)
PsForm(eps) = SurgeExp(eps, 6.8, 0.5, 7)
Ionisation(eps, A, λ) = A * SurgePower(eps, 13.6, 1, λ)
Excitation(eps, A, threshold, λ) = A * SurgePower(eps, threshold, 1.5, λ)

FakeAbsorb(eps) = 1e10 * (eps < 6.8)


# Note that the formula is calculated from incident energy, but the function is
# called with the energy after subtracting the ionisation threshold
function HydrogenSharingFuncDist(eps,ion=13.6,R=rand())
    m1 = -2.1 + 1

    # Work in ionisation energies
    eps = eps/ion + 1
    eps2 = (0.5^m1 * (1-R) + (eps - 0.5)^m1*R)^(1/m1) - 0.5

    return eps2*ion
end

using QuadGK
using Dierckx
HydrogenSharingMean(eps,ion=13.6) = quadgk(x -> HydrogenSharingFuncDist(eps,ion,x), 0, 1)[1]

ForwardAngle(en, R) = 1.

function CollFreqs(params, gas, ptype, Q, cs_params)
    temperature = params.temperature[1,1]

    IonFunc = let A=cs_params[:ion_A], λ=cs_params[:ion_λ]
        x->Ionisation(x,A,λ)
    end
    ExcFunc = let A=cs_params[:exc_A], threshold=cs_params[:exc_threshold], λ=cs_params[:exc_λ]
        x->Excitation(x,A,threshold,λ)
    end

    elastic = CreateCollFreq(params, gas, ptype, "elastic", CFS_ELASTIC(), GENCF_MANUAL(Elastic), temperature=temperature)
    psformation = CreateCollFreq(params, gas, ptype, "ps", CFS_LOSS(), GENCF_MANUAL(PsForm), threshold=6.8)
    exc = CreateCollFreq(params, gas, ptype, "excitation", CFS_INELASTIC(), GENCF_MANUAL(ExcFunc), temperature=temperature, threshold=cs_params[:exc_threshold])
    fake_loss = CreateCollFreq(params, gas, ptype, "fake", CFS_LOSS(), GENCF_MANUAL(FakeAbsorb))

    ion_sharing_ratio = 0.
    if Q == :hydrogen_dist
        ion_sharing_style = ISS_FUNC(HydrogenSharingFuncDist)
    elseif Q == :hydrogen_mean
        eps_list = 10 .^ linspace(log10(13.6), 6, 1001)
        mean_list = HydrogenSharingMean.(eps_list)
        spline = Spline1D(eps_list, mean_list, k=1)

        ion_sharing_style=ISS_FUNC(spline)
    else
        ion_sharing_ratio = Q
        ion_sharing_style=ISS_FRACTION()
    end
    ion = CreateCollFreq(params, gas, ptype, "ionisation", CFS_IONISATION(), GENCF_MANUAL(IonFunc), new_ptype=-1, ion_sharing_style=ion_sharing_style, ion_sharing_ratio=ion_sharing_ratio, threshold=13.6, angle_dist_cum=ForwardAngle)

    [elastic, psformation, ion, exc, fake_loss]
end

export SetupParams
@add_kwonly function SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q ; temperature=0., init_eps=1e4, eps_spread=:gauss_twiceth, recoil=true)
    if eps_spread == :gauss_twiceth
        init_deps = 6.8*2
    elseif eps_spread == :gauss_tenth
        init_deps = init_eps / 10.
    else
        error("Unknown spread type")
    end

    # For convenience:
    ion_A = Float64(ion_A)
    ion_λ = Float64(ion_λ)
    exc_A = Float64(exc_A)
    exc_threshold = Float64(exc_threshold)
    exc_λ = Float64(exc_λ)

    p = PARAMS(init_style=SwarmMC.PIS_SEPARABLE(SwarmMC.PIVS_GAUSSIAN_ENERGY(init_eps, init_deps)))
    
    cs_params = Dict(:ion_A => ion_A,
                     :ion_λ => ion_λ,
                     :exc_A => exc_A,
                     :exc_threshold => exc_threshold,
                     :exc_λ => exc_λ)

    # Just overwrite recoil if necessary
    # if Q ∈ [:hydrogen_dist, :hydrogen_mean] && recoil
    #     @warn "Setting recoil to false for hydrogen scattering"
    #     recoil = false
    # end

    if Q isa Real
        Q = Float64(Q)
    end

    if recoil
        m0 = amu / p.mass_unit
    else
        m0 = Inf
        temperature = 0.
    end

    n0 = 1e-6

	gas = GAS{:gas}(m0)
    electron = PARTTYPE{:electron}()
    region = REGION{:region}()

    push!(p.gas_list, gas)
    push!(p.ptype_list, electron)
    push!(p.region_list, region)
    
    p.chars[:dens, :region, :gas] = n0
    p.chars[:temp, :region, :gas] = temperature
    p.chars[:cflist, :gas, :electron] = (args...) -> CollFreqs(args..., Q, cs_params)

    p.save_name = MakeSaveName("PsFormation", T=temperature, recoil=recoil, Q=Q, ion_A=ion_A, ion_lam=ion_λ, exc_A=exc_A, exc_thresh=exc_threshold, exc_lam=exc_λ, init=init_eps, spread=eps_spread)

    p.meas_types = (MEAS_TYPE(:t,:cum), MEAS_TYPE(:cum,:eps), MEAS_TYPE(:cum))

    p.t_grid = GRID(GRID_LOG(), 1e7, 1e10, 1001)
    p.eps_grid = GRID(GRID_LOG(), 1e-5, 1e5, 1001)
    # For fake loss - need to make sure 6.8 is in here.
    p.eps_grid = [p.eps_grid ; 6.8] |> sort

    p.min_log2_weight = -30
    #p.weight_reduction = 0.5

    p.adapt_noncons_style = ANS_NOTHING()
    p.generate_next_step_style = GNS_NOTHING()
    
    p
end

end
