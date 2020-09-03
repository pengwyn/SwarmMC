
using SwarmMC
using Generic

using PsFormation

if "Set" ∈ keys(ENV)
    @assert ENV["Set"] |> lowercase ∈ ["a","b","c","dist","benchmark"]
    @info "Using only set " * ENV["Set"]
end

#coarse_temperature = [0., 100., 293., 1e3, 1e4, 5e4, 1e5, 5e5, 1e6]
#coarse_temperature = [0., 293., 8000.]
coarse_temperature = [0.]
# coarse_init_eps = [1e3, 1e4, 1e5]
coarse_init_eps = [1e4]
# coarse_eps_spread = [:gauss_twiceth, :gauss_tenth]
# coarse_eps_spread = [:gauss_twiceth]
coarse_eps_spread = [:gauss_tenth]
# coarse_Q = [0.5, 0.8, 0.9, 0.99, 1.0, :hydrogen_dist, :hydrogen_mean]
# coarse_Q = [0.5, 0.8, 0.9, 0.99, 1.0, :hydrogen_dist]
# coarse_Q = [0.5, 1.0, :hydrogen_dist]
coarse_Q = [1.0, :hydrogen_dist]
# coarse_Q = [:hydrogen_dist, :hydrogen_mean]
# coarse_recoil = [true,false]
coarse_recoil = [true]



# Benchmark numbers
if get(ENV, "Set", "benchmark") |> lowercase == "benchmark"
    recoil = true
    init_eps = 1e4
    eps_spread = :gauss_tenth
    temperature = 0.

    # Set A
    ion_λ = 30.
    ion_A = 0.3
    exc_A = 0.
    exc_threshold = 0.
    exc_λ = 1.
    Q = 1.0

    @CheckTurns for Q in [0.5, 0.8, 0.9, 0.95, 0.99, 1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for ion_A in [0.1, 0.5, 1.0]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for ion_λ in [10., 20., 40.]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end


    # Set B

    ion_λ = 30.
    ion_A = 0.3
    exc_A = 0.5
    exc_threshold = 10.
    exc_λ = 12.
    Q = :force_check

    @CheckTurns for exc_A in [0.1,0.5,1.0,1.5],
                    Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_threshold in [20.,30.],
                    Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_λ in [10.,20.],
                Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    # Set C

    ion_λ = 30.
    ion_A = 0.3
    exc_A = 1.0
    exc_threshold = 2.
    exc_λ = 10.
    Q = :force_check

    @CheckTurns for exc_A in [0.2,0.6,1.0,1.4,2.0],
                    Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_threshold in [1.,3.,6.],
                    Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_λ in [6.,20.],
                    Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end
end




# Set A distributions
if get(ENV, "Set", "dist") |> lowercase == "dist"
    ion_λ = 30.
    ion_A = 0.3
    exc_A = 0.
    exc_threshold = 0.
    exc_λ = 1.
    temperature = 0.
    init_eps = 1e4
    eps_spread = :gauss_tenth
    recoil = false

    @CheckTurns for Q in [1.0, :hydrogen_dist]
        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end
end

# Set A
if get(ENV, "Set", "a") |> lowercase == "a"
    ion_λ = 30.
    ion_A = 0.3
    exc_A = 0.
    exc_threshold = 0.
    exc_λ = 1.

    @CheckTurns for Q in unique([0.1:0.05:0.5 ; 0.5:0.005:0.8 ; 0.8:0.01:1.0]),
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for ion_A in 0.1:0.05:1.5,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for ion_λ in 10:2:50,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end
end


# Set B
if get(ENV, "Set", "b") |> lowercase == "b"

    ion_λ = 30.
    ion_A = 0.3
    exc_A = 0.5
    exc_threshold = 10.
    exc_λ = 12.

    @CheckTurns for exc_A in 0.1:0.1:1.5,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_threshold in 6:0.5:30,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_λ in 5:1:25,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end
end

# Set C
if get(ENV, "Set", "c") |> lowercase == "c"

    ion_λ = 30.
    ion_A = 0.3
    exc_A = 1.0
    exc_threshold = 2.
    exc_λ = 10.

    @CheckTurns for exc_A in 0.2:0.2:2.0,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_threshold in 0.5:0.5:6.0,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for exc_λ in 4:2:20,
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        temperature in coarse_temperature,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

    @CheckTurns for temperature in [0.0:50.0:500.0 ; 1000.0:1000.0:10000.0 ; 2e4:1e4:1e5],
        Q in coarse_Q,
        init_eps in coarse_init_eps,
        eps_spread in coarse_eps_spread,
        recoil in coarse_recoil

        p = SetupParams(ion_A, ion_λ, exc_A, exc_threshold, exc_λ, Q, temperature=temperature, init_eps=init_eps, eps_spread=eps_spread, recoil=recoil)
        props = LoopMaxTime(p, 1000)

        Save(p, props)
    end

end
