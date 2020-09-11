
"""
BunchedPropagate(params, N ; return_part_list=false)

Run a simulation of `N` particles using `params`. Returns a `PROPS_OUT` object.

If you pass true to `return_part_list` then the initial and final particle lists
will also be returned.
"""
@xport function BunchedPropagate(params, N, verbose=TypeTrue() ; kwds...)
    params = Finalise(params)
    int = SetupIntegrator(params)
    BunchedPropagateInternal(params, N, int, verbose; kwds...)
end
    
function BunchedPropagateInternal(params, N, int, verbose=TypeTrue() ; return_part_list=TypeFalse())
    props_out = PROPS_OUT(params, N*uN)
    walltime_start = time()

    part_list = [GenerateInitParticle(params.init_style, params) for i in 1:N]
    todo_list = eltype(part_list)[]

    if IsTrue(return_part_list)
        init_part_list = deepcopy(part_list)
    end

    for step_ind in 0:length(params.t_grid)-1
        if step_ind == 0
            binstart = minimum(x -> x.time, part_list)
            binend = params.t_grid[1]
            if binend <= binstart
                continue
            end
        else
            binstart,binend = params.t_grid[step_ind:step_ind+1]
        end

        timestep = binend - binstart

        setfield!.(part_list, :t_bin, step_ind)

        cur_time = binstart
        # end_time = start_time + timestep
        num_substeps = max(1, ceil(Int, timestep / params.max_substepsize))

        for substep_ind in 1:num_substeps
            isempty(part_list) && @goto endofloops

            if IsTrue(verbose) && @printagain(10)
                timepercentage = cur_time / params.t_grid[end] * 100
                if any(isfinite.([EpsFromVel(x) for x in part_list]))
                    avgeps = mean(x for x in EpsFromVel.(part_list) if isfinite(x))
                else
                    avgeps = -1.
                end

                avgvz = mean(x.vel[3] for x in part_list)
                avgz = mean(x.pos[3] for x in part_list)
                totweight = sum(x.weight*exp2(x.log2fac) for x in part_list)
                if totweight == 0.0uN
                    totweight = "≈2^" * string(sum(x.log2fac for x in part_list))
                else
                    totweight = round(totweight,sigdigits=3)
                end
                println("Timestep $step_ind/$(length(params.t_grid)-1) ($(round(timepercentage,sigdigits=3))%) Substep: ($substep_ind/$num_substeps). Mean eps: $(round(avgeps,sigdigits=3)). Mean z: $(round(avgz,sigdigits=3)). Tot weight: $totweight")
            end

            substep_end_time = min(cur_time + params.max_substepsize, binend)
            orig_N = length(part_list)
            todo_list,part_list = part_list,todo_list

            @assert isempty(part_list)

            count_seen = 0
            while !isempty(todo_list)
                if IsTrue(verbose) && @printagain(10)
                    println("Particles to go $(length(todo_list))/$(orig_N).")
                end

                part = pop!(todo_list)
                extra_list = GoToTime(int, params, part, props_out, substep_end_time)
                
                if part.weight > 0uN
                    push!(part_list, part)
                    count_seen += 1
                end

                append!(todo_list, extra_list)

                if count_seen > params.runaway_threshold*N
                    @error "Runaway number of particles!" count_seen N orig_N
                    error("Runaway number of particles")
                end
            end

            @assert all(x->x.weight > 0uN, part_list)

            ratio = length(part_list) / orig_N

            if ratio < 0.1 || length(part_list) == 0
                if length(new_list) == 0
                    @info "Out of particles with no adapt style."
                    @goto endofloops
                end
            end

            cur_time = substep_end_time

            part_list = GenerateNextStep(params.generate_next_step_style, part_list, N)
        end

        # Instantaneous measurements
        for part in part_list
            # bins = (part.t_bin, part.r_bin, part.z_bin, part.eps_bin, 0)
            steady_state = part.time > params.t_grid[end]*params.steady_state_timefrac

            for meas in props_out.meas
                UpdateMeasInst!(meas, part, steady_state)
            end
        end
    end

    @label endofloops

    props_out.walltime = (time() - walltime_start) * u"s"

    if IsTrue(return_part_list)
        return props_out, part_list, init_part_list
    else
        return props_out
    end
end

######################################
# * GenerateNextStep
#------------------------------------

GenerateNextStep(style::GNS_NOTHING, part_list, target_N) = part_list

function GenerateNextStep(style::GNS_UPDATE_LOG2FAC, part_list, target_N)
    for x in part_list
        fac = floor(Int, log2(ustrip(uN,x.weight)))
        x.log2fac += fac
        x.weight /= exp2(fac)
    end
    part_list
end

using Random
function GenerateNextStep(style::GNS_DOUBLE, part_list::Vector{PARTICLE}, target_N::Int)
    if isempty(part_list)
        @debug "Tried to double an empty list"
        return part_list
    end

    while length(part_list) < target_N÷2
        foreach(x -> x.log2fac -= 1, part_list)
        append!(part_list, deepcopy(part_list))

        @debug "Doubled list" length(part_list) target_N sum(x -> x.weight, part_list)
    end

    # The opposite is also possible, but I need to randomly choose particles to throw away.
    # This could also be done to bring the number back to the target, but I will keep with doubling for the moment.
    while length(part_list) > target_N*2
        AllSame(v) = all(==(first(v)), v)
        @assert AllSame(getproperty.(part_list, :weight)) && AllSame(getproperty.(part_list, :log2fac)) "GNS_DOUBLE() with ionisation only works if the weights are the same. Otherwise, should be using GNS_REGEN_ALL() instead"

        # Randomly choose half of them to remove.
        init_len = length(part_list)

        shuffle!(part_list)

        final_len = init_len÷2
        part_list = part_list[1:final_len]

        fac = init_len / final_len
        log2fac = floor(Int, log2(fac))
        fac /= exp2(log2fac)
        for x in part_list
            x.log2fac += log2fac
            x.weight *= fac
        end
    end

    part_list
end

function GenerateNextStep(style::GNS_REGEN_ALL, part_list::Vector{PARTICLE}, target_N::Int)
    # PARTICLE FILTEREING???
    if isempty(part_list)
        return part_list
    end
    if length(part_list) == 1
        return part_list
    end

    for x in part_list
        fac = floor(Int, log2(ustrip(uN,x.weight)))
        x.log2fac += fac
        x.weight /= exp2(fac)
    end

    # Only regen if we have a skewed distribution of particles so that the
    # effective number is much less/greater than target_N.

    preweight = ustrip.(uN, getproperty.(part_list, :weight))

    log2w_list = getproperty.(part_list, :log2fac) + log2.(preweight)
    common_fac = maximum(log2w_list)
    log2w_list .-= common_fac

    # This is <w>^2 / <w>^2
    eff_num = sum(exp2.(log2w_list))^2 / sum(exp2.(log2w_list .+ 1))

    if eff_num > 1 && (target_N÷10 < eff_num < target_N*10)
        return part_list
    end

    # @info "Going to regen particle list" eff_num target_N
    
    log2_tot = log2(sum(exp2.(log2w_list)))
    orig_weight_tot = log2_tot + common_fac

    ratios = exp2.(log2w_list .- log2_tot)

    target_log2fac = log2_tot + common_fac - log2(target_N)
    to_set_log2fac = floor(Int, target_log2fac)
    to_set_weight = exp2(target_log2fac - to_set_log2fac) * uN

    new_list = PARTICLE[]
    w = StatsBase.weights(ratios)
    for i = 1:target_N
        old_ind = sample(eachindex(w), w)
        part = copy(part_list[old_ind])

        part.log2fac = to_set_log2fac
        part.weight = to_set_weight

        push!(new_list, part)
    end

    preweight = ustrip.(uN, getproperty.(new_list, :weight))
    log2w_list = getproperty.(new_list, :log2fac)
    common_fac = maximum(log2w_list)
    log2w_list .-= common_fac
    new_weight_tot = log2(sum(preweight .* exp2.(log2w_list))) + common_fac

    new_list
end


####################################
# * Parallel things
#----------------------------------

using Distributed


function ParallelBunchedPropagate(args...)
    global globparams
    Base.invokelatest(BunchedPropagate, globparams, args...)
end

function SetupParallelBunchedPropagate(p)
    global globparams
    globparams = Base.invokelatest(Finalise, p)
    nothing
end

"""
LoopMaxTime(params, N ; walltime=:auto)

Runs simulations repeatedly until `walltime` (a `Unitful` quantity) is reached.

If `walltime` is `:auto` then try to lookup the environment variable
`SWARMMC_WALLTIME` and interpret it as a time in minutes.

`N` and `params` are passed through to `BunchedPropagate`.
"""
@xport function LoopMaxTime(p, args... ; walltime = :auto, multiprocessing=true)
    if walltime == :auto
        if "SWARMMC_WALLTIME" in keys(ENV)
            walltime = parse(Float64, ENV["SWARMMC_WALLTIME"]) * u"minute"
        else
            walltime = 1.0 * u"minute"
        end
    end
    
    unittime() = time()*u"s"
    
    starttime = unittime()

    npercycle = 1
    numdone = 0

    props = nothing
    ss = nothing

    if multiprocessing
        @everywhere $SetupParallelBunchedPropagate($(deepcopy(p)))
    else
        p = Finalise(p)
    end

    while numdone == 0 || unittime() - starttime < walltime

        if multiprocessing
            thisprops = @distributed (+) for i in 1:npercycle*nworkers()
                ParallelBunchedPropagate(args...)
            end
        else
            thisprops = BunchedPropagate(p, args...)
        end

        if props == nothing
            props = thisprops
        else
            props += thisprops
        end

        numdone += npercycle
        timetaken = unittime() - starttime

        timepern = timetaken / numdone

        timeremaining = walltime - timetaken

        if @printagain()
            println("Have done $numdone cycles. Time to go: $(PrettyTime(ustrip(u"s",timeremaining))).")
        end
        
    end

    props
end


