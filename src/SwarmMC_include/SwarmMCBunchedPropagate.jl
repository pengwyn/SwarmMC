# This is the wrapper to finalise the params
@xport BunchedPropagate(params, N, verbose=TypeTrue() ; kwds...) = BunchedPropagateInternal(Finalise(params), N, verbose; kwds...)
function BunchedPropagateInternal(params, N, verbose=TypeTrue() ; return_part_list=TypeFalse())
    props_out = PROPS_OUT(params, N*uN)
    walltime_start = time()
    
    int = SetupIntegrator(params)

    part_list = [GenerateInitParticle(params.init_style, params) for i in 1:N]
    if IsTrue(return_part_list)
        init_part_list = deepcopy(part_list)
    end

    # For debugging/early failure
    num_redos_fakenoncons = 0

    for step_ind in 0:length(params.t_grid)-1
        if isempty(part_list)
            @warn "Exiting because of empty part_list - this shouldn't happen anymore"
            break
        end

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

        if IsTrue(verbose) && @printagain(10)
            timepercentage = binstart / params.t_grid[end] * 100
            print("Timestep $step_ind/$(length(params.t_grid)-1) ($(round(timepercentage,sigdigits=3))%) ")
        end

        setfield!.(part_list, :t_bin, step_ind)

        abort = false
        recommendation,new_list = BunchedOneStep(int, params, N, part_list, binstart, timestep, props_out, verbose=verbose)

        part_list = new_list

        if recommendation == :double_substeps
            error("Not allowing this any more")
        elseif recommendation == :finished
            abort = true
        elseif recommendation == :success
        else
            error("Unknown recommendation $recommendation.")
        end

        if abort
            break
        end

    end

    props_out.walltime = (time() - walltime_start) * u"s"

    if IsTrue(return_part_list)
        return props_out, part_list, init_part_list
    else
        return props_out
    end
end

function BunchedOneStep(int, params, target_N, part_list, start_time, timestep, props_out ; verbose=TypeTrue()::TypeBool)
    cur_time = start_time
    end_time = start_time + timestep
    num_substeps = max(1, ceil(Int, timestep / params.max_substepsize))


    for substep_ind in 1:num_substeps
        if isempty(part_list)
            break
        end

        if IsTrue(verbose) && @printagain(10)
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
            println("Substep: ($substep_ind/$num_substeps). Mean eps: $(round(avgeps,sigdigits=3)). Mean z: $(round(avgz,sigdigits=3)). Tot weight: $totweight")
        end

        substep_end_time = min(cur_time + params.max_substepsize, end_time)
        todo_list = copy.(part_list)
        result,new_list = BunchedSubStep(int, params, target_N, todo_list, substep_end_time, props_out, verbose=verbose)

        before_cull = new_list
        new_list = filter(x->x.weight > 0uN, new_list)

        ratio = length(new_list) / length(part_list)

        if result == :got_too_many
            if params.adapt_noncons_style == ANS_NOTHING()
                error("No way to adapt for too many particles with $(params.adapt_noncons_style).")
            elseif params.adapt_noncons_style == ANS_FAKE_NONCONS()
                error("Not implemented")
            else
                error("Unknown ANS")
            end
        elseif result != :success
            error("Unknown result $result")
        end

        if ratio < 0.1 || length(new_list) == 0
            if params.adapt_noncons_style == ANS_NOTHING()
                if length(new_list) == 0
                    @debug "Out of particles with no adapt style."
                    return :finished,new_list
                end
                # Otherwise just keep going
            elseif params.adapt_noncons_style == ANS_FAKE_NONCONS()
                error("Not implemented")
            else
                error("Unknown ANS")
            end
        end

        cur_time = substep_end_time

        if length(new_list) == 0
            @show length(part_list)
        end

        before_N = length(new_list)
        part_list = GenerateNextStep(params.generate_next_step_style, new_list, target_N)
        after_N = length(part_list)
    end

    for part in part_list
        # bins = (part.t_bin, part.r_bin, part.z_bin, part.eps_bin, 0)
        steady_state = part.time > params.t_grid[end]*params.steady_state_timefrac

        for meas in props_out.meas
            UpdateMeasInst!(meas, part, steady_state)
        end
    end

    return :success,part_list
end

function BunchedSubStep(int, params, target_N, todo_list, target_time, props_out ; verbose=TypeTrue()::TypeBool)
    new_list = PARTICLE[]
    init_todo_N = length(todo_list)

    count_alive = 0
    while !isempty(todo_list)
        if IsTrue(verbose) && @printagain(10)
            println("Particles to go $(length(todo_list))/$(init_todo_N).")
        end

        part = pop!(todo_list)
        extra_list = GoToTime(int, params, part, props_out, target_time)
            
        if part.weight > 0uN
            push!(new_list, part)
            count_alive += 1
        end

        append!(todo_list, extra_list)

        if count_alive > params.runaway_threshold*target_N
            @info "Runaway number of particles!" count_alive target_N init_todo_N
            # Return new_list anyway just in case it helps
            return :got_too_many,new_list
        end
    end
    return :success,new_list
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


