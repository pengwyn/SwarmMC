

using Generic

using CurrentTrace

#cell_width = 1e-4
cell_width = 2e-3

# @CheckTurns for as_loss in [false, true],
# #@CheckTurns for as_loss in [true],
#                 delay_type in [:fixed_dt, :exp, :fractional],
#                 #delay_type in [:fractional],
#                 delay_median in [1e6, 1e7, 1e5],
#                 trap_mag in [20., 1., 0.01, 0.],
#                 return_eps in [0., :thermal, 1.],
# 				Eorig in [92.3, 52.9]

# @CheckTurns for as_loss in [false],
#                 delay_type in [:fractional],
#                 delay_median in [1e6],
#                 trap_mag in [1.],
#                 return_eps in [0., :thermal, 1.],
#                 Eorig in [92.3, 52.9]

# @CheckTurns for as_loss in [true],
#                 delay_type in [:fixed_dt],
#                 delay_median in [1e5],
#                 trap_mag in [0.1, 1., 10.],
#                 return_eps in [0.],
# 				Eorig in [92.3, 52.9, 9.23, 923.]

# @CheckTurns for as_loss in [true],
#                 delay_type in [:fixed_dt],
#                 delay_median in [1e5],
#                 trap_mag in [0.1],
#                 return_eps in [0.],
# 				Eorig in [9.23, 92.3, 923.]

@CheckTurns for as_loss in [false],
                delay_type in [:fractional],
                delay_median in [1e6,1e7],
                trap_mag in [1.0],
                return_eps in [1.],
				Eorig in [92.3]
    
    if as_loss
        if !(return_eps == 0.) || !(delay_type == :fixed_dt) || !(delay_median == 1e5)
            continue
        end
    end

    if !as_loss && trap_mag == 0.
        continue
    end
    
    p = SetupParams(solvation_as_loss=as_loss,
                    delay_type=delay_type,
                    delay_median=delay_median,
                    trap_mag=trap_mag,
                    cell_width_SI=cell_width,
                    return_eps=return_eps,
                    Eorig=Eorig)

    props = LoopMaxTime(p, 10)

    Save(p, props)
end
