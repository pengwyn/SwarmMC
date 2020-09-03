import Maxwell
import MagneticField
using SwarmMC

function FuncWrapper()
    N = 10
    meas_type_combs = [(MEAS_TYPE(:t),),
                       #(MEAS_TYPE(:cum),),
                       #(MEAS_TYPE(:ss,:cum),),
                       (MEAS_TYPE(:ss,:cum,:eps),),
                       (MEAS_TYPE(:t,:cum), MEAS_TYPE(:eps,)),
                       (MEAS_TYPE(:eps,:t,:z),),
                       #(MEAS_TYPE(:cum,:eps,:t),)
                       ]
                       

    programs = [("Maxwell", () -> Maxwell.SetupParams(293.)),
                ("MagneticField", () -> MagneticField.SetupParams(1000.)),
                ]

    for meas_types in meas_type_combs,
        (name,program) in programs

        #params = Maxwell.SetupParams(293.)
        params = program()
        params.meas_types = meas_types
        params = FinaliseParams(params)

        BunchedPropagate(params, 2)
        
        timing = @timed BunchedPropagate(params, N)

        time = timing[2]
        bytes = timing[3]

        open("speed_results.txt", "a") do file
            write(file, "$(now()), $name(N=$N), $(meas_types), $(round(time,signif=3)), $(signif(bytes / 1024^3, 3)) GB\n")
        end
    end
end

FuncWrapper()
