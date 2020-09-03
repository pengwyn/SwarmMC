
using DanUtils
using NessRobsonLoss

@CheckTurns for (lossp,alist) = [(0.5, [0, 1e-5, 1e-1, 10]),
                                 (-0.5, [1e-5, 0.1, 1]),
                                 (-1, [1e-6, 5e-3])]
    for lossa = alist,
        Temp = [293.]   
        #for ans = [ANS_SUBSTEPS(), ANS_FAKE_NONCONS(), ANS_NOTHING()],
        for ans = [ANS_SUBSTEPS(), ANS_NOTHING()],
            #gns = [GNS_DOUBLE(), GNS_REGEN_ALL(), GNS_NOTHING(), GNS_UPDATE_LOG2FAC()],
            gns = [GNS_DOUBLE(), GNS_REGEN_ALL(), GNS_NOTHING()],
            numpart = [1, 10, 1000],
            doreduction = [false, true],
            split_fake = (ans == ANS_FAKE_NONCONS() ? [false, true] : [false])

            if gns == GNS_DOUBLE()
                if doreduction
                    continue
                end
            elseif gns == GNS_REGEN_ALL()
            else # This is UPDATE_LOG2FAC or NOTHING
                if ans == ANS_FAKE_NONCONS() && doreduction
                    continue
                elseif !doreduction # This is SUBSTEPS or NOTHING
                    continue
                end
            end


            if ans == ANS_FAKE_NONCONS() && numpart < 10
                continue
            end

            if doreduction
                weight_reduction = 0.5
            else
                weight_reduction = 1.0
            end

            local p
            try
                println("$lossp, $lossa, $ans, $gns, $numpart, $doreduction, $split_fake")

                p = SetupParams(lossa, lossp, ans, gns, weight_reduction, split_fake, Temp)
                p.save_name = p.save_name * ":N=$(numpart)"
                p = SwarmMC.Finalise(p)
                props = LoopMaxTime(p, numpart)

                Save(p, props)
            catch exc
                print_with_color(:red, "Got an exception ($(string(exc))) with $(p.save_name)", bold=true)
                if !isa(exc, RemoteException) || !isa(exc.captured.ex, ErrorException)
                    rethrow(exc)
                end
            end
        end
    end
end
