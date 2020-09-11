import Hardsphere

@testset "Temperature $tmtr" for tmtr in [nothing,293*u"K"]
    params = Hardsphere.SetupParams(tmtr)
    params = Finalise(params)
    props = BunchedPropagate(params,10)

    quants = Quants(params, props, include_errors=false)

    if tmtr == 293u"K"
        @test 0.65eV < quants[:energy] < 0.85eV
        @test 4000u"m/s" < quants[:bulk_W] < 7000u"m/s"
    elseif tmtr == nothing
        @test 0.65eV < quants[:energy] < 0.85eV
        @test 4000u"m/s" < quants[:bulk_W] < 7000u"m/s"
    end
end
