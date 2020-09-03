import Hardsphere

@testset "Temperature $tmtr" for tmtr in [nothing,293*u"K"]
    params = Hardsphere.SetupParams(tmtr)
    params = Finalise(params)
    props = BunchedPropagate(params,10)

    quants = Quants(params, props, include_errors=false)

    if tmtr == 293u"K"
        @test 0.72eV < quants[:energy] < 0.77eV
        @test 5000u"m/s" < quants[:bulk_W] < 6000u"m/s"
    elseif tmtr == nothing
        @test 0.70eV < quants[:energy] < 0.75eV
        @test 5000u"m/s" < quants[:bulk_W] < 6000u"m/s"
    end
end
