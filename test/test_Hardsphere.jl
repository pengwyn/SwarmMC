import Hardsphere

@testset "Temperature $tmtr" for tmtr in [nothing,293*u"K"]
    params = Hardsphere.SetupParams(tmtr)
    params = Finalise(params)
    props = BunchedPropagate(params,thorough_run ? 50 : 1)

    quants = Quants(params, props, include_errors=false)

    if thorough_run
        if tmtr == 293u"K"
            @test 0.73eV < quants[:energy] < 0.79eV
            @test 5500u"m/s" < quants[:bulk_W] < 6000u"m/s"
        elseif tmtr == nothing
            @test 0.70eV < quants[:energy] < 0.76eV
            @test 5500u"m/s" < quants[:bulk_W] < 6000u"m/s"
        end
    else
        if tmtr == 293u"K"
            @test 0.3eV < quants[:energy] < 1.4eV
            @test 2000u"m/s" < quants[:bulk_W] < 9000u"m/s"
        elseif tmtr == nothing
            @test 0.3eV < quants[:energy] < 1.4eV
            @test 2000u"m/s" < quants[:bulk_W] < 9000u"m/s"
        end
    end
        
end
