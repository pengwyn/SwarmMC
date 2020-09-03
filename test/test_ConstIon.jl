import ConstIon

@testset "GNS $gns, $ion_style, $ionweight" for gns in [GNS_REGEN_ALL(), GNS_DOUBLE(), GNS_UPDATE_LOG2FAC()],
    ion_style in [ISS_AFE(), ISS_FRACTION()],
    ionweight in [true,false]

    params = Finalise(ConstIon.SetupParams(;gns, ion_style, ionweight))
    props = BunchedPropagate(params,10)
end