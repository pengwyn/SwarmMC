
include("NessRobsonLoss.jl")
using .NessRobsonLoss

@CheckTurns for (lossp,alist) = [(0.5, [0,1e-5,1e-4,1e-3,1e-2,1e-1,1,10]*u"Å^2/sqrt(eV)"),
                      (-0.5, [0,1e-5,1e-2,1]*u"Å^2*sqrt(eV)"),
                      (-1, [0,1e-6,1e-5,1e-4,3e-4,5e-4,7e-4,9e-4,1e-3,2e-3,5e-3]*u"Å^2*eV")]
    for lossa = alist
        # tmtr = 0.
        tmtr = 293.0u"K"

		gns_style = GNS_REGEN_ALL()
		weight_reduction = 0.5

		# gns_style = GNS_DOUBLE()
		# weight_reduction = 1.0

		num_part = 100
		split_fake = false

        @show lossp, lossa

        p = SetupParams(; lossa, lossp, gns_style, weight_reduction,
                        # split_fake,
                        tmtr)
        p.save_name = p.save_name * ":N=$(num_part)"
        p = SwarmMC.Finalise(p)
        props = LoopMaxTime(p, num_part)

        Save(p, props)
    end
end
