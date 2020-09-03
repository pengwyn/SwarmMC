module SwarmMC

##############

using Requires

# Reexporting for convenience
using Reexport

@reexport using AndExport
using AutoParameters
using MsgWrap

include("extmodules/DanUtils.jl")
@reexport using .DanUtils

include("extmodules/Constants.jl")
@reexport using .Constants

# # This is a little odd - I don't know how to do this properly?
Unitful.register(Constants)

using Statistics
using LinearAlgebra
using EllipsisNotation
using DelimitedFiles
using StatsBase

using StaticArrays
using ArgCheck
using UnPack

using Dierckx: Spline1D
include("extmodules/UnitfulDierckx.jl")
using .UnitfulDierckx
using QuadGK

using JLD2, FileIO

# FIXME: Need to replace all of these with searchsortedlast
include("extmodules/BisectInterp.jl")
using .BisectInterp

# Need these for ease of JLD2 loading (hopefully not anymore with ObjectSaving.jl)
# export StaticArrays, BisectInterp, Dierckx

###############

load_dir = "SwarmMC_include/"

include(load_dir * "SwarmMCCommon.jl")
include(load_dir * "SwarmMCErrors.jl")
include(load_dir * "SwarmMCProps.jl")
include(load_dir * "SwarmMCCollFreq.jl")
include(load_dir * "SwarmMCCollision.jl")
include(load_dir * "SwarmMCEvolution.jl")
include(load_dir * "SwarmMCMeasurements.jl")
include(load_dir * "SwarmMCInitStyles.jl")
include(load_dir * "SwarmMCFinalisation.jl")
include(load_dir * "SwarmMCBunchedPropagate.jl")
include(load_dir * "SwarmMCSave.jl")
include(load_dir * "SwarmMCAnalyse.jl")
include(load_dir * "SwarmMCLoaders.jl")

###############

function __init__()
    Unitful.register(Constants)
    merge!(Unitful.basefactors, localunits)

    @require NonNegLeastSquares="b7351bd1-99d9-5c5d-8786-f205a815c4d7" include("SwarmMCLoaderAngleDist.jl")
end


end
