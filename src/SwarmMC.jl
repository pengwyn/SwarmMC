module SwarmMC

##############

using Requires

# Reexporting for convenience
using Reexport

@reexport using AndExport
using AutoParameters
using MsgWrap

include("extmodules/DanUtilsInternal.jl")
@reexport using .DanUtilsInternal

include("extmodules/ConstantsInternal.jl")
@reexport using .ConstantsInternal

# # This is a little odd - I don't know how to do this properly?
Unitful.register(ConstantsInternal)

using Statistics
using LinearAlgebra
using EllipsisNotation
using DelimitedFiles
using StatsBase

using StaticArrays
using ArgCheck
using UnPack

using Dierckx: Spline1D
include("extmodules/UnitfulDierckxInternal.jl")
using .UnitfulDierckxInternal
using QuadGK

using JLD2, FileIO

# FIXME: Need to replace all of these with searchsortedlast
include("extmodules/BisectInterpInternal.jl")
using .BisectInterpInternal

# Need these for ease of JLD2 loading (hopefully not anymore with ObjectSaving.jl)
# export StaticArrays, BisectInterp, Dierckx

###############

load_dir = "SwarmMC_include/"

include(joinpath(load_dir, "SwarmMCCommon.jl"))
include(joinpath(load_dir, "SwarmMCErrors.jl"))
include(joinpath(load_dir, "SwarmMCProps.jl"))
include(joinpath(load_dir, "SwarmMCCollFreq.jl"))
include(joinpath(load_dir, "SwarmMCCollision.jl"))
include(joinpath(load_dir, "SwarmMCEvolution.jl"))
include(joinpath(load_dir, "SwarmMCMeasurements.jl"))
include(joinpath(load_dir, "SwarmMCInitStyles.jl"))
include(joinpath(load_dir, "SwarmMCFinalisation.jl"))
include(joinpath(load_dir, "SwarmMCBunchedPropagate.jl"))
include(joinpath(load_dir, "SwarmMCSave.jl"))
include(joinpath(load_dir, "SwarmMCAnalyse.jl"))
include(joinpath(load_dir, "SwarmMCLoaders.jl"))

###############

function __init__()
    Unitful.register(ConstantsInternal)
    merge!(Unitful.basefactors, localunits)

    @require NonNegLeastSquares="b7351bd1-99d9-5c5d-8786-f205a815c4d7" include(joinpath(load_dir, "SwarmMCLoaderAngleDist.jl"))
    @require MAT="23992714-dd62-5051-b70f-ba57cb901cac" include(joinpath(load_dir,"SwarmMCLoaderMAT.jl"))
end


end
