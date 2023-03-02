module RomeoDFT

using Reexport
@reexport using DFControl
@reexport using RemoteREPL
using RemoteHPC
using RemoteHPC: @timeout, suppress
using DFWannier
using JLD2: JLD2
using Distances
using Dates
using LinearAlgebra
using EulerAngles
using ThreadPools
using LoggingExtras
using Optim
using Sockets
using SnoopPrecompile
using ProgressMeter
using PrettyTables

@reexport using Overseer


const AnglesType = Angles{Float64,2,DFWannier.MagneticVector{Float64, Vector{Float64}}}
const ColinMatrixType = DFW.ColinMatrix{Float64, Matrix{Float64}}
const MagneticVectorType = DFW.MagneticVector{Float64, Vector{Float64}}

const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "RomeoDFT") :
                   abspath(Base.DEPOT_PATH[1], "config", "RomeoDFT")

const SEARCHERS_DIR = Ref("")

config_path(p...) = joinpath(CONFIG_DIR, p...)
searchersdir(p...) = joinpath(SEARCHERS_DIR[], p...)

include("utils.jl")
include("states.jl")

include("jobs.jl")

@enum TrialOrigin RandomMixed EulerAngleMixed LinearMixed PostProcess IntersectionMixed Unknown 
@enum MixingMode RandomMixing EulerAngleMixing LinearMixing UnknownMixing

include("database.jl")

include("Systems/core.jl")
include("Systems/postprocessing.jl")
include("Systems/firefly.jl")
include("Systems/intersection.jl")
include("Systems/structural.jl")
include("search.jl")
include("orchestrator.jl")
include("client.jl")
include("analysis.jl")
include("CLI/cli.jl")

export Searcher, connect_orchestrator
export ground_state, unique_states

# @precompile_all_calls begin
#     tn = tempname()
#     l = Searcher(tn)
#     Entity(l, ServerInfo("", "", "", 5))
#     save(l)
#     load(l)
#     rm(tn, recursive=true)

#     isalive(local_server())
# end

using Requires
using TOML

set_searchers_dir() = 
    SEARCHERS_DIR[] = TOML.parsefile(config_path("config.toml"))["searchers_directory"]

function __init__()
    set_searchers_dir()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        @eval include("Systems/bandsplotter.jl")
        @require LaTeXStrings="b964fa9f-0449-5b57-a5c2-d3ea65f4040f" @eval include("plotting.jl")
    end
end

end # RomeoDFT
