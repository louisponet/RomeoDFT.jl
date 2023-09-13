module RomeoDFT

using Reexport
# @recompile_invalidations begin
    @reexport using DFControl
    import DFControl: Calculation
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
    using ProgressMeter
    using PrettyTables
    using TOML
    using TimerOutputs
    using REPL.TerminalMenus

    @reexport using Overseer
    using Overseer: AbstractEntity
# end

const AnglesType         = Angles{Float64,2,DFWannier.MagneticVector{Float64, Vector{Float64}}}
const ColinMatrixType    = DFW.ColinMatrix{Float64, Matrix{Float64}}
const MagneticVectorType = DFW.MagneticVector{Float64, Vector{Float64}}

const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "RomeoDFT") :
                   abspath(Base.DEPOT_PATH[1], "config", "RomeoDFT")

config_path(p...)   = joinpath(CONFIG_DIR, p...)
searchers_dir(p...) = joinpath(TOML.parsefile(config_path("config.toml"))["searchers_directory"], p...)

include("utils.jl")
include("states.jl")
const StateType = State{Float64, MagneticVectorType, ColinMatrixType}

include("jobs.jl")

@enum TrialOrigin RandomMixed EulerAngleMixed LinearMixed PostProcess IntersectionMixed ModelOptimized Unknown 
@enum MixingMode RandomMixing EulerAngleMixing LinearMixing UnknownMixing

include("components.jl")
include("search.jl")
include("Systems/core.jl")
include("Systems/postprocessing.jl")
# include("Systems/firefly.jl")
include("Systems/intersection.jl")
include("Systems/structural.jl")
include("Systems/electrides.jl")
include("Systems/hp.jl")
include("Systems/model.jl")
include("stages.jl")
include("orchestrator.jl")
include("client.jl")
include("analysis.jl")
include("CLI/cli.jl")
include("logging.jl")

export Searcher, connect_orchestrator
export ground_state, unique_states

# @compile_workload begin
#     tn = tempname()
#     l = Searcher(tn)
#     Entity(l, ServerInfo("", "", "", 5))
#     save(l)
#     load(l)
#     rm(tn, recursive=true)

#     isalive(local_server())
# end

using Requires

function __init__()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" begin
        @eval include("Systems/bandsplotter.jl")
        @require LaTeXStrings="b964fa9f-0449-5b57-a5c2-d3ea65f4040f" @eval include("plotting.jl")
    end
end

end # RomeoDFT
