module v0_5
using ..Overseer
using ..v0_4
using ..RomeoDFT: State, Unknown, TrialOrigin, MixingMode, UnknownMixing, Structure, Calculation, DFWannier
using ..DFControl: Projection
using ..Structures
using ..Calculations
using ..v0_1
using ..v0_2
import ..RomeoDFT: gencalc, State

"""
    Trial
Holds the trial [`State`](@ref), i.e. where the system will be constrained towards and the origin of the state i.e. from random trials or other generating method.
"""

@component struct Trial
    state::State{Float64, DFWannier.MagneticVector{Float64, Vector{Float64}}, DFWannier.ColinMatrix{Float64, Matrix{Float64}}}
    origin::TrialOrigin
end

Trial(s::State) = Trial(s, Unknown)
Base.convert(::Type{Trial}, x::v0_4.Trial) = Trial(x.state, Unknown)

"""
    BaseCase

Signals the calculation that represents the base \"vanilla\" DFT + U case that is used to determine
for example the number of electrons on the magnetic ions.
"""
@component struct BaseCase end

"""
    TimingInfo

Stores some timing information.
"""
@component mutable struct TimingInfo
    starttime::Float64
    cur_time::Float64
    pending::Float64
    running::Float64
    scf_time::Float64
    postprocessing::Float64
    current_running::String
    current_runtime::Float64
    current_filesize::Float64
end
TimingInfo() = TimingInfo(0, 0, 0, 0, 0, 0, "", 0.0, 0.0)
function Base.convert(::Type{TimingInfo}, x::v0_4.TimingInfo)
    return TimingInfo(x.starttime, x.cur_time, x.pending, x.running, x.scf_time,
                      x.postprocessing, x.current_running, x.current_runtime, 0.0)
end

raw"""
    Simulation

The parameters of the FireFly algorithm.
Update is done through formula (see https://pubs.rsc.org/en/content/articlelanding/2019/cp/c9cp03618k):
```math
x_i^{t+1} = x_i^{t} \sum\limits_j \beta e^{-\gamma r_{ij}^2} \left( x_j^t - x_i^t) + \alpha_t \varepsilon_t
```
"""
@pooled_component mutable struct Simulation
    template_structure::Structure
    template_calculation::Calculation
    current_generation::Int
    current_best::Entity
    n_fireflies::Int
    mixing_mode::MixingMode
    
    #FireFly parameters
    γ::Float64
    α::Float64
    β::Float64
end
function Base.convert(::Type{Simulation}, x::v0_1.Simulation)
    return Simulation(x.template_structure, x.template_calculation, x.current_generation, x.current_best,
                      x.n_fireflies, UnknownMixing, x.γ, x.α, x.β)
end

"""
    ServerInfo

Component that holds where a job will run and which environments and executables to use.
"""
@pooled_component Base.@kwdef mutable struct ServerInfo
    server::String
    pw_exec::String
    environment::String
    priority::Int = 5
end
function Base.convert(::Type{ServerInfo}, x::v0_2.ServerInfo)
    return ServerInfo(server = x.server, pw_exec=x.pw_exec, environment=x.environment)
end

"""
    IntersectionSearcher

Searcher that is part of a search over all intersections (midpoints) between previously found metastable states.
"""
@pooled_component Base.@kwdef mutable struct IntersectionSearcher
    mindist::Float64 = 0.25    
end


"""
    Intersection

Intersection between two previously found metastable states.
"""
@component struct Intersection
    parent1::Entity
    parent2::Entity
end

"Identifies whether a entity is a unique state."
@pooled_component Base.@kwdef mutable struct Unique
    thr::Float64 = 1e-2
    bands::Bool = false
end

"Identifies whether a entity is a unique state."
@pooled_component Base.@kwdef struct RandomSearcher
    nsearchers::Int=10
end

"Identifies whether a entity is a unique state."
@component struct Done
    cleaned::Bool
end


"""
    BandsSettings

Component that holds the settings to generate the bands calculation input and eventual band plot.
"""
@pooled_component Base.@kwdef mutable struct BandsSettings
    kpoints::Union{Vector{NTuple{4,Float64}},Int} = 20
    ymin::Float64 = -5.0
    ymax::Float64 = 5.0
end
Base.convert(::Type{BandsSettings}, x::v0_4.BandsSettings) =
    BandsSettings(x.kpoints, x.ymin, x.ymax)
    

"""
    WannierSettings

Component that holds settings for the wannier calculation.
"""
@pooled_component Base.@kwdef mutable struct WannierSettings
    plot_wannier::Bool = false
    projections::Dict{Symbol,Vector{Projection}}
    dos_ratio::Float64 = 0.2
end
Base.convert(::Type{WannierSettings}, x::v0_4.WannierSettings) =
    WannierSettings(x.plot_wannier, x.projections, x.dos_ratio)
    
    
"""
    ProjwfcSettings

Component that holds settings for the projwfc calculation.
"""
@pooled_component Base.@kwdef mutable struct ProjwfcSettings
    Emin::Float64 = 20
    Emax::Float64 = 10
    deltaE::Float64 = 0.1
    dos_ratio::Float64 = 0.2
end
Base.convert(::Type{ProjwfcSettings}, x::v0_4.ProjwfcSettings) =
    ProjwfcSettings(x.Emin, x.Emax, x.deltaE, x.dos_ratio)

"""
    NSCFSettings

Component that holds settings for the NSCF calculation.
"""
@pooled_component Base.@kwdef mutable struct NSCFSettings
    kpoints::NTuple{3,Int} = (6, 6, 6)
end
Base.convert(::Type{NSCFSettings}, x::v0_4.NSCFSettings) =
    NSCFSettings(x.kpoints)

end
