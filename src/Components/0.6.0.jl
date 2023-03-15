module v0_6
using ..Overseer
using ..v0_5
using ..RomeoDFT: MagneticVectorType, ColinMatrixType, State
using ..v0_2
"""
    IntersectionSearcher

Searcher that is part of a search over all intersections (midpoints) between previously found metastable states.
- `mindist`: minimum Euclidean distance to other trials for an intersection to be accepted.
- `max_intersections_per_generation`: maximum number of new trials for each generation. Most "distant" Trials will be used. 
"""
@pooled_component Base.@kwdef mutable struct IntersectionSearcher
    mindist::Float64 = 0.25
    max_intersections_per_generation::Int = 100
end
function Base.convert(::Type{IntersectionSearcher}, x::v0_5.IntersectionSearcher)
    return IntersectionSearcher(mindist=x.mindist)
end

"""
    StopCondition

This controls when to stop the search.
When the number of unique new states per trial is below the `unique_ratio` for `n_generations` consecutive generations,
the stop condition is met and the search will go in finalizing mode i.e. finish running trials and cleanup.
"""
Base.@kwdef mutable struct StopCondition
    unique_ratio::Float64 = 0.1
    n_generations::Int = 3
end

"""
    Results

Component holding the important results of SCF calculations.
"""
@component struct Results
    state::State{Float64, MagneticVectorType, ColinMatrixType}
    constraining_steps::Int
    closest_to_target::Float64
    total_energy::Float64
    Hubbard_energy::Float64
    niterations::Int
    converged::Bool
    fermi::Float64
    accuracy::Float64
end
function Base.convert(::Type{Results}, x::v0_2.Results)
    return Results(x.state, x.constraining_steps, x.closest_to_target, x.total_energy, x.Hubbard_energy, x.niterations, x.converged, x.fermi, x.converged ? 0.0 : typemax(Float64))
end

"""
    Log
Used for storing logs to an entity.
"""
@component struct Log
    logs::Vector{String}
end
Log() = Log(String[])

function Base.show(io::IO, log::Log)
    print(io, typeof(log), "($(length(log.logs)) logs)")
end
function Base.show(io::IO, ::MIME"text/plain", log::Log)
    println(io, typeof(log), " with $(length(log.logs)) logs:")
    for (i, l) in enumerate(log.logs)
        println(io, "[$i] $l")
    end
end


end
