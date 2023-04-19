module v0_1
using ..Overseer
using ..RomeoDFT
using ..RomeoDFT: State, DFWannier

@component mutable struct Timer
    prev_t::Float64
    dt::Float64
end

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
end
TimingInfo() = TimingInfo(0, 0, 0, 0, 0, 0)

@component mutable struct SimJob
    local_dir::String
    remote_dir::String
    submitted::Bool
end

"""
    Trial

Holds the trial [`State`](@ref), i.e. where the system will be constrained towards. 
"""
@component struct Trial
    state::State
end

@component struct Results
    state::State{Float64, DFWannier.MagneticVector{Float64, Vector{Float64}}, DFWannier.ColinMatrix{Float64, Matrix{Float64}}}
    constraining_steps::Int
    closest_to_target::Float64
    total_energy::Float64
    Hubbard_energy::Float64
    niterations::Int
    converged::Bool
end

@pooled_component struct ServerInfo
    server::String
    pw_exec::String
    environment::String
end

"""
    Generation

Stores which generation of the search an entity belongs to.
"""
@component struct Generation
    generation::Int
end

"""
    Simulation

The parameters of the FireFly algorithm.
"""
@pooled_component mutable struct Simulation
    template_structure::Structure
    template_calculation::Calculation
    current_generation::Int
    current_best::Entity
    n_fireflies::Int
    γ::Float64
    α::Float64
    β::Float64
end

# These mirror to some degree the RemoteHPC.JobState
"""
    Submit

Signals whether a job should be submitted.
"""
@component struct Submit    end
@component struct Submitted end
@component struct Running   end
@component struct Completed end
@component struct Pulled end

@pooled_component Base.@kwdef mutable struct RelaxSettings
    force_convergence_threshold::Float64  = 1e-3
    energy_convergence_threshold::Float64 = 1e-4
    ion_dynamics::String                  = "bfgs"
    cell_dynamics::String                 = "bfgs"
    symmetry::Bool                        = true
    variable_cell::Bool                   = true
end

"""
    RelaxResults

Holds the results of a relaxation.
"""
@component struct RelaxResults
    n_steps::Int
    total_force::Float64
    final_structure::Structure
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
Base.push!(l::Log, args...) = push!(l.logs, args...)
Base.prepend!(l::Log, args...) = prepend!(l.logs, args...)
Base.append!(l::Log, args...) = append!(l.logs, args...)
Base.push!(l::Log, l2::Log) = push!(l.logs, l2.logs)
Base.prepend!(l::Log, l2::Log) = prepend!(l.logs, l2.logs)
Base.append!(l::Log, l2::Log) = append!(l.logs, l2.logs)

"""
    ShouldRerun

Signals that a job should be reran.
`data_to_pop` denotes the [`Components`](@ref Component) that should be removed from the [`Entity`](@ref).
"""
@component struct ShouldRerun
    data_to_pop::Set{DataType}
end
ShouldRerun(args::DataType...) = ShouldRerun(Set(args))
Base.push!(s::ShouldRerun, d::DataType) = push!(s.data_to_pop, d)

"""
    Rerun

Saves how many times a job has been reran.
"""
@component struct Rerun
    count::Int
end

"""
    HPSettings

Holds the settings for a HP calculation.
"""
@pooled_component Base.@kwdef struct HPSettings
    nq::NTuple{3, Int} = (2,2,2)
    conv_thr_chi::Float64 = 1e-6
    find_atpert::Int = 1
    U_conv_thr::Float64 = 0.1
end

"""
    HPResults

Holds the results of a HP calculations.
"""
@component struct HPResults
    U::Vector
end

"""
    Child

Signals that an [`Entity`](@ref) has a child entity, e.g. as part of postprocessing.
"""
@component struct Child
    child::Entity
end

"""
    Parent

Signals that an [`Entity`](@ref) has a Parent from which it's derived.
"""
@component struct Parent
    parent::Entity
end

@component struct RelaxChild
    child::Entity
end

end
