module v0_2
using ..Overseer
using ..DFControl.Jobs
using ..v0_1
using ..RomeoDFT: State, DFWannier, MagneticVectorType, ColinMatrixType
using ..RomeoDFT.DFControl: Projection, Calculation, QE, Structure

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
end
function Base.convert(::Type{Results}, x::v0_1.Results)
    return Results(x.state,
                   x.constraining_steps,
                   x.closest_to_target,
                   x.total_energy,
                   x.Hubbard_energy,
                   x.niterations,
                   x.converged,
                   0.0)
end

"""
    ServerInfo

Component that holds where a job will run and which environments and executables to use.
"""
@pooled_component struct ServerInfo
    server::String
    pw_exec::String
    environment::String
    single_node_environment::String
end
function Base.convert(::Type{ServerInfo}, x::v0_1.ServerInfo)
    return ServerInfo(x.server, x.pw_exec, x.environment,
                      x.server ∈ ("daint", "theospc22") ? "default" : "normal_1node")
end

@component mutable struct SimJob
    local_dir::String
    remote_dir::String
    jobstate::Jobs.JobState
end
# x is supposed to be SimJob of previous version
function Base.convert(::Type{SimJob}, x::v0_1.SimJob)
    return SimJob(replace(x.local_dir, "v0.1" => "0.1.0"), x.remote_dir, RemoteHPC.Unknown)
end

"""
    BandsResults

Stores the path to a bandstructure plot.
"""
@component struct BandsResults
    plot_file::String
end

"""
    BandsSettings

Component that holds the settings to generate the bands calculation input and eventual band plot.
"""
@component Base.@kwdef mutable struct BandsSettings
    kpoints::Union{Vector{NTuple{4,Float64}},Int} = 20
    ymin::Float64 = -5.0
    ymax::Float64 = 5.0
    created::Bool = false
end

"""
    NSCFSettings

Component that holds settings for the NSCF calculation.
"""
@component Base.@kwdef mutable struct NSCFSettings
    kpoints::NTuple{3,Int} = (6, 6, 6)
    created::Bool = false
end

"""
    ProjwfcSettings

Component that holds settings for the projwfc calculation.
"""
@component Base.@kwdef mutable struct ProjwfcSettings
    Emin_rel::Float64 = 20
    Emax_rel::Float64 = 10
    deltaE::Float64 = 0.1
    dosratio::Float64 = 0.2
    created::Bool = false
end

"""
    WannierSettings

Component that holds settings for the wannier calculation.
"""
@component Base.@kwdef mutable struct WannierSettings
    plot_wannier::Bool = false
    projections::Dict{Symbol,Vector{Projection}}
    dos_ratio::Float64 = 0.2
    created::Bool = false
end

"""
    WannierResults

Component that holds the important results of a wannier calculation.
"""
@component Base.@kwdef struct WannierResults
    hami_path::String
    wfuncs_path::String
end

"""
    Archived

Component that is used to archive entities that are not active.
"""
@component mutable struct Archived
    archive_path::String
    isarchived::Bool
end

"""
    Template

Component that holds the template structure and calculation that will be used during the search.
"""
@pooled_component struct Template
    structure::Structure
    calculation::Calculation{QE}
end

@component struct SCFSettings
    replacement_flags::Dict
end

"""
    RelaxResults

Holds the results of a relaxation.
"""
@component struct RelaxResults
    n_steps::Int
    total_force::Float64
    final_structure::Structure
    diff::Float64
end

Base.convert(::Type{RelaxResults}, x::v0_1.RelaxResults) = RelaxResults(x.n_steps, x.total_force, x.final_structure, 0.0)


end
