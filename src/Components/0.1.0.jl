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

@component struct Parent
    parent::Entity
end
@component struct RelaxChild
    child::Entity
end

@component struct RelaxResults
    n_steps::Int
    total_force::Float64
    final_structure::Structure
end    

end
