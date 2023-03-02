module v0_4
using ..Overseer
using ..v0_1
using ..v0_2
using ..v0_3
using ..DFControl: Projection
using ..Structures
using ..Calculations
import ..RomeoDFT: gencalc, State

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
end
TimingInfo() = TimingInfo(0, 0, 0, 0, 0, 0, "", 0.0)
function Base.convert(::Type{TimingInfo}, x::v0_1.TimingInfo)
    return TimingInfo(x.starttime, x.cur_time, x.pending, x.running, x.scf_time,
                      x.postprocessing, "", 0.0)
end

"""
    Error

Component that can be used to signal an error has occurred.
"""
@component struct Error
    err::Exception
    stack::String
end
Error(exc::Exception, trace::StackTraces.StackTrace) = Error(exc, join(trace, "\n"))
Error(str::String) = Error(ErrorException(str), stacktrace()[2:end])
Base.convert(::Type{Error}, x::v0_3.Error) = Error(ErrorException(x.msg), "")
function Base.show(io::IO, err::Error)
    (print(io, "Error: "); showerror(io, err.err); print(io, "\nTrace:\n"); println(io,
                                                                                    err.stack))
end

"""
    BandsSettings

Component that holds the settings to generate the bands calculation input and eventual band plot.
"""
@component Base.@kwdef mutable struct BandsSettings
    kpoints::Union{Vector{NTuple{4,Float64}},Int} = 20
    ymin::Float64 = -5.0
    ymax::Float64 = 5.0
end
Base.convert(::Type{BandsSettings}, x::v0_2.BandsSettings) =
    BandsSettings(x.kpoints, x.ymin, x.ymax)
    
function gencalc(job, settings::BandsSettings)
    if settings.kpoints isa Int
        kpath = filter(x -> all(y -> y in (1 / 3, 0.5, 0.0), x[1:3]),
                       Structures.high_symmetry_kpath(job.structure,
                                                      settings.kpoints))
        kpath[end] = (kpath[end][1:3]..., 1.0)
        return Calculations.gencalc_bands(job["scf"], kpath)
    else
        return Calculations.gencalc_bands(job["scf"], kpoints)
    end
end

"""
    NSCFSettings

Component that holds settings for the NSCF calculation.
"""
@component Base.@kwdef mutable struct NSCFSettings
    kpoints::NTuple{3,Int} = (6, 6, 6)
end
Base.convert(::Type{NSCFSettings}, x::v0_2.NSCFSettings) =
    NSCFSettings(x.kpoints)
gencalc(job, settings::NSCFSettings) = 
    Calculations.gencalc_nscf(job[1], settings.kpoints)

"""
    ProjwfcSettings

Component that holds settings for the projwfc calculation.
"""
@component Base.@kwdef mutable struct ProjwfcSettings
    Emin::Float64 = 20
    Emax::Float64 = 10
    deltaE::Float64 = 0.1
    dos_ratio::Float64 = 0.2
end
Base.convert(::Type{ProjwfcSettings}, x::v0_2.ProjwfcSettings) =
    ProjwfcSettings(x.Emin_rel, x.Emax_rel, x.deltaE, x.dosratio)
gencalc(job, settings::ProjwfcSettings) = 
    Calculations.gencalc_projwfc(job[1], settings.Emin, settings.Emax, settings.deltaE)
 

"""
    WannierSettings

Component that holds settings for the wannier calculation.
"""
@component Base.@kwdef mutable struct WannierSettings
    plot_wannier::Bool = false
    projections::Dict{Symbol,Vector{Projection}}
    dos_ratio::Float64 = 0.2
end
Base.convert(::Type{WannierSettings}, x::v0_2.WannierSettings) =
    WannierSettings(x.plot_wannier, x.projections, x.dos_ratio)
    
function gencalc(job, wsettings::WannierSettings)
    for (elsym, projs) in wsettings.projections
        for a in job.structure[elsym]
            a.projections = projs
        end
    end
    return Calculations.gencalc_wan(job, wsettings.dos_ratio; plot_wannier = wsettings.plot_wannier)
end 

@component struct Hybrid end

"""
    Trial
Holds the trial [`State`](@ref), i.e. where the system will be constrained towards and the origin of the state i.e. from random trials or other generating method.
"""

@component struct Trial
    state::State
    origin::String
end

Trial(s::State) = Trial(s, "unknown")
Base.convert(::Type{Trial}, x::v0_1.Trial) = Trial(x.state, "unknown")

function Base.show(io::IO, t::Trial)
    println(io, t.state)
    println(io, "origin: $(t.origin)")
end
    
end
