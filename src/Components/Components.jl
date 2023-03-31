abstract type PostProcessSettings end

const COMP_FILES = filter(x->x != "Components.jl", readdir(@__DIR__))

const DATABASE_VERSION = maximum(x -> VersionNumber(splitext(x)[1]), COMP_FILES)
include("0.1.0.jl")
include("0.2.0.jl")
include("0.3.0.jl")
include("0.4.0.jl")
include("0.5.0.jl")
include("0.6.0.jl")

const VERSION_MODS = [v"0.1.0" => v0_1,
                      v"0.2.0" => v0_2,
                      v"0.3.0" => v0_3,
                      v"0.4.0" => v0_4,
                      v"0.5.0" => v0_5,
                      v"0.6.0" => v0_6]

versions() = first.(VERSION_MODS)

function type_mods(first::VersionNumber = VersionNumber(0, 1),
                   last::VersionNumber = DATABASE_VERSION)
    out = Dict{Symbol,Vector{Module}}()
    for (v, mod) in filter(x -> x[1] <= last, VERSION_MODS)
        types = filter(x -> (t = getfield(mod, x); t isa DataType && string(x)[1] != '#'),
                       names(mod; all = true))
        for t in types
            if !haskey(out, t)
                out[t] = [mod]
            elseif v <= first
                out[t][1] = mod
            else
                push!(out[t], mod)
            end
        end
    end
    return out
end

# Register current type versions
const TYPE_MODS = type_mods()

const Template             = RomeoDFT.v0_2.Template
const WannierResults       = RomeoDFT.v0_2.WannierResults
const Done                 = RomeoDFT.v0_5.Done
const Child                = RomeoDFT.v0_1.Child
const Submit               = RomeoDFT.v0_1.Submit
const Results              = RomeoDFT.v0_6.Results
const SCFSettings          = RomeoDFT.v0_3.SCFSettings
const Rerun                = RomeoDFT.v0_1.Rerun
const Completed            = RomeoDFT.v0_1.Completed
const Pulled               = RomeoDFT.v0_1.Pulled
const IntersectionSearcher = RomeoDFT.v0_6.IntersectionSearcher
const StopCondition        = RomeoDFT.v0_6.StopCondition
const SimJob               = RomeoDFT.v0_3.SimJob
const Bin                  = RomeoDFT.v0_6.Bin
const Simulation           = RomeoDFT.v0_5.Simulation
const Unique               = RomeoDFT.v0_5.Unique
const Submitted            = RomeoDFT.v0_1.Submitted
const BaseCase             = RomeoDFT.v0_5.BaseCase
const BandsResults         = RomeoDFT.v0_2.BandsResults
const Generation           = RomeoDFT.v0_1.Generation
const BandsSettings        = RomeoDFT.v0_5.BandsSettings
const TimingInfo           = RomeoDFT.v0_5.TimingInfo
const Log                  = RomeoDFT.v0_1.Log
const RelaxSettings        = RomeoDFT.v0_1.RelaxSettings
const Error                = RomeoDFT.v0_4.Error
const Hybrid               = RomeoDFT.v0_4.Hybrid
const Archived             = RomeoDFT.v0_2.Archived
const Trial                = RomeoDFT.v0_5.Trial
const Parent               = RomeoDFT.v0_1.Parent
const ProjwfcSettings      = RomeoDFT.v0_5.ProjwfcSettings
const RelaxChild           = RomeoDFT.v0_1.RelaxChild
const NSCFSettings         = RomeoDFT.v0_5.NSCFSettings
const Intersection         = RomeoDFT.v0_5.Intersection
const RelaxResults         = RomeoDFT.v0_2.RelaxResults
const HPResults            = RomeoDFT.v0_1.HPResults
const FlatBands            = RomeoDFT.v0_3.FlatBands
const Running              = RomeoDFT.v0_1.Running
const Timer                = RomeoDFT.v0_1.Timer
const ServerInfo           = RomeoDFT.v0_5.ServerInfo
const ShouldRerun          = RomeoDFT.v0_1.ShouldRerun
const WannierSettings      = RomeoDFT.v0_5.WannierSettings
const RandomSearcher       = RomeoDFT.v0_5.RandomSearcher
const HPSettings           = RomeoDFT.v0_1.HPSettings

export Template            
export WannierResults      
export Done                
export Child               
export Submit              
export Results             
export SCFSettings         
export Rerun               
export Completed           
export Pulled              
export IntersectionSearcher
export StopCondition       
export SimJob              
export Bin                 
export Simulation          
export Unique              
export Submitted           
export BaseCase            
export BandsResults        
export Generation          
export BandsSettings       
export TimingInfo          
export Log                 
export RelaxSettings       
export Error               
export Hybrid              
export Archived            
export Trial               
export Parent              
export ProjwfcSettings     
export RelaxChild          
export NSCFSettings        
export Intersection        
export RelaxResults        
export HPResults           
export FlatBands           
export Running             
export Timer               
export ServerInfo          
export ShouldRerun         
export WannierSettings     
export RandomSearcher      
export HPSettings          
