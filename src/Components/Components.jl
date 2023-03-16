abstract type PostProcessSettings end

const COMP_FILES = filter(x->x != "Components.jl", readdir(@__DIR__))

const DATABASE_VERSION = maximum(x -> VersionNumber(splitext(x)[1]), COMP_FILES) 
const VERSION_MODS = map(COMP_FILES) do x
    return VersionNumber(splitext(x)[1]) => include(joinpath(@__DIR__, x))
end

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
const Results              = @eval getfield(TYPE_MODS[:Results][end], :Results)
const RelaxSettings        = @eval getfield(TYPE_MODS[:RelaxSettings][end], :RelaxSettings)
const Results              = @eval getfield(TYPE_MODS[:Results][end], :Results)
const SCFSettings          = @eval getfield(TYPE_MODS[:SCFSettings][end], :SCFSettings)
const Intersection         = @eval getfield(TYPE_MODS[:Intersection][end], :Intersection)
const Error                = @eval getfield(TYPE_MODS[:Error][end], :Error)
const RelaxResults         = @eval getfield(TYPE_MODS[:RelaxResults][end], :RelaxResults)
const Completed            = @eval getfield(TYPE_MODS[:Completed][end], :Completed)
const Pulled               = @eval getfield(TYPE_MODS[:Pulled][end], :Pulled)
const FlatBands            = @eval getfield(TYPE_MODS[:FlatBands][end], :FlatBands)
const Hybrid               = @eval getfield(TYPE_MODS[:Hybrid][end], :Hybrid)
const Archived             = @eval getfield(TYPE_MODS[:Archived][end], :Archived)
const IntersectionSearcher = @eval getfield(TYPE_MODS[:IntersectionSearcher][end], :IntersectionSearcher)
const Running              = @eval getfield(TYPE_MODS[:Running][end], :Running)
const SimJob               = @eval getfield(TYPE_MODS[:SimJob][end], :SimJob)
const Trial                = @eval getfield(TYPE_MODS[:Trial][end], :Trial)
const Timer                = @eval getfield(TYPE_MODS[:Timer][end], :Timer)
const Template             = @eval getfield(TYPE_MODS[:Template][end], :Template)
const WannierResults       = @eval getfield(TYPE_MODS[:WannierResults][end], :WannierResults)
const Parent               = @eval getfield(TYPE_MODS[:Parent][end], :Parent)
const Simulation           = @eval getfield(TYPE_MODS[:Simulation][end], :Simulation)
const StopCondition        = @eval getfield(TYPE_MODS[:StopCondition][end], :StopCondition)
const Done                 = @eval getfield(TYPE_MODS[:Done][end], :Done)
const ServerInfo           = @eval getfield(TYPE_MODS[:ServerInfo][end], :ServerInfo)
const Unique               = @eval getfield(TYPE_MODS[:Unique][end], :Unique)
const Submitted            = @eval getfield(TYPE_MODS[:Submitted][end], :Submitted)
const WannierSettings      = @eval getfield(TYPE_MODS[:WannierSettings][end], :WannierSettings)
const BaseCase             = @eval getfield(TYPE_MODS[:BaseCase][end], :BaseCase)
const BandsResults         = @eval getfield(TYPE_MODS[:BandsResults][end], :BandsResults)
const Generation           = @eval getfield(TYPE_MODS[:Generation][end], :Generation)
const ProjwfcSettings      = @eval getfield(TYPE_MODS[:ProjwfcSettings][end], :ProjwfcSettings)
const RandomSearcher       = @eval getfield(TYPE_MODS[:RandomSearcher][end], :RandomSearcher)
const BandsSettings        = @eval getfield(TYPE_MODS[:BandsSettings][end], :BandsSettings)
const TimingInfo           = @eval getfield(TYPE_MODS[:TimingInfo][end], :TimingInfo)
const NSCFSettings         = @eval getfield(TYPE_MODS[:NSCFSettings][end], :NSCFSettings)
const Log                  = @eval getfield(TYPE_MODS[:Log][end], :Log)
const Submit               = @eval getfield(TYPE_MODS[:Submit][end], :Submit)
const ShouldRerun          = @eval getfield(TYPE_MODS[:ShouldRerun][end], :ShouldRerun)
const Rerun                = @eval getfield(TYPE_MODS[:Rerun][end], :Rerun)
const HPSettings           = @eval getfield(TYPE_MODS[:HPSettings][end], :HPSettings)
const HPResults            = @eval getfield(TYPE_MODS[:HPResults][end], :HPResults)
const Child                = @eval getfield(TYPE_MODS[:Child][end], :Child)
const Parent               = @eval getfield(TYPE_MODS[:Parent][end], :Parent)
const RelaxChild           = @eval getfield(TYPE_MODS[:RelaxChild][end], :RelaxChild)
const Bin           = @eval getfield(TYPE_MODS[:Bin][end], :Bin)

for (t, mod) in TYPE_MODS
    mt = getfield(mod[end], t)
    v = VersionNumber(replace(split(string(mod[end]), ".")[end], "_" => "."))
    p = joinpath(@__DIR__, string(v) * ".jl")
    @eval export $t
    @eval Base.pathof(::Type{$t}) = $p
end
