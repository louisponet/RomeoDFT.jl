using RemoteHPC
using DFControl: Calculation
const DATABASE_VERSION = VersionNumber(0, 6)

function version_convert(::Type{T}, x) where {T}
    return error("No version conversion function for type $T from $(typeof(x))")
end

const VERSION_MODS = map(readdir(joinpath(@__DIR__, "Components"))) do x
    return VersionNumber(splitext(x)[1]) => include(joinpath(@__DIR__, "Components", x))
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
for (t, mod) in type_mods()
    mt = getfield(mod[end], t)
    v = VersionNumber(replace(split(string(mod[end]), ".")[end], "_" => "."))
    p = joinpath(@__DIR__, "Components", string(v) * ".jl")
    @eval const $t = $mt
    @eval export $t
    @eval Base.pathof(::Type{$t}) = $p
end

function version_convert(v, type_mods, id)
    cur_type = getfield(last(type_mods)[id], first(type_mods))
    if id == length(last(type_mods))
        return convert(cur_type, v)
    else
        return version_convert(convert(cur_type, v), type_mods, id + 1)
    end
end

Base.convert(::Type{T}, s) where T<:State = State(s.occupations)
Base.convert(::Type{T}, s::T) where T<:State = State(s.occupations)

function RemoteHPC.load(rootdir::String, l::AbstractLedger)
    ledger_version = jldopen(joinpath(rootdir, "ledger.jld2"), "r") do f
        if haskey(f, "version")
            return f["version"]
        else
            return VersionNumber(0, 1)
        end
    end

    @assert ledger_version in versions() "Unknown version: $ledger_version"

    type2mod = type_mods()

    alltypes = DataType[]
    for (k, mods) in type2mod
        for m in mods
            push!(alltypes, getfield(m, k))
        end
    end
    typemap = Dict([replace(string(t), "RomeoDFT" => "Occupations") => t for t in alltypes])
    for t in keys(type2mod)
        typemap["Occupations." * string(t)] = getfield(RomeoDFT, t)
    end
    typemap["Occupations.TrialOrigin"] = RomeoDFT.TrialOrigin
    typemap["Occupations.MixingMode"] = RomeoDFT.MixingMode
    for sys in filter(x -> isdefined(RomeoDFT, x) && getfield(RomeoDFT, x) isa DataType && getfield(RomeoDFT, x) <: System, names(RomeoDFT, all=true))
        typemap["Occupations."*string(sys)] = getfield(RomeoDFT, sys)
    end
    
    ledger = jldopen(joinpath(rootdir, "ledger.jld2"), "r", typemap=typemap) do f
        l.sleep_time = get(f, "sleep_time", l.sleep_time)
        l.mode = get(f, "mode", l.mode)
        l.mode = get(f, "searcher_stage", l.mode)
        l.finished = get(f, "finished", l.finished)
        return f["ledger"]
    end
    for c in keys(components(l))
        Overseer.ensure_component!(ledger, c)
    end
    if l.mode == :searching
        l.mode = :search
    end
    set_searcher_stages!(ledger, l.mode)

    for (t, mods) in type2mod
        # do sequential update
        for i in 1:length(mods)-1
            old_t = getfield(mods[i], t)
            new_t = getfield(mods[i+1], t)
            if old_t != new_t && old_t ∈ ledger
                Overseer.ensure_component!(ledger, new_t)
                newcomp = ledger[new_t]
                for e in @entities_in(ledger, old_t)
                    newcomp[e] = version_convert(e[old_t], (t, mods[2:end]), 1)
                end
                delete!(components(ledger), old_t)
                if newcomp isa PooledComponent
                    Overseer.make_unique!(newcomp)
                end
            end
        end
    end
    l.ledger = ledger
    return l
end

function RemoteHPC.save(rootdir::String, l::AbstractLedger)
    lp = joinpath(rootdir, "ledger.jld2")
    if ispath(lp)
        cp(lp, joinpath(rootdir, "ledger_bak.jld2"); force = true)
    end
    return JLD2.jldsave(joinpath(rootdir, "ledger.jld2");
                        ledger     = Overseer.ledger(l),
                        version    = DATABASE_VERSION,
                        sleep_time = l.sleep_time,
                        finished   = l.finished,
                        mode       = l.mode)
end
