function flatbands(all_out)
    first = all_out[:total_magnetization][end] > 0 ? :up : :down
    second = first == :up ? :down : :up
    bands = all_out[:bands]
    outbands = Vector{Float64}(undef,
                               2 * length(bands[first]) * length(bands[first][1].eigvals))
    icur = 1
    @inbounds for ib in 1:length(bands[:up])
        for s in (first, second)
            b = bands[s][ib]
            for e in b.eigvals
                outbands[icur] = e
                icur += 1
            end
        end
    end
    return outbands
end

"""
    sssp_distance(bands1, bands2, fermi)

Calculates
```math
\\sqrt{\\frac{\\sum_{mk} f_{mk} (\\varepsilon^1_{mk} - \\varepsilon^2_{mk} + \\Delta)^2}{\\sum_{mk}f_{mk}}}
```
for `bands1` and `bands2` flattened bandstructures. `fermi` is used to crudely determine ``f_{mk}`` and ``\\Delta`` is optimized to minimize the distance. 
"""
function sssp_distance(bands1, bands2, fermi)
    function obj(Δ)
        n = 0
        d = 0.0
        for (b1, b2) in zip(bands1, bands2)
            if b1 <= fermi && b2 <= fermi
                d += (b1 - b2 + Δ)^2
                n += 1
            end
        end
        return sqrt(d / n)
    end
    return optimize(x -> obj(x[1]), [0.0]).minimum
end

function process_Hubbard(hubbard, target=nothing)
    nhubbard = length(hubbard)
    states = Vector{State{Float64, DFWannier.MagneticVector{Float64, Vector{Float64}}, DFWannier.ColinMatrix{Float64, Matrix{Float64}}}}(undef, nhubbard)
    Threads.@threads for i in 1:nhubbard
        states[i] = State(hubbard[i])
    end
    target = target === nothing ? states[1] : target
    tstates = filter(x -> length(x) > 0, states)
    if length(tstates) > 5
        dists = map(x -> Euclidean()(target, x), tstates[5:end])

        mindist, minid = findmin(dists)
        minid += 5
    else
        dists = Float64[]
        mindist = typemax(Float64)
        minid = 0
    end

    state = states[end]
    (; dists, mindist, minid, state)
end

"""
    ResultsProcessor

Looks through entities with completed [`SimJobs`](@ref SimJob) and processes and stores the results. If the whole workflow of calculations for a
particular entity was finished it also cleans up the serverside directory.
It creates [`Result`](@ref) and [`BandsResults`](@ref) components for the correct entities.
"""
struct ResultsProcessor <: System end
Overseer.requested_components(::ResultsProcessor) = (Error, Unique, FlatBands)

function results_from_output(res::Dict, basecase=false)
    if haskey(res, :Hubbard) && !isempty(res[:Hubbard][1])
        try
            out = process_Hubbard(res[:Hubbard])
            if !haskey(res, :Hubbard_iterations) && !basecase
                minid = res[:scf_iteration][end]
                mindist = typemax(Float64)
                state = State()
            else
                minid = get(res, :Hubbard_iterations, out.minid)
                mindist = out.mindist
                state = out.state
            end
        catch 
            mindist = typemax(Float64)
            minid   = 0
            state   = State()
        end
    else
        mindist = typemax(Float64)
        minid   = 0
        state   = State()
    end
    hub_energy     = haskey(res, :Hubbard_energy) ? res[:Hubbard_energy][end] : typemax(Float64)
    total_energy   = haskey(res, :total_energy)   ? res[:total_energy][end]   : typemax(Float64)
    accuracy       = haskey(res, :accuracy)       ? res[:accuracy][end]       : typemax(Float64)
    scf_iterations = haskey(res, :scf_iteration)  ? res[:scf_iteration][end]  : 0
    fermi          = haskey(res, :fermi)          ? res[:fermi]               : 0.0
    converged      = res[:converged]
    return Results(state, minid, mindist, total_energy, hub_energy, scf_iterations,
                                       converged, fermi, accuracy)
end

function hubbard_outputdata(j; calcs = map(x->x.name, j.calculations), kwargs...)
    hub_parse_func(d, x, y) = d[:Hubbard_iterations] = d[:scf_iteration][end]
    
    parse_funcs = Dict([n => ["Hubbard_conv true" => hub_parse_func] for n in calcs])
    return outputdata(j; calcs = calcs, extra_parse_funcs=parse_funcs, kwargs...)
end

function Overseer.update(::ResultsProcessor, m::AbstractLedger)
    @error_capturing_threaded for e in @safe_entities_in(m, Pulled && SimJob && TimingInfo)
        if !ispath(joinpath(e.local_dir, "scf.out"))
            continue
        end
        
        curt = now()
        
        j = local_load(Job(e.local_dir))
        o = hubbard_outputdata(j, calcs = ["scf"])
        
        if !isempty(o)
            res = o["scf"]
            results = results_from_output(res, oldest_parent(m, e) in m[BaseCase])
            
            m[e] = results
            
            e.scf_time = haskey(res, :timing) ? Dates.tons(res[:timing][end].cpu) / 1000 : e.running
            
            if haskey(res, :bands) && haskey(res, :total_magnetization)

                bands = flatbands(res)
                m[e] = FlatBands(bands)

                fermi = res[:fermi]
                
                up   = res[:bands].up
                down = res[:bands].down

                n_up_conduction   = count(x->maximum(x.eigvals) > fermi, up)
                n_down_conduction = count(x->maximum(x.eigvals) > fermi, down)

                if n_up_conduction == 0 || n_down_conduction == 0
                    log(e, "There were no conduction bands for a spin channel, increasing nbnd and rerunning.")
                    suppress() do
                        new_template = deepcopy(m[Template][e])
                        c = new_template.calculation
                        if haskey(c, :nbnd)
                            c[:nbnd] = ceil(Int, 1.2 * c[:nbnd])
                        else
                            c[:nbnd] = ceil(Int, res[:n_KS_states] * 1.2)
                        end
                        should_rerun(m, e, new_template)
                    end
                end
            end
            
        end
        e.postprocessing += Dates.datetime2unix(now()) - Dates.datetime2unix(curt)
    end
end

struct BandsPlotter <: System end
Overseer.update(::BandsPlotter, ::AbstractLedger) = nothing

"""
    UniqueExplorer

Takes states newly found by the [`FireFly`](@ref) simulation that are unique and generates the [`SCFSettings`](@ref), [`BandsSettings`](@ref),
[`NSCFSettings`](@ref) and [`ProjwfcSettings`](@ref) for the post processing workflow.
After this, the flies themselves are [`Archived`](@ref).
"""
struct UniqueExplorer <: System end
Overseer.requested_components(::UniqueExplorer) = (Archived, Intersection, NSCFSettings, ProjwfcSettings, BandsSettings, HPSettings, Child)

function Overseer.update(::UniqueExplorer, m::AbstractLedger)
    if isempty(m[Unique])
        return
    end
    new_states = 0
    # First unique entity holds all the settings that are associated
    # with all the states found to be unique, i.e. relaxing projwfc etc
    unique_e = entity(m[Unique], 1)
    unique_c = m[Unique][unique_e]
    postprocess = any(x -> unique_e in m[x], (BandsSettings, NSCFSettings, ProjwfcSettings, RelaxSettings, HPSettings))
    postprocess_children = any(x -> unique_e in m[x], (RelaxSettings, HPSettings)) # because structure will change
    
    rescomp   = m[Results]
    bandscomp = m[FlatBands]

    @error_capturing for e in @safe_entities_in(m, Pulled && Results && !Unique && !Parent && !Done && !ShouldRerun)

        if !e.converged || isempty(e.state.occupations)
            m[e] = Done(false)
            continue
        end
        
        found = false
        ebands = bandscomp[e].bands
        @sync for e2 in @safe_entities_in(m, Unique && FlatBands && Results)
            Threads.@spawn begin
                found && return
                
                momdiffs = e2.state.magmoms .- e.state.magmoms
                sum(abs, momdiffs) >= unique_c.thr && return
                
                if unique_c.bands
                    if sssp_distance(ebands, e2[FlatBands].bands, e.fermi) >= unique_c.thr
                        return
                    end
                elseif Euclidean()(e2.state, e.state) >= unique_c.thr
                    return
                end
                found = true
            end
        end
        
        if found
            e ∉ m[BaseCase] && pop!(bandscomp, e)
            m[e] = Done(false)
            continue
        end
        
        m[Unique][e] = unique_e
        
        if !postprocess
            m[e] = Done(false)
            continue
        end
        
        if postprocess_children
            
            pp_e = create_postprocess_child!(m, e)
            m[e] = Done(false)
        else
            pp_e = e
        end
        
        new_states += 1
    end
    
    if new_states != 0
        @debugv 1 "Found $new_states new unique states."
    end
end

abstract type PostProcessSystem <: System end

struct BandsSystem <: PostProcessSystem end
struct NSCFSystem <: PostProcessSystem end
struct ProjwfcSystem <: PostProcessSystem end
struct PPSystem <: PostProcessSystem end

calcname(::BandsSystem)   = "bands"
calcname(::NSCFSystem)    = "nscf"
calcname(::ProjwfcSystem) = "projwfc"
calcname(::PPSystem)      = "pp"

Settings_T(::BandsSystem)   = BandsSettings
Settings_T(::NSCFSystem)    = NSCFSettings
Settings_T(::ProjwfcSystem) = ProjwfcSettings
Settings_T(::PPSystem)      = PPSettings
    
function Overseer.update(s::PostProcessSystem, m::AbstractLedger)
    cname = calcname(s)
    sT = Settings_T(s)
    @error_capturing for e in @safe_entities_in(m, SimJob && sT)
        if any(x->x.name == cname, e.job.calculations)
            continue
        end
        suppress() do
            add_calc!(e.job, e[sT])
        end
    end
    @error_capturing for e in @safe_entities_in(m, SimJob && Pulled)
        if !ispath(joinpath(e.local_dir, "$cname.out"))
            continue
        end
        if e.job.calculations[end].name == cname
            m[e] = Done(false)
        end
    end
end
