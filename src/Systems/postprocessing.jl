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
    total_energy   = haskey(res, :total_energy) ? res[:total_energy][end] : typemax(Float64)
    accuracy   = haskey(res, :accuracy) ? res[:accuracy][end] : typemax(Float64)
    scf_iterations = haskey(res, :scf_iteration) ? res[:scf_iteration][end] : 0
    fermi          = haskey(res, :fermi) ? res[:fermi] : 0.0
    converged = res[:converged]
    return Results(state, minid, mindist, total_energy, hub_energy, scf_iterations,
                                       converged, fermi, accuracy)
end

function hubbard_outputdata(j; calcs = map(x->x.name, j.calculations), kwargs...)
    parse_func(d, x, y) = d[:Hubbard_iterations] = d[:scf_iteration][end]
    dict = Dict([n => ["Hubbard_conv true" => parse_func] for n in calcs])
    return outputdata(j; calcs = calcs, extra_parse_funcs=dict, kwargs...)
end

function Overseer.update(::ResultsProcessor, m::AbstractLedger)
    @debugv 2 "[START] ResultsProcessor"
    reslock = ReentrantLock()

    if isempty(m[Unique])
        unique_e = Entity(m, Unique())
    else
        unique_e = entity(m[Unique], 1)
    end
    unique_c = m[Unique][1]
    base_c = m[BaseCase] 
    @sync for e in @entities_in(m, SimJob && TimingInfo && (Trial || Template) && !Results)
        Threads.@spawn if ispath(joinpath(e.local_dir, "scf.out"))
            curt = now()
            j = local_load(Job(e.local_dir))
            o = Dict()
            try
                o = hubbard_outputdata(j, calcs = ["scf"]) 
                if !isempty(o)
                    e.job.calculations[1].run = false
                    res = o["scf"]
                    results = results_from_output(res, e in base_c)
                    e.scf_time = haskey(res, :timing) ? Dates.tons(res[:timing][end].cpu) / 1000 : e.running
                    lock(reslock)
                    try
                        m[e] = results
                    catch err
                        m[e] = Error(err, stacktrace(catch_backtrace()))
                    finally
                        unlock(reslock)
                    end
                    if haskey(res, :bands) && haskey(res, :total_magnetization)
                        bands = flatbands(res)
                        lock(reslock)
                        try
                            m[e] = FlatBands(bands)
                        catch err
                            m[e] = Error(err, stacktrace(catch_backtrace()))
                        finally
                            unlock(reslock)
                        end
                    end
                    lock(reslock)
                    try
                        state = results.state
                        if results.converged && !isempty(state.occupations) && e in m[FlatBands]
                            if unique_c.bands && e in m[FlatBands]
                                ebands = m[FlatBands][e].bands
                                if !any(x -> x.e != e.e &&
                                             sum(abs, m[Results][x].state.magmoms .- state.magmoms) <= unique_c.thr &&
                                             sssp_distance(ebands, x[FlatBands].bands, results.fermi) <= unique_c.thr,
                                         @entities_in(m, FlatBands && Unique))
                                         
                                    m[Unique][e] = unique_e
                                end
                            elseif !any(x -> e.e != x.e && Euclidean()(x.state, state) <= unique_c.thr, @entities_in(m, Results && Unique))
                                m[Unique][e] = unique_e
                            end
                        end
                        if e ∉ m[Unique] && !isempty(state.occupations)
                            e in m[FlatBands] && e ∉ m[BaseCase] && pop!(m[FlatBands], e)
                            m[e] = Done(false)
                        end
                        
                    catch err
                        m[e] = Error(err, stacktrace(catch_backtrace()))
                    finally
                        unlock(reslock)
                    end
                end
                e.postprocessing += Dates.datetime2unix(now()) - Dates.datetime2unix(curt)
            catch err
                lock(reslock)
                m[e] = Error(err, stacktrace(catch_backtrace()))
                unlock(reslock)
            end
        end
    end
    @debugv 2 "[STOP] ResultsProcessor"
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
Overseer.requested_components(::UniqueExplorer) = (Archived,Intersection, NSCFSettings, ProjwfcSettings, BandsSettings)

function Overseer.update(::UniqueExplorer, m::AbstractLedger)
    @debugv 2 "[START] UniqueExplorer"
    results = m[Results]
    generation = m[Generation]
    simulation = m[Simulation]
    nsettings = m[NSCFSettings]
    curgen = maximum(x->x.generation, m[Generation], init=0)
    new_states = 0
    for e in @entities_in(m, Results && SimJob && Unique)
        postprocess = false
        if !isempty(m[BandsSettings])
            band_e = entity(m[BandsSettings], length(m[BandsSettings]))
            m[BandsSettings][e] = band_e
            postprocess = true
        end
        if !isempty(m[NSCFSettings])
            nscf_e = entity(m[NSCFSettings], length(m[NSCFSettings]))
            m[NSCFSettings][e] = nscf_e
            postprocess = true
        end
        if !isempty(m[ProjwfcSettings])
            projwfc_e = entity(m[ProjwfcSettings], length(m[ProjwfcSettings]))
            m[ProjwfcSettings][e] = projwfc_e
            postprocess = true
        end
        if !postprocess
            m[e] = Done(false)
        end
        new_states += 1
    end
    simn = simname(m)
    if new_states != 0
        @debugv 1 "Found $new_states new unique states."
    end
    @debugv 2 "[STOP] UniqueExplorer"
end

