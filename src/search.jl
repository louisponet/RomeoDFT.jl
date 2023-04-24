"""
    Searcher(name)

Represents a global search with label `name`. 

## Data
The `name` also references the sub directory where all the data will be stored, which can be retrieved using `dirname`.
This data consists of:
- `$(DATABASE_VERSION)/ledger.jld2`: stores all the search data that is kept in memory during the search. See [`load`](@ref) and [`save`](@ref)
- `job_backups`: a directory storing all the input and output data of the various calculations ran during the search
- `log.log`: file that keeps a log of everything that happens during the search. The amount of information logged can be tuned using different levels passed to `--verbosity <n>` flag of `romeo searcher create` or `romeo searcher start`
- `HTTP.log`: file holding all HTTP related logs
- `report.out`: a report generated when a [`Searcher`](@ref) is finished

## CommandLine
Most interaction with a [`Searcher`](@ref) can be done through the command line using the `romeo searcher` subcommands (use `-h` to show more info on commands and flags).
The most important ones are:
- `romeo searcher create`: interactive wizard to initialize and start a new [`Searcher`](@ref)
- `romeo searcher start`: starts a [`Searcher`](@ref)
- `romeo searcher stop`: stops a [`Searcher`](@ref)
- `romeo searcher status`: shows a status overview
- `romeo searcher unload`: unloads a [`Searcher`](@ref) from the orchestrator
- `romeo searcher log`: opens the `log.log` file in the editor specified by the `\$EDITOR` environment variable
- `romeo searcher plot`: creates a standard [`State`](@ref) plot
- `romeo searcher mode`: sets the mode, see [Modes](@ref) for more information

## Life Cycle
The lifecycle of a [`Searcher`](@ref) can be summarized by the following steps:
1. Potential geometry and U parameter optimization when `--hp-base` or `--relax-base` are passed to `romeo searcher create`
2. Initial unconstrained `scf` calculation to determine the number of electrons in the target manifold.
3. Random search with `--nrand` trials
4. Based on the unique [`States`](@ref State) found in `3.` (determined by `--unique-thr`) perform midpoint intersections (`n = (n1 + n2)/2`) between each pair of unique states
5. repeat `4.` for each newly found unique state
6. if no new intersections can be created and the stopcondition is not yet reached go to step `3.`
7. if the stop condition determined by `--stopping-unique-ratio` and `--stopping-n-generations` is reached switch to `postprocess` mode
8. if all jobs are finished and the stop condition is still satisfied, stop the [`Searcher`](@ref) and print results to `report.out`, otherwise switch back to `3.`

## Modes
There several execution modes that a [`Searcher`](@ref) can be in:
- `search`: search mode characterized by performing intersection and random search, and postprocessing on unique [`States`](@ref State).
- `postprocess`: when the stopcondition is satisfied, no new intersections or random trials will be generated, and only `postprocessing` actions will be taken
- `cleanup`: no new trials or postprocessing jobs will be created and submitted, only cleanup of the currently running jobs will be performed, after which the [`Searcher`](@ref) will finish
- `manual`: only used when performing manual tasks and you don't want any searching or postprocessing to happen automatically

These can be set using `romeo searcher mode`. 
"""
Base.@kwdef mutable struct Searcher <: AbstractLedger
    rootdir    ::String
    loop       ::Union{Nothing,Task}           = nothing
    stop       ::Bool                          = false
    finished   ::Bool                          = false
    sleep_time ::Float64                       = 30
    loop_error ::Bool                          = false
    mode       ::Symbol                        = :postprocess # What is the searcher doing i.e. searching or merely postprocessing
    ledger     ::Ledger                        = Ledger((mode == :postprocess ? [core_stage()] : search_stages())...)
    locks      ::Dict{DataType, ReentrantLock} = Dict{DataType, ReentrantLock}()
    timer      ::TimerOutput                   = TimerOutput()
    verbosity  ::Int                           = 0
end

mode(m::Searcher) = m.mode

function Searcher(dir::AbstractString; kwargs...)
    dir = isdir(dir) ? dir : searchers_dir(dir)
    Searcher(; rootdir = abspath(dir), kwargs...)
end

function Base.show(io::IO, l::Searcher)
    println(io, "Searcher:")
    println(io, "rootdir: $(l.rootdir)")
    if Simulation ∈ l
        count = 1
        for (s, es) in pools(l[Simulation])
            println(io, "current_best: ", s.current_best)
            println(io, "current_generation: ", s.current_generation)
            count += 1
        end
    end
    println(io, "unique states: $(length(l[SCFSettings]))")

    show(io, l.ledger)
    return 
end

function searcher_name(s::AbstractString)
    splt = strip.(split(s, searchers_dir()), '/')
    if length(splt) > 1
        return splt[end]
    else
        return s
    end
end
searcher_name(m::Searcher) = searcher_name(m.rootdir)

Base.dirname(m::Searcher) = m.rootdir
Base.joinpath(m::Searcher, p...)                         = joinpath(realpath(dirname(m)), p...)
Base.joinpath(s::Server, m::Searcher, p...)              = RemoteHPC.islocal(s) ? joinpath(m, p...) : abspath(s, joinpath("RomeoDFT", searcher_name(m), p...))
Base.joinpath(s::Server, m::Searcher, e::AbstractEntity) = joinpath(s, m, "$(Entity(e).id)")
Base.joinpath(m::Searcher, e::AbstractEntity, p...) = joinpath(m, "job_backups", "$(Entity(e).id)", p...)
local_dir(m::Searcher, e::AbstractEntity) = joinpath(m, e)

function simname(m)
    fdir = searcher_name(m)
    return replace(fdir, "/" => "_")
end

all_servers(m::Searcher) = [map(x -> Server(x.server), m[ServerInfo].data); local_server()]

function average_runtime(m::Searcher)
    successful = filter(x->x.converged, @entities_in(m, Results && TimingInfo))
    if length(successful) == 0
        return 0
    else
        return sum(x->x.scf_time, successful) / length(successful)
    end
end

##### OVERSEER Functionality
Overseer.ledger(l::Searcher) = l.ledger

function Base.getindex(l::Searcher, ::Type{T}) where {T}
    lck = get!(l.locks, T, ReentrantLock())
    SafeLoggingComponent(l.ledger[T], l.ledger[Log], lck)
end

function Overseer.Entity(l::Searcher, args...)
    return Entity(l.ledger, Log(), args...)
end

function Overseer.update(l::Searcher)
    LoggingExtras.withlevel(LoggingExtras.Debug; verbosity=l.verbosity) do
        @debugv 1 "[START] Updating" 
        for s in l.ledger.stages
            update(s, l)
        end
        @debugv 1 "[STOP] Updating"
    end
    l.timer = TimerOutputs.flatten(l.timer)
    return nothing
end

function Overseer.update(stage::Stage, l::Searcher)
    # Steps in a stage get executed in sequence, but if
    # a step is a vector they are threaded
    for step in stage.steps
        if step isa Vector
            Threads.@threads for t in step
                if t isa System
                    # to = TimerOutput()
                    @debugv 4 "[START] $t"
                    # @timeit to "$(typeof(t))" update(t, l)
                    update(t, l)
                    @debugv 4 "[STOP] $t"
                    # merge!(l.timer, to)
                else
                    update(t, l)
                end
            end
        else
            if step isa System
                # to = TimerOutput()
                @debugv 4 "[START] $step"
                # @timeit to "$(typeof(step))" update(step, l)
                update(step, l)
                @debugv 4 "[STOP] $step"
                # merge!(l.timer, to)
            else
                update(step, l)
            end
        end
    end
end

####### RemoteHPC
storage_dir(l::Searcher) = joinpath(l.rootdir, string(DATABASE_VERSION))

function RemoteHPC.save(l::Searcher)
    return (mkpath(storage_dir(l)); save(storage_dir(l), l))
end

fieldnames_to_save(::Searcher) = filter(x -> x ∉ (:rootdir, :loop, :ledger, :locks, :timer), fieldnames(Searcher))

function RemoteHPC.save(rootdir::String, l::Searcher)
    lp = joinpath(rootdir, "ledger.jld2")
    if ispath(lp)
        cp(lp, joinpath(rootdir, "ledger_bak.jld2"); force = true)
    end


    JLD2.jldopen(joinpath(rootdir, "ledger.jld2"), "w") do f
    
        comp_group = JLD2.Group(f, "components")
        for (T, c) in components(l)
            if !isempty(c)
                comp_group["$T"] = c
            end
        end
        
        for name in fieldnames_to_save(l)
            f["$name"] = getfield(l, name)
        end

        f["entities"] = entities(l)
        f["free_entities"] = Overseer.free_entities(l)
    end
    return joinpath(rootdir, "ledger.jld2")
end

function RemoteHPC.load(l::Searcher; version = nothing)
    if version === nothing
        versions = []
        for d in readdir(l.rootdir)
            
            !isdir(joinpath(l.rootdir, d)) || "ledger.jld2" ∉ readdir(joinpath(l.rootdir, d)) && continue
            
            try
                push!(versions, VersionNumber(d))
            catch
                nothing
            end
        end
        version = maximum(versions)
    end
    return load(joinpath(l.rootdir, string(version)), l)
end

function RemoteHPC.load(rootdir::String, l::AbstractLedger)
    ledger_version = JLD2.jldopen(joinpath(rootdir, "ledger.jld2"), "r") do f
        if haskey(f, "version")
            return f["version"]
        else
            return VersionNumber(0, 1)
        end
    end

    @assert ledger_version in versions() "Unknown version: $ledger_version"

    alltypes = DataType[]
    for (k, mods) in TYPE_MODS
        for m in mods
            push!(alltypes, getfield(m, k))
        end
    end
    typemap = Dict([replace(string(t), "RomeoDFT" => "Occupations") => t for t in alltypes])
    for t in keys(TYPE_MODS)
        typemap["Occupations." * string(t)] = getfield(RomeoDFT, t)
    end
    typemap["Occupations.TrialOrigin"] = RomeoDFT.TrialOrigin
    typemap["Occupations.MixingMode"] = RomeoDFT.MixingMode
    for sys in filter(x -> isdefined(RomeoDFT, x) && getfield(RomeoDFT, x) isa DataType && getfield(RomeoDFT, x) <: System, names(RomeoDFT, all=true))
        typemap["Occupations."*string(sys)] = getfield(RomeoDFT, sys)
    end
   
    saved_components = JLD2.jldopen(joinpath(rootdir, "ledger.jld2"), "r", typemap=typemap) do f
        for name in fieldnames_to_save(l)
            curfield = getfield(l, name)
            setfield!(l, name, get(f, "$name", curfield))
        end
        
        if haskey(f, "components")
            compdict = Dict{DataType, Overseer.AbstractComponent}()
            for cname in keys(f["components"])
                c = f["components"][cname]
                compdict[eltype(c)] = c
            end
            l.ledger.entities = get(f, "entities", l.ledger.entities)
            l.ledger.free_entities = get(f, "free_entities", l.ledger.free_entities)
            return compdict
        else
            tl = f["ledger"]
            l.ledger.entities = tl.entities
            l.ledger.free_entities = tl.free_entities
            return components(tl)
        end
    end

    Overseer.ledger(l).components = saved_components
    
    if l.mode == :searching
        l.mode = :search
    end

    current_components = Overseer.ledger(l).components
    for (t, mods) in TYPE_MODS
        # do sequential update
        mod_types = getfield.(mods, t)
        curid = findfirst(T -> haskey(current_components, T), mod_types)
        if curid === nothing || curid == length(mod_types)
            continue
        end
        
        old_T = mod_types[curid]
        new_T = mod_types[end]
        Overseer.ensure_component!(l, new_T)
                
        newcomp = components(l)[new_T]
        for e in @entities_in(current_components[old_T])
            newcomp[e] = version_convert(e[old_T], mod_types, curid)
        end
                
        if newcomp isa PooledComponent
            Overseer.make_unique!(newcomp)
        end

        delete!(current_components, old_T)
    end
    
    set_searcher_stages!(l, l.mode)
    Overseer.ensure_component!(l, Log)

    for e in valid_entities(l)
        if e ∉ l.ledger[Log]
            l.ledger[Log][e] = Log()
        end
    end
    
    return l
end

load_archived(l::Searcher) = load_archived!(load(l))

function load_archived!(l::Searcher)
    for e in @safe_entities_in(l[Archived])
        if ispath(e.archive_path)
            try
                all = JLD2.load(e.archive_path)["components"]
                for a in all
                    l[e] = a
                end
                pop!(l[Archived], e)
            catch
                @warn "Something went wrong for $(e.e).\n$(stacktrace(catch_backtrace()))"
            end
        end
    end
    return l
end

### SEARCHING
function set_mode!(l::Searcher, mode::Symbol)
    set_searcher_stages!(l.ledger, mode)
    l.mode = mode
    return l
end

function setup_scf(scf_file, supercell;
                      Hubbard_maxstep = 100,
                      Hubbard_mixing_beta = 0.4,
                      Hubbard_strength = 1.0,
                      Hubbard_conv_thr = 0.1,
                      electron_maxstep = 500, kwargs...)
    @assert ispath(scf_file) ArgumentError("scf file not found")

    template = DFC.FileIO.qe_parse_calculation(scf_file)

    calc = Calculation{QE}(; name = "scf", exec = Exec(; name = "pw", path = "pw.x"),
                           flags = template.flags, data = template.data)
    calc[:system][:Hubbard_maxstep]     = Hubbard_maxstep
    calc[:system][:Hubbard_mixing_beta] = Hubbard_mixing_beta
    calc[:system][:Hubbard_strength]    = Hubbard_strength
    calc[:system][:Hubbard_conv_thr]    = Hubbard_conv_thr
    calc[:calculation] = "scf"
    calc[:verbosity] = "high"
    calc[:electron_maxstep] = electron_maxstep
    for (k, v) in kwargs
        calc[k] = v
    end
    
    if any(!isequal(1), supercell)
        kpts = DFC.data(calc, :k_points).data
        set_kpoints!(calc, (ceil.(Int, kpts[1:3] ./ supercell)..., kpts[4:6]...))
    end
    return calc
end

function setup_structure(structure_file, supercell, primitive)
    @assert ispath(structure_file) ArgumentError("Structure file not found")
    
    str = splitext(structure_file)[end] == ".in" ?
          DFC.FileIO.qe_parse_calculation(structure_file)[end] :
          Structure(structure_file)
    str = primitive ? Structures.find_primitive(str) : str

    if length(supercell) != 3
        error("wrong supercell specified")
    end
    if any(!isequal(1), supercell)
        str = Structures.create_supercell(str, (supercell .- 1)...)
    end

    
    loc_server = local_server()
    pseudosets = load(loc_server, PseudoSet(""))
    while isempty(pseudosets)
        @info "No PseudoSets found on Server $(loc_server.name), please configure one now..."
        RemoteHPC.configure()
        pseudosets = load(loc_server, PseudoSet(""))
    end
    pseudo_choice = request("Select pseudoset:", RadioMenu(pseudosets))
    if pseudo_choice < 0
        return
    end

    set_pseudos!(str, load(loc_server, PseudoSet(pseudosets[pseudo_choice])))

    atsyms = unique(map(x -> string(x.name), str.atoms))
    choices = request("Select magnetic elements:", MultiSelectMenu(atsyms))

    mag = (1e-5, -1e-5)
    magcount = 1
    for c in choices
        U = RemoteHPC.ask_input(Float64, "Set U for element $(atsyms[c])")
        for a in filter(x->x.name == Symbol(atsyms[c]), str.atoms)
            a.dftu.U = U
            a.magnetization = [0.0, 0.0, mag[mod1(magcount, 2)]]
            magcount += 1
        end
    end

    magatsyms = map(x-> Symbol(atsyms[x]), collect(choices))
    sort!(str.atoms, by = x -> x.name in magatsyms ? 0 : typemax(Float64)) 

    @info "Structure that will be used:"
    display(str)
    return str
end

function setup_ServerInfo()
    server_names = load(local_server(), Server())
    server_choice = request("Select Server to run on:", RadioMenu(server_names))
     
    s = Server(server_names[server_choice])
    execs = load(s, Exec(; path = "pw.x"))
    while isempty(execs)
        @info "No pw executables found on Server $(s.name), please configure one now."
        RemoteHPC.configure()
        execs = load(s, Exec(; path = "pw.x"))
    end
    
    println("Server $(s.name):")
    
    pw_exec = execs[request("Please select pw exec", RadioMenu(execs))]
    envs = load(s, Environment())
    while isempty(envs)
        @info "No Environment found on Server $(s.name), please configure one now."
        RemoteHPC.configure()
        envs = load(s, Environment())
    end
    multi_env = envs[request("Please multi node Environment", RadioMenu(envs))]
    priority = RemoteHPC.ask_input(Int, "Please set priority", 5)

    return ServerInfo(s.name, pw_exec, multi_env, priority)
end

function setup_search(name, scf_file, structure_file=scf_file;
                      nflies = 10,
                      mixing = EulerAngleMixing,
                      γ = 1.0,
                      α = 0.5,
                      β = 0.5,
                      sleep_time = 30,
                      primitive = false,
                      supercell = [1,1,1],
                      unique_thr=1e-2,
                      mindist = 0.25,
                      stopping_unique_ratio  = 0.1,
                      stopping_n_generations = 3,
                      relax_unique=false,
                      relax_base=false,
                      force_convergence_threshold=1e-3,
                      energy_convergence_threshold=1e-4,
                      ion_dynamics="bfgs",
                      cell_dynamics="bfgs",
                      symmetry=true,
                      variable_cell=true,
                      
                      hp_base = false,
                      hp_unique = false,
                      hp_nq = (2,2,2),
                      hp_conv_thr_chi = 1e-4,
                      hp_find_atpert = 2,
                      hp_U_conv_thr = 0.1,
                      
                      kwargs...)
                      
    dir = searchers_dir(name)
    
    if ispath(joinpath(dir, string(DATABASE_VERSION), "ledger.jld2"))
        choice = request("Overwrite previous results?", RadioMenu(["no", "yes"]))
        if choice < 0
            return
        elseif choice == 1
            l = load(Searcher(dir))
            l.sleep_time = sleep_time
            return l
        else
            rm(dir, recursive=true)
        end
    else
        mkpath(dir)
    end

    calc = setup_scf(scf_file, supercell; kwargs...)
    str = setup_structure(structure_file, supercell, primitive)

    magats = filter(ismagnetic, str.atoms)
    suppress() do
        calc[:system][:Hubbard_conv_thr] *= length(magats)
    end

    l = Searcher(; rootdir = dir, sleep_time = sleep_time)

    sim_e = Entity(l, setup_ServerInfo(),RandomSearcher(nflies),
                           Template(deepcopy(str), deepcopy(calc)),
                           IntersectionSearcher(mindist, 100),
                           StopCondition(stopping_unique_ratio, stopping_n_generations),
                           Generation(1))
    # BaseCase simulation entity
    base_e = Entity(l, BaseCase(),
                       Template(deepcopy(str),
                       deepcopy(calc)),
                       Generation(1))
                       
    unique_e = Entity(l, Unique(unique_thr, true), Generation(1))
    
    relset = RelaxSettings(force_convergence_threshold, energy_convergence_threshold, ion_dynamics, cell_dynamics, symmetry, variable_cell)
    if relax_unique
        l[unique_e] = relset
    end
    if relax_base
        l[base_e] = relset
    end

    hpset = HPSettings(hp_nq, hp_conv_thr_chi, hp_find_atpert, hp_U_conv_thr)
    if hp_unique
        l[unique_e] = hpset
    end
    if hp_base
        l[base_e] = hpset
    end
    
    ishybrid = haskey(calc, :input_dft) || haskey(calc, :exxdiv_treatment)
    
    if ishybrid
        l[base_e] = Hybrid()
        l[sim_e] = Hybrid()
    end

    set_mode!(l, :search)
    save(l)
    return l
end

"Writes out a final report to `<rootdir>/report.out` after the searcher has stopped."
function final_report(l::Searcher)
    structure_name = Structures.name(l[Template][1].structure)

    n_trials    = 0
    n_converged = 0
    n_unique    = 0
    for e in @entities_in(l, Results && (Intersection || RandomSearcher))
        n_trials += 1

        if e.converged
            n_converged += 1
        end
        
        if e in l[Unique]
            n_unique += 1
        end
    end

    unique_thr = round(l[Unique][1].thr, digits=2)
    
    open(joinpath(l.rootdir, "report.out"), "w") do f
        write(f, "Global search for system $structure_name finished.\n")
        write(f, "Generations:         $(maximum_generation(l))\n")
        write(f, "Trials:              $n_trials\n")
        write(f, "Converged:           $n_converged\n")
        write(f, "Unique (thr = $unique_thr): $n_unique\n")
        write(f, "Unique/Trials:       $(n_unique/n_trials)\n\n")
        write_groundstate(f, l)
        plot_states(f, l)
        plot_evolution(f, l; color=false)
    end
end

function write_groundstate(io::IO, l::Searcher)
    if any(x->x.converged, l[Results])
        verify_groundstates!(l)
        
        groundstate = ground_state(l)
        
        println(io, "Groundstate: ")
        println(io, "\t$(Entity(groundstate)), Generation($(groundstate.generation))")
        println(io, "\tEnergy:  $(groundstate.total_energy) Ry")
        
        s = groundstate[Results].state
        magmoms = join(string.(round.(s.magmoms, digits=3)), " ")
        println(io, "\tMagmoms: ", magmoms)
    else
        println(io, "No Results yet")
    end
    println(io)
end

function plot_evolution(io::IO, l::Searcher; color=true)
    stop, n_unique, n_total = stop_check(l)
    if length(n_unique) > 1
        p1 = Main.UnicodePlots.lineplot(n_unique./n_total, title = "Unique/Trial", xlabel = "Generation", color=:white)
        minima = min_energy_per_generation(l) 
        p2 = minima[1] == 0 ? Main.UnicodePlots.lineplot(2:length(minima), view(minima,2:length(minima)), title = "Minimum energy [Ry]", xlabel = "Generation", color=:white) :
                              Main.UnicodePlots.lineplot(1:length(minima), view(minima,1:length(minima)), title = "Minimum energy [Ry]", xlabel = "Generation", color=:white)
        if StopCondition in l && !isempty(l[StopCondition])
            Main.UnicodePlots.hline!(p1, l[StopCondition][1].unique_ratio, color=:cyan, name="stop ratio")
        end
        if color
            println(io, string(p1; color=true))
        else
            println(io, p1)
        end
        println(io)
        Main.UnicodePlots.hline!(p2, minima[end], color=:cyan, name="minimum")
        if color
            println(io, string(p2; color=true))
        else
            println(io, p2)
        end
        println(io)
    end
end

function plot_states(io::IO, l::Searcher)
    es = filter(x -> x.converged, @entities_in(l, Results && Template))
    if !isempty(es)
        energies = relative_energies(es)
        magmoms  = map(x->sum(x.state.magmoms), es)
        p1 = Main.UnicodePlots.scatterplot(magmoms, energies, title = "States", xlabel = "Total m", ylabel = "E rel [eV]", color=:white, marker="+")
        println(io, p1)
        println(io)
    end
end

function status(io::IO, l::Searcher)
    header = "$(searcher_name(l)) @ Generation($(maximum(x->x.generation, l[Generation], init=0)))"
    horstring = "+" * "-"^(length(header)+2) * "+"
    println(io, horstring)
    println(io, "| ", header, " |")
    println(io, horstring)
    println(io)
    if l.loop_error
        status = "ERRORED"
    elseif l.loop === nothing || istaskdone(l.loop)
        status = l.finished ? "FINISHED" : "STOPPED"
    else
        if l.mode == :postprocess
            status = "POSTPROCESSING"
        elseif l.mode == :search
            status = "SEARCHING"
        elseif l.mode == :cleanup
            status = "CLEANUP"
        else
            status = "UNKNOWN"
        end
    end
    println(io, "Status:        $status")
    println(io, "Unique states: $(length(l[Unique]))")
    println(io, "Total Trials:  $(length(@entities_in(l, Results && !Parent)))")
    
    println(io)
    write_groundstate(io, l)
    plot_states(io, l)
    plot_evolution(io, l)
    es = sort(collect(@entities_in(l, SimJob)), by = x -> x.e.id)
    if !isempty(es)
        println(io, "Current simulation jobs:")
        curgen = isempty(l[Generation]) ? 0 : maximum(x->x.generation, l[Generation].data)
        es_str = ["ID";map(e -> string(e.e.id), es)]
        states = Vector{String}(undef, length(es))
        @sync for (ie, e) in enumerate(es)
            Threads.@spawn begin
                states[ie] = string(state(l[SimJob][e].job))
            end
        end
            
        status_str = ["Status"; map(enumerate(es)) do (ie, e)
            if e in l[Done]
                return "Cleanup"
            elseif e in l[Error]
                return "Errored"
            elseif e in l[SimJob]
                if e in l[Submit]
                    return "Awaiting Submission"
                else
                    return states[ie]
                end
            else
                return "Awaiting Creation" 
            end
        end]
        servers_str = ["Server"; map(es) do e
            if e in l[SimJob]
                return l[SimJob][e].job.server
            else
                return "Unknown"
            end
        end]
        longest_e = maximum(length, es_str)
        longest_status = maximum(length, status_str)
        longest_server = maximum(length, servers_str)
        horizontal_line = "+" * "-"^(longest_e+2) * "+" * "-"^(longest_status+2) * "+" * "-"^(longest_server + 2) * "+"

        pad(str, longest) = str * " "^(longest - length(str))
        println(io, horizontal_line)
        println(io, "| ", pad(es_str[1], longest_e), " | ", pad(status_str[1], longest_status), " | ", pad(servers_str[1], longest_server), " |")
        println(io, horizontal_line)
        for i = 2:length(es_str)
            println(io, "| ", pad(es_str[i], longest_e), " | ", pad(status_str[i], longest_status), " | ", pad(servers_str[i], longest_server), " |")
        end
        println(io, horizontal_line)
        println(io)
    end
    err_es = @entities_in(l, SimJob && Error)
    if !isempty(err_es)
        println(io, "Errored simulation jobs:")
        for e in err_es
            println(io, "$(e.e):\n", e[Error])
        end
        println(io)
    end
    c = 0
    for e in @entities_in(l, SimJob && (Parent || NSCFSettings || BandsSettings || ProjwfcSettings))
        c += 1
    end
    if c != 0
        println(io, "Running post processing jobs: $c")
        println(io)
    end
    c = 0
    for e in @entities_in(l, Done && (Parent || NSCFSettings || BandsSettings || ProjwfcSettings))
        c += 1
    end
    if c != 0
        println(io, "Finished post processing jobs: $c")
        println(io)
    end
    println(io, "Runner task:")
    println(io, l.loop)
   
end

function add_search_entity!(m::AbstractLedger, search_e::Overseer.AbstractEntity, components...)
    e = Entity(m)
    for c in components
        comp = m[typeof(c)]
        if comp isa PooledComponent 
            comp[e] = search_e
        else
            comp[e] = c
        end
    end
    m[Template][e] = search_e
    return e
end

"""
    set_status!(m::Searcher, e::AbstractEntity, s)

Sets the status of an [`Entity`](@ref). This is used by the various systems that are part of
the workflow to make an [`Entity`](@ref) flow through it in the correct way.
"""
function set_status!(m, e, s::T) where {T}
    if e ∉ m[T]
        for t in (Submit, Submitted, Running, Completed, Pulled)
            trypop!(m[t], e)
        end
        m[e] = s
    end
end

function set_status!(m, e, s::RemoteHPC.JobState)
    c = nothing
    if s == RemoteHPC.Submitted 
        c = Submitted()
    elseif isparseable(s)
        c = Completed()
    elseif s == RemoteHPC.Running
        c = Running()
    elseif s == RemoteHPC.Pending
        c = Running()
    end
    if c !== nothing
        set_status!(m, e, c)
    end
end

"""
    should_rerun(m::Searcher, e::AbstractEntity), datatypes...)

Signals that an entity should rerun. This is mainly used in self-consistent kind of calculations such as
the `vc-relax` - `hp` loop. It will be picked up on by the [`Rerunner`](@ref) system that will
remove the [`Entity`](@ref) from the `datatypes` components.  
"""
function should_rerun(m, e, datatypes...)
    rer = m[ShouldRerun]
    if e in rer
        for d in datatypes
            push!(rer[e], d)
        end
    else
        m[e] = ShouldRerun(datatypes...)
    end
end


function stop(l::Searcher)
    if l.loop !== nothing && !istaskdone(l.loop)
        l.stop = true
        while l.stop && !istaskdone(l.loop)
            sleep(0.1)
        end
    end
    l.stop=false
end
 
function RemoteHPC.start(l; kwargs...)
    stop(l) # Should be noop 
    Main.Revise.revise()
    Base.invokelatest() do
        l.loop = Threads.@spawn loop(l; kwargs...)
    end
end

"""
    loop(l::Searcher; verbosity=l.verbosity, sleep_time=l.sleep_time)

The main loop that executes the global searching.
"""
function loop(l::Searcher; verbosity=l.verbosity, sleep_time=l.sleep_time)
    if verbosity >= 0
        l.verbosity = verbosity
    end
    l.sleep_time = sleep_time
    http_logpath = joinpath(l.rootdir, "HTTP.log")
    logpath = joinpath(l.rootdir, "log.log")

    prevledger = joinpath(storage_dir(l), "ledger.jld2")
    if ispath(prevledger)
        cp(prevledger, joinpath(storage_dir(l),"$(now())_bak.jld2"))
    end
    
    logger = RemoteHPC.TimestampLogger(RemoteHPC.TeeLogger(RemoteHPC.NotHTTPLogger(RemoteHPC.TimeBufferedFileLogger(logpath, interval=0.0001)),
                                                           RemoteHPC.HTTPLogger(RemoteHPC.TimeBufferedFileLogger(http_logpath, interval=0.0001))))
    simn = simname(l)
    with_logger(logger) do
        LoggingExtras.withlevel(LoggingExtras.Debug; verbosity=l.verbosity) do
                                   
            while !isalive(local_server())
                @error "Local server not alive"
                sleep(1)
            end
            
            save(l)
            if length(l[Results]) == 0
                @debug "Starting search for global minimum state."
            else
                @debug "Restarting search for global minimum state."
            end 
        end
    end
    if isempty(l[StopCondition])
        Entity(l, StopCondition())
    end
    while !l.stop
        curt = now()
        if ispath(logpath) && filesize(logpath) > 1e8
            write(logpath, "")
        end
        if ispath(http_logpath) && filesize(http_logpath) > 1e8
            write(http_logpath, "")
        end
        with_logger(logger) do
            if all(isalive, all_servers(l))
                try 
                    RemoteHPC.@timeout 600 update(l)
                    l.loop_error = false
                catch e
                    RemoteHPC.log_error(e)
                    l.loop_error = true
                end
            end
            secs = round(Int, l.sleep_time)
            ms = round(Int, (l.sleep_time - secs)*1000)
            while !l.stop &&  round(now() - curt, Second) < Second(secs) + Millisecond(ms)
                sleep(0.5)
            end
            save(l)
        end
    end
    @debug "Stopping loop"
    @info "Stopping all pending and Submitted jobs."
    lck = ReentrantLock()
    @sync for e in @entities_in(l, SimJob)
        Threads.@spawn begin
            if state(e.job) ∈ (RemoteHPC.Pending, RemoteHPC.Submitted)
                abort(e.job)
                lock(lck)
                set_status!(l, e, Submit())
                unlock(lck)
            end
        end
    end
    l.loop = nothing
    return l.stop = false
end

"""
    isolate_entity(l::Searcher, e::AbstractEntity, from_scratch=true)

Creates a new [`Searcher`](@ref) containing only the chosen [`Entity`](@ref). 
"""
function isolate_entity(l::Searcher, e::AbstractEntity, from_scratch=true)
    rootdir = l.rootdir * "_isolate_$(Entity(e).id)"
    if ispath(rootdir)
        rm(rootdir, recursive=true)
    end
    
    out_l = Searcher(rootdir)
    
    set_mode!(out_l, mode(l))
    out_l.verbosity = l.verbosity
    out_l.sleep_time = l.sleep_time
    
    Entity(out_l, l[Entity(1)]...)
    if !from_scratch
        Entity(out_l, l[e]...)
    else
        new_e = Entity(out_l)
        for c in (Template, Trial, RelaxSettings, HPSettings, ProjwfcSettings, NSCFSettings, BandsSettings, Parent, Generation)
            if c in l[e]
                out_l[new_e] = l[e][c]
            end
        end
    end

    
    save(out_l)
    return out_l
end

"""
    create_child!(l::AbstractLedger, parent_entity, components...)

Creates a child from the `parent_entity`, using its [`Generation`](@ref), the state in its [`Results`](@ref) for [`Trail`](@ref), its [`Template`](@ref),
and any other `Component` in `components`.
"""
function create_child!(l::AbstractLedger, parent_entity, components...)
    gen = parent_entity in l[Generation] ? l[Generation][parent_entity] : Generation(maximum_generation(l))

    ppe = Entity(l, components...)
    
    for c in l[parent_entity]
        cT = typeof(c)

        comp = l[cT]
        if ppe in comp || !(cT <: PostProcessSettings)
            continue
        end

        if comp isa PooledComponent
            comp[ppe] = parent_entity
        else
            comp[ppe] = deepcopy(c)
        end
    end
    
    l[ppe] = Parents(Entity(parent_entity))

    if ppe ∉ l[Trial]
        l[ppe] = parent_entity in l[Results] && l[Results][parent_entity].converged ?
                 Trial(l[Results][parent_entity].state, PostProcess) :
                 Trial(l[Trial][parent_entity].state, PostProcess)
    end

    if ppe ∉ l[Generation]
        l[ppe] = gen
    end

    if ppe ∉ l[Template]
        l[ppe] = deepcopy(l[Template][parent_entity])
    end

    if parent_entity ∈ l[Children]
        push!(l[Children][parent_entity], ppe)
    else
        l[Children][parent_entity] = Children(ppe)
    end
        
    return ppe
end

"""
    create_postprocess_child!(l::AbstractLedger, parent_entity, components...)

First calls [`create_child`](@ref), then will look at the components held by the [`Unique`](@ref) entity to decide what other postprocessing components should be added."""
function create_postprocess_child!(l, parent_entity, comps...)
    
    ppe = create_child!(l, parent_entity, comps...)
    unique_e = l[entity(l[Unique], 1)]
    
    for c in components(unique_e)
        if c isa PooledComponent && eltype(c) != Unique && ppe ∉ c
            c[ppe] = unique_e
        end
    end
    return ppe
end

function recurse_parents(f::Function, l::AbstractLedger, e::AbstractEntity, inclusive=false)
    parent_comp = l[Parents]
    if inclusive
        f(e)
    end

    if e ∉ parent_comp
        return
    end
    
    for c in parent_comp[e].parents
        f(c)
        recurse_parents(f, l, c)
    end
end

function recurse_children(f::Function, l::AbstractLedger, e::AbstractEntity, inclusive=false)
    child_comp = l[Children]
    if inclusive
        f(e)
    end
    
    if e ∉ child_comp
        return
    end
    
    for c in child_comp[e].children
        f(c)
        recurse_children(f, l, c)
    end
end

function youngest_child(l, e)
    out = Entity(e)
    recurse_children(l, e) do child
        out = child.id > out.id ? child : out
    end
    return out
end

function oldest_parent(l, e)
    out = Entity(e) 
    recurse_parents(l, e) do parent
        out = parent.id < out.id ? parent : out
    end
    return out
end


function all_children_done(l, parent)
    done_comp = l[Done]
    out = true
    recurse_children(l, parent, true) do child
        if child ∉ done_comp
            out &= false
        end
    end
    return out
end
   
function gather_logs(l, e)
    
    all_logs = Pair{Int, Log}[]
    
    log_comp = l[Log]
    
    recurse_parents(l, e) do parent
        push!(all_logs, parent.id => log_comp[parent])
    end
    
    recurse_children(l, e, true) do child
        push!(all_logs, child.id => log_comp[child])
    end

    sort!(all_logs, by = x -> first(x))

    full_log = Log()
    for (id, log) in all_logs
        append!(full_log, log)
    end
    
    return full_log
end
