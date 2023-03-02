using JLD2: jldopen
using REPL.TerminalMenus

dft_energy(res) = res.total_energy - res.Hubbard_energy

function core_stage()
    return Stage(:core, [JobCreator(), [JobSubmitter(),Cleaner()], [JobMonitor(), OutputPuller()], ResultsProcessor(), RelaxProcessor(), ErrorCorrector(), UniqueExplorer(), Relaxor(), Stopper()])
end

function cleanup_stage()
    return Stage(:cleanup, [JobCreator(), JobMonitor(), OutputPuller(), ResultsProcessor(), RelaxProcessor(), UniqueExplorer(), Cleaner(), Stopper()])
end

function intersection_stage()
    return Stage(:intersection, [Intersector(), RandomTrialGenerator()])
end
function firefly_stage()
    return Stage(:firefly, [FireFly(), PostFireflyExplorer(), Archiver()])

end

function search_stages()
    return [Stage(:main, [intersection_stage(), core_stage()])]
end

Base.@kwdef mutable struct Searcher <: AbstractLedger
    rootdir::String
    loop::Union{Nothing,Task} = nothing
    stop::Bool = false
    finished::Bool = false
    sleep_time::Float64 = 30
    loop_error::Bool = false
    mode::Symbol = :postprocess # What is the searcher doing i.e. searching or merely postprocessing
    ledger::Ledger = Ledger((mode == :postprocess ? [core_stage()] : search_stages())...)
end

function Searcher(dir::AbstractString; kwargs...)
    dir = isdir(dir) ? dir : searchers_dir(dir)
    Searcher(; rootdir = abspath(dir), kwargs...)
end
Overseer.ledger(l::Searcher) = l.ledger

function set_searcher_stages!(l::AbstractLedger, s::Symbol)
    if s == :postprocess
        l.stages = [core_stage()]
    elseif s == :search
        l.stages = search_stages()
    elseif s == :cleanup
        l.stages = [cleanup_stage()]
    else
        error("Searcher stage $s not recognized...")
    end
    prepare(l)
end

set_searcher_stages!(l::Searcher) = set_searcher_stages!(l.ledger, l.mode)

function set_mode!(l::Searcher, mode::Symbol)
    set_searcher_stages!(l.ledger, mode)
    l.mode = mode
    return l
end

storage_dir(l::Searcher) = joinpath(l.rootdir, string(DATABASE_VERSION))
function RemoteHPC.load(l::Searcher;
                                 version = nothing)
    if version === nothing
        versions = []
        for d in readdir(l.rootdir)
            if isdir(joinpath(l.rootdir, d)) 
                try
                    push!(versions, VersionNumber(d))
                catch
                    nothing
                end
            end
        end
        version = maximum(versions)
    end
    return load(joinpath(l.rootdir, string(version)), l)
end

function searcher_name(s::AbstractString)
    splt = strip.(split(s, searchers_dir()), '/')
    if length(splt) > 1
        return splt[end]
    else
        return s
    end
end
searcher_name(m::AbstractLedger) = searcher_name(m.rootdir)

function local_dir(m::AbstractLedger, e) 
    return abspath(joinpath(m.rootdir, "job_backups","$(Entity(e).id)"))
end

function remote_dir(m::AbstractLedger, e) 
    fdir = searcher_name(m)
    return joinpath("RomeoDFT", fdir, "$(Entity(e).id)")
end

function simname(m)
    fdir = searcher_name(m)
    return replace(fdir, "/" => "_")
end

function RemoteHPC.save(l::Searcher)
    return (mkpath(storage_dir(l)); save(storage_dir(l), l))
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

function energies(l)
    out = [zeros(i.n_fireflies, i.current_generation) for i in l[Simulation].data]
    for e in @entities_in(l, Generation && Simulation && Results)
        pid = Overseer.pool(l[Simulation], e)
        gen = e.generation
        o = out[pid]
        o[findfirst(i -> iszero(o[i, gen]), 1:size(o, 1)), gen] = dft_energy(e)
    end
    return out
end

function states(l)
    out = []
    for r in @entities_in(l[Results])
        if !r.converged
            continue
        end
        if !any(x -> sum(abs.(x.state.magmoms .- r.state.magmoms)) < 1e-2 ||
                    Euclidean()(x.state, r.state) < 1e-2, out)
            push!(out, r)
        end
    end
    return out
end

function setup_scf(scf_file, supercell)
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

function setup_structure(structure_file, primitive, supercell)
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

    
    local_server = Server(gethostname())
    pseudosets = load(local_server, PseudoSet(""))
    while isempty(pseudosets)
        @info "No PseudoSets found on Server $(local_server.name), please configure one now..."
        RemoteHPC.configure()
        pseudosets = load(local_server, PseudoSet(""))
    end
    pseudo_choice = request("Select pseudoset:", RadioMenu(pseudosets))
    if pseudo_choice < 0
        return
    end

    set_pseudos!(str, load(local_server, PseudoSet(pseudosets[pseudo_choice])))

    atsyms = unique(map(x -> string(x.name), str.atoms))
    choices = request("Select magnetic elements:", MultiSelectMenu(atsyms))

    for c in choices
        U_str = ask_input("Set U for element $(atsyms[c]): ")
        U = parse(Float64, U_str)
        for a in filter(x->x.name == Symbol(atsyms[c]), str.atoms)
            a.dftu.U = U
            a.magnetization = [0.0, 0.0, 1e-5]
        end
    end

    @info "Structure that will be used:"
    display(str)
    return str
end

function setup_ServerInfo()
    server_names = load(local_server, Server())
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
    priority = parse(Int, ask_input("Please set priority [Int] (use 5 for default): "))

    return ServerInfo(s.name, pw_exec, multi_env, priority)
end

function setup_search(name, scf_file, structure_file=scf_file;
                      nflies = 10,
                      mixing = EulerAngleMixing,
                      γ = 1.0,
                      α = 0.5,
                      β = 0.5,
                      sleep_time = 30,
                      Hubbard_maxstep = 100,
                      Hubbard_mixing_beta = 0.4,
                      Hubbard_strength = 1.0,
                      Hubbard_conv_thr = 0.1,
                      primitive = false,
                      supercell = [1,1,1],
                      electron_maxstep = 500,
                      verbosity = 0,
                      unique_thr=1e-2,
                      mindist = 0.25,
                      stopping_unique_ratio = 0.1,
                      stopping_n_generations =3,
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
        end
    else
        mkpath(dir)
    end

    calc = setup_scf(scf_file, supercell)
    str = setup_structure(structure_file, supercell, primitive) 

    l = Searcher(; rootdir = dir, sleep_time = sleep_time)

    Entity(l, setup_ServerInfo())
   

    ishybrid = haskey(calc, :input_dft) || haskey(calc, :exxdiv_treatment)
    
    # BaseCase simulation entity
    base_e = Entity(l, BaseCase(),
                       Template(deepcopy(str),
                       deepcopy(calc)),
                       Generation(1))
                       
    if ishybrid
        l[base_e] = Hybrid()
    end

    sim_entity = Entity(l, RandomSearcher(nflies),
                           Template(deepcopy(str), deepcopy(calc)),
                           Unique(unique_thr, true),
                           IntersectionSearcher(mindist, 100),
                           StopCondition(stopping_unique_ratio, stopping_n_generations),
                           Generation(1))
                           
    if ishybrid
        l[sim_entity] = Hybrid()
    end

    set_mode!(l, :search)
    save(l)
    return orchestrator_eval("start_searcher(\"$(l.rootdir)\"; verbosity=$(verbosity), sleep_time=Int($(l.sleep_time)))")
end

function ground_state(l; by = x -> x.total_energy)
    min_entity = Entity(0)
    min_energy = typemax(Float64)
    for e in @entities_in(l, Results && FlatBands)
        if e.converged
            test = by(e)
            if test < min_energy
                min_energy = test
                min_entity = e.e
            end
        end
    end
    return l[min_entity]
end

function unique_states(es; thr = 1e-4, bands = true)
    iu = ones(length(es))
    @inbounds for i in 1:length(es)-1
        e1 = es[i]
        if iu[i] == 1
            Threads.@threads for j in i+1:length(es)
                if iu[j] == 1
                    e2 = es[j]
                    dist = bands ? sssp_distance(e1.bands, e2.bands, e1.fermi) :
                           Euclidean()(e1.state, e2.state)
                    if dist < thr
                        iu[j] = 0
                    end
                end
            end
        end
    end
    return es[findall(isequal(1), iu)]
end

function unique_states(l::AbstractLedger; kwargs...)
    return unique_states(filter(x -> x.converged,
                                  @entities_in(l[Results] && l[FlatBands])); kwargs...)
end

function final_report(l::Searcher)
    open(joinpath(l.rootdir, "report.out"), "w") do f
        write(f, "Global search for system $(Structures.name(l[Template][1].structure)) finished.\n")
        write(f, "Generations:         $(maximum(x->x.generation, l[Generation], init=0))\n")
        write(f, "Trials:              $(length(l[Results]))\n")
        write(f, "Converged:           $(length(filter(x->x.converged, l[Results])))\n")
        write(f, "Unique (thr = $(round(l[Unique][1].thr, digits=2))): $(length(l[Unique])-1)\n")
        write(f, "Unique/Trials:       $(length(l[Unique])/length(l[Results]))\n\n")
        write_groundstate(f, l)
        plot_states(f, l)
        plot_evolution(f, l; color=false)
    end
end

function write_groundstate(io::IO, l::Searcher)
    if any(x->x.converged, l[Results])
        groundstate = ground_state(l)
        println(io, "Groundstate: ")
        if groundstate in l[Generation]
            println(io, "\t$(groundstate.e), Generation($(groundstate.generation))")
        end
        println(io, "\tEnergy:  $(groundstate.total_energy) Ry")
        s = groundstate[Results].state
        println(io, "\tMagmoms: ", join(string.(round.(s.magmoms, digits=3)), " "))
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
    es = filter(x -> x.converged, @entities_in(l[Results] && l[Template]))
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
    println(io, "Total Trials:  $(length(@entities_in(l, Results)))")
    
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
    for e in @entities_in(l, SimJob && NSCFSettings)
        c += 1
    end
    if c != 0
        println(io, "Running post processing jobs: $c")
        println(io)
    end
    c = 0
    for e in @entities_in(l, Unique && Results && Done)
        c += 1
    end
    if c != 0
        println(io, "Finished post processing jobs: $c")
        println(io)
    end
    println(io, "Runner task:")
    println(io, l.loop)
   
end

function stop(l::Searcher)
    if l.loop !== nothing && !istaskdone(l.loop)
        l.stop = true
        while l.stop && !istaskdone(l.loop)
            sleep(0.1)
        end
    end
end
 
function RemoteHPC.start(l; kwargs...)
    stop(l) # Should be noop 
    Main.Revise.revise()
    Base.invokelatest() do
        l.loop = Threads.@spawn loop(l; kwargs...)
    end
end

function loop(l::Searcher; verbosity=0, sleep_time=l.sleep_time)
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
        LoggingExtras.withlevel(LoggingExtras.Debug; verbosity=verbosity) do
                                   
            while !isalive(local_server())
                @error "Local server not alive"
                sleep(1)
            end
            
            @debugv 2 "[START] Saving" 
            save(l)
            @debugv 2 "[STOP] Saving"
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
            LoggingExtras.withlevel(LoggingExtras.Debug; verbosity=verbosity) do
                if isalive(local_server())
                    try 
                        @debugv 1 "[START] Updating" 
                        RemoteHPC.@timeout 600 update(l)
                        l.loop_error = false
                        @debugv 1 "[STOP] Updating" 
                    catch e
                        RemoteHPC.log_error(e)
                        l.loop_error = true
                    end
                end
                while !l.stop &&  round(now() - curt, Second) < Second(l.sleep_time)
                    sleep(0.5)
                end
                @debugv 2 "[START] Saving" 
                save(l)
                @debugv 2 "[STOP] Saving" 
            end
        end
    end
    @debug "Stopping loop"
    @info "Stopping all pending and Submitted jobs."
    @sync for e in @entities_in(l, SimJob)
        Threads.@spawn begin
            if state(e.job) ∈ (RemoteHPC.Pending, RemoteHPC.Submitted)
                abort(e.job)
                l[e] = Submit()
            end
        end
    end
    l.loop = nothing
    return l.stop = false
end

