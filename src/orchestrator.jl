mutable struct SearchOrchestrator
    searchers::Dict{String, Searcher}
    paused::Bool
    stop::Bool
end
const ORCHESTRATOR = Ref{SearchOrchestrator}()

# TODO: actually print it
loaded_searchers() = collect(keys(ORCHESTRATOR[].searchers))

function list_searchers()
    ls = loaded_searchers()
    to_text = sprint() do io
        for l in ls
            println(io, l)
        end
    end
    return Text(to_text)
end

function load_searcher(p::AbstractString)
    l = load(Searcher(p))
    n = searcher_name(l)
    searchers = ORCHESTRATOR[].searchers
    if haskey(searchers, n)
        throw(ArgumentError("Already loaded searcher at $p"))
    end
    searchers[n] = l
    @info "Loaded searcher and assigned it label $n"
    return n
end

function unload_searcher(n::AbstractString)
    o = ORCHESTRATOR[]
    if !haskey(o.searchers, n)
        throw(ArgumentError("No searcher loaded with label $n"))
    end
    stop(o.searchers[n])
    save(o.searchers[n])
    delete!(o.searchers, n)
    @info "Unloaded searcher $n"
    return n
end

function start_searcher(n::AbstractString; kwargs...)
    o = ORCHESTRATOR[]
    if !haskey(o.searchers, n)
        n = load_searcher(n)
    end
    l = o.searchers[n]
    if l.loop !== nothing && (l.stop || !istaskdone(l.loop))
        throw(ArgumentError("Searcher $n already started"))
    end
    l.finished = false
    start(l; kwargs...)
    @info "Started searcher $n"
    return n
end

function stop_searcher(n::AbstractString)
    o = ORCHESTRATOR[]
    if !haskey(o.searchers, n)
        throw(ArgumentError("No searcher loaded with label $n"))
    end
    l = o.searchers[n]
    stop(l)
    @info "Stopped searcher $n"
    return n
end

function stop_searchers(o::SearchOrchestrator)
    # Gracefully try to stop all loaded searches
    stop_tasks = Pair{String, Task}[]
    for (k, l) in o.searchers
        @debug "Stopping search $k"
        push!(stop_tasks, k => Threads.@spawn stop(l))
    end
    tries = 0
    while tries < 25
        sleep(10)
        not_done = filter(x -> !istaskdone(x[2]), stop_tasks)
        if isempty(not_done)
            break
        end
        @debug "Searches not stopped yet: $(join(map(x->x[1], not_done), " - "))"
        tries += 1
    end
end

function load_searcher(f::Function, n::AbstractString)
    o = ORCHESTRATOR[]
    Main.Revise.revise()
    if haskey(o.searchers, n)
        Base.invokelatest(f, o.searchers[n])
    else
        l = load(Searcher(n))
        Base.invokelatest(f, l)
        save(l)
    end
end

function status()
    o = ORCHESTRATOR[]
    lines = String[]
    
    header = ["Searcher", "Generation", "Status", "Unique", "Trials"]

    data = Matrix{String}(undef, length(o.searchers), 5)
    names = sort(collect(keys(o.searchers)))
    
    for (i, n) in enumerate(names)
        l = o.searchers[n]
        data[i, 1] = n
        data[i, 2] = string(maximum(x->x.generation, l[Generation], init=0))
        if l.loop_error
            data[i, 3] = "ERRORED"
        elseif l.loop === nothing || istaskdone(l.loop)
            data[i, 3] = l.finished ? "FINISHED" : "STOPPED"
        else
            data[i, 3] = uppercase(string(l.mode))
        end
        data[i, 4] = string(length(l[Unique])-1)
        data[i, 5] = string(length(@entities_in(l, Results && !Parent)))
    end

    maxlen = displaysize(stdout)[2]
    horizontal_length = maxlen + 2 + 2 # 2 spaces 2 |
    horizontal_line = "+" * "-"^(maxlen + 2) * "+"

    title = "Global Search Orchestrator"
    mid = ceil(Int, (maxlen + 2 - length(title))/2)

    
    highlighter = Highlighter((d, i, j) -> j == 3 && d[i, j] == "ERRORED", crayon"red")
    
    to_text = sprint() do io
        println(io, horizontal_line)
        println(io, "|" * " "^(maxlen + 2) * "|")
        println(io, "|" * " "^mid * title * " "^(maxlen -  mid - length(title) + 2) * "|")
        println(io, "|" * " "^(maxlen + 2) * "|")
        println(io, horizontal_line)
        pretty_table(io, data, header = header, highlighters = (highlighter,))
    end
    return Text(to_text)
end

function status(n::AbstractString)
    o = ORCHESTRATOR[]
    if !haskey(o.searchers, n)
        throw(ArgumentError("No searcher loaded with label $n"))
    end
    l = o.searchers[n]
    return Text(sprint(io -> status(io, l)))
end

function unload_finished()
    ls = collect(values(ORCHESTRATOR[].searchers))
    finished_ls = filter(x -> x.finished, ls)
    ns = map(x->searcher_name(x), finished_ls)
    if !isempty(ns)
        @debug "Unloading $(length(ns)) finished searchers: $(join(ns, " - "))"
        for (i, n) in enumerate(ns)
            unload_searcher(n)
        end
    end
end

function plot_states(name::AbstractString, outfile::AbstractString; kwargs...)
    load_searcher(name) do l
        p = plot_states(l; title = name, kwargs...)
        savefig(p, outfile)
    end
    return "States of $name plotted to $outfile"
end

function set_mode!(name::AbstractString, mode)
    m = Symbol(mode)
    stop_searcher(name)
    try
        prev_mode = :none
        load_searcher(name) do l
            prev_mode = l.mode
            set_mode!(l, m)
            @info "Switched from $prev_mode to $(l.mode)"
        end
    catch err
        @error exception=err
    finally
        start_searcher(name)
    end
end

## SERVERS
# TODO: should also list stuff that is not usedin a searcher
function list_servers()
    names = String[]
    for l in values(ORCHESTRATOR[].searchers)
        for s in l[ServerInfo].data
            s.server ∉ names && push!(names, s.server)
        end
        local_s = local_server()
        local_s.name ∉ names && push!(names, local_s.name)
    end
    data = Matrix{String}(undef, length(names), 3)
    for (i, n) in enumerate(names)
        data[i, 1] = n
        s = Server(n)
        data[i, 2] = string(isalive(s))
        data[i, 3] = string(RemoteHPC.version(s))
    end
    to_text = sprint() do io
        pretty_table(io, data, header = ["NAME", "ALIVE", "VERSION"])
    end
    return Text(to_text)
end

## ORCHESTRATOR CONTROL
function stop_orchestrator()
    ORCHESTRATOR[].stop = true
end

function pause_orchestrator()
    if !ORCHESTRATOR[].paused
        @info "Pausing Orchestrator, stopping searchers..."
        stop_searchers(ORCHESTRATOR[])
        ORCHESTRATOR[].paused = true
        @info "Orchestrator paused"
    else
        @error "Orchestrator was already paused"
    end
end

function resume_orchestrator()
    if ORCHESTRATOR[].paused
        @info "Resuming Orchestrator, restarting searchers..."
        for n in keys(ORCHESTRATOR[].searchers)
            try
                start_searcher(n)
            catch
                nothing
            end
        end
        @info "Orchestrator resumed"
    else
        @error "Orchestrator was already paused"
    end
end

function run_orchestrator(;verbosity=0)
    mkpath(config_path())
    logger = RemoteHPC.TimestampLogger(RemoteHPC.TeeLogger(RemoteHPC.NotHTTPLogger(RemoteHPC.TimeBufferedFileLogger(config_path("log.log"), interval=0.0001)),
                                                           RemoteHPC.HTTPLogger(RemoteHPC.TimeBufferedFileLogger(config_path("HTTP.log"), interval=0.0001))))
    with_logger(logger) do
        LoggingExtras.withlevel(LoggingExtras.Debug; verbosity=verbosity) do
        try
            # RemoteREPL
            portpath = config_path("port")
            port, server = listenany(Sockets.localhost, 27754)
            write(portpath, "$port")
            @debug "Started Server at $port"
            t = Threads.@spawn serve_repl(server, on_client_connect = sess -> sess.in_module = RomeoDFT)

            o = ORCHESTRATOR[] = SearchOrchestrator(Dict{String, Searcher}(), false, false)
            
            loaded_searchers_path = config_path("loaded_searchers.txt")
            if ispath(loaded_searchers_path)
                rootdirs = readlines(loaded_searchers_path)
                @debug "Found $(length(rootdirs)) previously loaded Searchers."
                for r in rootdirs
                    try
                        n = load_searcher(r)
                        start_searcher(n)
                        sleep(1)
                    catch e
                        RemoteHPC.log_error(e)
                    end
                end
            end
            
            while !o.stop
                ls = collect(values(o.searchers))
                running_ls = filter(x -> x.loop !== nothing && !istaskdone(x.loop), ls)
                crashed_ls = filter(x -> x.loop !== nothing && istaskfailed(x.loop), ls)
                
                msg = "$(length(ls)) loaded searchers; $(length(running_ls)) running"
                if !isempty(crashed_ls)
                    ns = map(x->simname(x), crashed_ls)
                    msg *= "; $(length(crashed_ls)) crashed: $(join(ns, " - "))"
                end
                write(loaded_searchers_path, join(map(x -> x.rootdir, ls), "\n"))
                @debugv 2 msg
                unload_finished()
                sleep(1)
            end
            
            @debug "Stopping Server."
            stop_searchers(o)
            @debug "K bye then"
        catch e
            RemoteHPC.log_error(e)
        end
    end
    end
end
