function orchestrator_port()
    pp = config_path("port")
    return parse(Int, read(pp, String))
end
function orchestrator_eval(str)
    try
        RemoteREPL.remote_eval(Sockets.localhost, orchestrator_port(), str)
    catch
        error("Orchestrator not running, start it first with `start_orchestrator()` or `romeo orchestrator start`.")
    end
end

function connect_orchestrator()
    try
        port = orchestrator_port()
        connect_repl(Sockets.localhost, port)
        orchestrator_eval("loaded_searchers()")
    catch
        error("Orchestrator not running, start it first with `start_orchestrator()`.")
    end
end

function start_orchestrator()
    try
        orchestrator_eval("loaded_searchers()")
    catch
        t = mtime(config_path("log.log"))
        julia_cmd = VERSION.minor >= 9 ? [joinpath(Sys.BINDIR, "julia"), "--heap-size-hint=2.5G"] : [joinpath(Sys.BINDIR, "julia")]
        
        r = run(Cmd(Cmd([julia_cmd..., "--project=$(config_path())", "--startup-file=no","-t", "auto", "-e",  "using Revise; using Plots; using LaTeXStrings; using UnicodePlots; using RomeoDFT; RomeoDFT.run_orchestrator()", "&>", config_path("log.log")]), detach=true, ignorestatus=true), wait=false)
        p = ProgressMeter.ProgressUnknown("Starting orchestrator. This takes a while the first time..."; spinner=true)

        while process_running(r) && mtime(config_path("log.log")) == t
            ProgressMeter.next!(p, spinner="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
            sleep(0.1)
        end
        ProgressMeter.finish!(p, spinner = process_running(r) ? '✓' : '✗')
    end
end

function client_stop_orchestrator()
    try
        id = parse(Int, orchestrator_eval("getpid()").content)
        orchestrator_eval("stop_orchestrator()")
        p = ProgressMeter.ProgressUnknown("Waiting for orchestrator shutdown..."; spinner=true)
        isalive = ccall(:uv_kill, Cint, (Cint, Cint), id, 0) == 0
        while isalive
            ProgressMeter.next!(p, spinner="⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏")
            isalive = ccall(:uv_kill, Cint, (Cint, Cint), id, 0) == 0
            sleep(0.1)
        end
        ProgressMeter.finish!(p)
    catch
        error("Orchestrator not running.")
    end
end

"""
    orchestrator_submit(searcher)

Save and submit `searcher` to the running [Orchestrator](@ref).
"""
function orchestrator_submit(l; verbosity=0)
    save(l)
    orchestrator_eval("start_searcher(\"$(l.rootdir)\"; verbosity=$(verbosity), sleep_time=Int($(l.sleep_time)))")
end
