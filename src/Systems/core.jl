"""
    JobCreator

Takes care of creating jobs from [`SCFSettings`](@ref), new [`Trials`](@ref Trial), or [`BaseCase`](@ref).
It creates a [`SimJob`](@ref) component for each and a [`Submit`](@ref) component to signal that it
should be submitted.
"""
struct JobCreator <: System end
function Overseer.requested_components(::JobCreator)
    return (SimJob, Submit, Template, SCFSettings, NSCFSettings, BandsSettings,
            BandsResults, Hybrid,
            ProjwfcSettings, Simulation, Trial, Generation, Results, WannierResults,
            WannierSettings, TimingInfo, BaseCase, Submitted, Completed, Running)
end

function Overseer.update(::JobCreator, m::AbstractLedger)
    info = m[SearcherInfo][1]
    max_new = info.max_concurrent_trials - info.n_running_calcs

    tot_new = 0
    lck = ReentrantLock()
    info.n_pending_calcs = length(m[Submit])
    # initial job
    @error_capturing_threaded for e in @entities_in(m,
                                                    Template &&
                                                    (BaseCase || SCFSettings || Trial) &&
                                                    !SimJob && !Done)
        if tot_new > max_new && e ∉ m[BaseCase]
            lock(lck) do
                info.n_pending_calcs += 1
            end
            continue
        else
            lock(lck) do
                tot_new += 1
            end
        end

        scf_calc = deepcopy(e.calculation)
        # Some setup here that's required
        scf_calc.run = true
        delete!(scf_calc, :disk_io)

        suppress() do
            set_name!(scf_calc, "scf")

            if e in m[Trial]
                scf_calc[:system][:Hubbard_occupations] = generate_Hubbard_occupations(m[Trial][e].state,
                                                                                       e.structure)
            end

            scf_calc[:system][:Hubbard_conv_thr] = 1e-12
            if e in m[Intersection]
                scf_calc[:system][:Hubbard_maxstep] = 10000
            end

            if e in m[SCFSettings]
                for (f, v) in e.replacement_flags
                    if v isa Dict
                        for (f2, v2) in v
                            scf_calc[f][f2] = v2
                        end
                    else
                        scf_calc[f] = v
                    end
                end
            end

            if e in m[BaseCase]
                delete_Hubbard!(scf_calc)
            end
        end
        if e in m[Hybrid]
            c1 = deepcopy(scf_calc)
            c2 = deepcopy(scf_calc)
            suppress() do
                c1[:system][:Hubbard_conv_thr] = 1e-9
                c1[:system][:Hubbard_maxstep] = 100000
                c1[:restart_mode] = "from_scratch"
                set_name!(c1, "scf_1")

                return c2[:restart_mode] = "restart"
            end
            calcs = Calculation[c1, c2]
        else
            calcs = Calculation[scf_calc]
        end
        job = Job("$(simname(m))_id_$(Entity(e).id)",
                  deepcopy(e.structure),
                  calcs;
                  dir = local_dir(m, e), server = local_server().name,
                  environment = "default")

        m[e] = SimJob(job.dir, "", job)
        set_status!(m, e, Submit())
    end
end

function ensure_pseudos_uploaded!(m::AbstractLedger)
    isempty(m[Template]) && return

    # This can by definition only happen if we're restarting from previous finished run, so a template structure should exist
    str  = m[Template][1].structure
    calc = m[Template][1].calculation

    servers = all_servers(m)
    for s in servers
        pseudodir = joinpath(s, m, "pseudos")
        if !isalive(s) || ispath(s, pseudodir)
            continue
        end

        mkpath(s, pseudodir)

        for a in str.atoms
            ppath = joinpath(pseudodir, "$(a.element.symbol).UPF")

            if ispath(s, ppath)
                continue
            end

            pseudo = a.pseudo
            pserver = Server(pseudo.server)

            if !isempty(pseudo.path) && isalive(pserver)
                write(s, ppath, read(pserver, pseudo.path))

            elseif !isempty(pseudo.pseudo)
                write(s, ppath, pseudo.pseudo)

            else
                error("Some pseudos could not be saved, set them first")
            end
        end
    end
end

"""
    JobSubmitter

Submits entities with [`SimJob`](@ref) and [`Submit`](@ref) components.
"""
struct JobSubmitter <: System end
function set_server_pseudos!(j::Job, s::Server, m::Searcher)
    for a in j.structure.atoms
        a.pseudo.path = joinpath(s, m, "pseudos/$(a.element.symbol).UPF")
        a.pseudo.server = s.name
    end
end

Overseer.requested_components(::JobSubmitter) = (Running, Submitted)

function Overseer.update(::JobSubmitter, m::AbstractLedger)
    info = m[SearcherInfo][1]
    # We find the server to submit the next batch to by finding the one
    # with the pool with the least entities
    sinfo = m[ServerInfo]
    ensure_pseudos_uploaded!(m)

    # TODO: For now only 1 server
    server_info = sinfo[1]
    server_entity = entity(sinfo, 1)
    server = Server(server_info.server)

    # First we check whether server is alive and
    # all previously submitted simjobs are in a pending or running state
    if !isempty(@entities_in(m, Submitted && !Error))
        return
    end

    submit_comp = m[Submit]
    pvec = sortperm(Overseer.indices(submit_comp).packed; rev = true)
    permute!(submit_comp, pvec)

    lck = ReentrantLock()

    @error_capturing_threaded for e in @safe_entities_in(m, Submit && SimJob)
        sinfo[e] = server_entity

        e.job.server = server_info.server
        already_submitted = state(server, e.remote_dir) in
                            (RemoteHPC.Running, RemoteHPC.Pending, RemoteHPC.Submitted)

        if already_submitted
            log(e, "was already submitted somehow...")
            if e ∉ m[TimingInfo]
                curt = Dates.datetime2unix(now())
                m[e] = TimingInfo(curt, curt, 0.0, 0.0, 0.0, 0.0, "", 0.0, 0.0)
            end
            set_status!(m, e, Running())
            continue
        end

        # So we don't have to always do abspath(server, remote_dir)
        e.remote_dir = joinpath(server, m, e)

        for c in filter(x -> x.run, e.job.calculations)
            if eltype(c) == QE
                ex      = load(server, Exec(server_info.pw_exec))
                ex.path = joinpath(dirname(ex.path), exec(c.exec))
                c.exec  = ex

            elseif eltype(c) == Wannier90
                c.exec = load(server, Exec("wannier90"))
            end
        end

        e.job.environment = server_info.environment

        set_server_pseudos!(e.job, local_server(), m)
        local_save(e.job, e.local_dir)
        e.job.dir = e.remote_dir

        set_server_pseudos!(e.job, server, m)

        suppress() do
            priority = e in m[NSCFSettings] || e in m[BaseCase] ? server_info.priority + 1 :
                       server_info.priority
            submit(e.job; fillexecs = false, versioncheck = false,
                          priority = priority)
        end

        if e ∉ m[TimingInfo]
            curt = Dates.datetime2unix(now())
            m[e] = TimingInfo(curt, curt, 0.0, 0.0, 0.0, 0.0, "", 0.0, 0.0)
        end

        set_status!(m, e, Submitted())
        lock(lck) do
            info.n_running_calcs += 1
        end
    end
end

"""
    JobMonitor

Monitors the progress of running entities, i.e. with [`SimJob`](@ref) and [`TimingInfo`](@ref) components.
It updates the timing information in [`TimingInfo`](@ref) and potentially aborts and resubmits stalled jobs.
"""
struct JobMonitor <: System end

function isparseable(s::RemoteHPC.JobState)
    return s ∈
           (RemoteHPC.Completed, RemoteHPC.Failed, RemoteHPC.Completing, RemoteHPC.Timeout)
end

function Overseer.update(::JobMonitor, m::AbstractLedger)
    run_check_time = average_runtime(m)
    run_check_time = run_check_time == 0 ? 1800 : run_check_time

    function maybe_rerun(e, logmsg)
        if e in m[Rerun] && m[Rerun][e].count >= 3
            m[e] = Error(e, "Reran already 3 times and still no dice...")
        else
            log(e, logmsg)
            should_rerun(m, e)
        end
    end

    lck = ReentrantLock()
    @error_capturing_threaded for e in @safe_entities_in(m,
                                                         SimJob && TimingInfo && !Completed &&
                                                         !Submit && !Pulled)
        if !any(x -> x.run, e.job.calculations) && e ∉ m[Submitted]
            continue
        end

        server = Server(e.job.server)
        curt   = Dates.datetime2unix(now())
        prevt  = e.cur_time
        dt     = curt - prevt

        s = state(e.job)
        set_status!(m, e, s)

        if s == RemoteHPC.Pending
            e.pending += dt

        elseif s == RemoteHPC.Running
            e.running += dt

            cur_running   = Client.last_running_calculation(e.job).name
            prev_filesize = e.current_filesize
            cur_filesize  = filesize(server, joinpath(e.remote_dir, e.job[cur_running].outfile))

            e.current_filesize = Float64(cur_filesize)
            if cur_running != e.current_running || cur_filesize != prev_filesize
                e.current_runtime = 0.0
            else
                e.current_runtime += dt
            end

            e.current_running = cur_running

            # If filesize didn't change for 30min we abort
            if e.current_runtime > run_check_time && prev_filesize == e.current_filesize
                abort(e.job)
                for c in e.job.calculations
                    if c.name == cur_running
                        break
                    end
                    c.run = false
                end

                maybe_rerun(e,
                            "Aborted and resubmitted job during: $(Client.last_running_calculation(e.job).name).")
                m[SearcherInfo][1].n_running_calcs -= 1
            end

        elseif s == RemoteHPC.Failed
            ofile = joinpath(e.remote_dir, Client.last_running_calculation(e.job).outfile)
            if ispath(server, e.remote_dir) && filesize(server, ofile) == 0.0
                maybe_rerun(e, "Job failed.")
            end
            m[SearcherInfo][1].n_running_calcs -= 1

        elseif s in (RemoteHPC.NodeFail, RemoteHPC.Cancelled)
            maybe_rerun(e,
                        "Job Cancelled during: $(Client.last_running_calculation(e.job).name).")

        elseif isparseable(s)
            set_status!(m, e, Completed())
            lock(lck) do
                m[SearcherInfo][1].n_running_calcs -= 1
            end
        end
        e.cur_time = curt
    end
end

"""
    Cleaner

Removes temporary files of entities with a finished [`SimJob`](@ref).
"""
struct Cleaner <: System end
Overseer.requested_components(::Cleaner) = (Done, SimJob)

function Overseer.update(::Cleaner, m::AbstractLedger)
    @error_capturing_threaded for e in @safe_entities_in(m, Done && SimJob)
        if e.cleaned
            continue
        end
        s = Server(e.job.server)

        if e.remote_dir != e.local_dir && ispath(s, e.remote_dir)
            log(e, "Cleaning $(e.remote_dir) on Server $(s.name)")
            rm(s, e.remote_dir)
        end
        e[Done] = Done(true)
    end

    @error_capturing_threaded for e in @safe_entities_in(m, SimJob && Error)
        s = Server(e.job.server)

        opath = joinpath(e.remote_dir, "outputs")
        if ispath(s, opath)
            log(e, "Cleaning $opath on Server $(s.name)")
            rm(s, opath)
        end
    end
end

"""
    OutputPuller

Retrieves the remote files and stores them in the local jobdirs
"""
struct OutputPuller <: System end
Overseer.requested_components(::OutputPuller) = (Done, SimJob, Completed, Pulled)
function Overseer.update(::OutputPuller, m::AbstractLedger)
    @error_capturing_threaded for e in @safe_entities_in(m, SimJob && Completed)
        curt = Dates.datetime2unix(now())

        server = Server(e.job.server)

        should_resubmit = false

        running_calcs = filter(x -> x.run, e.job.calculations)
        RemoteHPC.pull(e.job, e.local_dir; calcs = running_calcs)

        for c in running_calcs
            ofile = joinpath(e.local_dir, c.outfile)
            if !ispath(ofile) || filesize(ofile) == 0
                log(e, "$(c.name) output file is missing or empty, resubmitting.")
                c.run = true
                should_resubmit = true
                ispath(ofile) && rm(ofile)
            else
                c.run = false
            end
        end

        if should_resubmit
            should_rerun(m, e)
        else
            set_status!(m, e, Pulled())
        end

        m[e].postprocessing += Dates.datetime2unix(now()) - curt
    end
end

"""
    SimJobRemover

Checks whether everything is done and if so removes the [`SimJob`](@ref) from an [`Entity`](@ref).
"""
struct SimJobRemover <: System end

function Overseer.update(::SimJobRemover, m::AbstractLedger)
    for e in @safe_entities_in(m, Done && SimJob)
        if e.cleaned
            pop!(m[SimJob], e)
        end
    end
end

"""
    Rerunner

Handles rerunning of jobs, removing the required [`Components`] as given in the [`ShouldRerun`](@ref) component.
"""
struct Rerunner <: System end

Overseer.requested_components(::Rerunner) = (ShouldRerun, Rerun)

function Overseer.update(::Rerunner, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, ShouldRerun)
        
        from_scratch = any(x -> x isa Template, e.components_to_replace)

        if from_scratch
            m[e] = Done(false)

            # TODO always assumes postprocessing and does not keep previous results
            create_postprocess_child!(m, e, e.components_to_replace...)
        else
            for d in e.components_to_replace
                m[e] = d
            end
            
            set_status!(m, e, Submit())
            
            trypop!(m[Done], e)
            
            m[e] = Rerun(e in m[Rerun] ? m[Rerun][e].count + 1 : 1)
        end

        pop!(m[ShouldRerun], e)
    end
end

"""
    ErrorCorrector

Tries to correct some common errors.
"""
struct ErrorCorrector <: System end

function Overseer.update(::ErrorCorrector, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, Results && !ShouldRerun && !Done)
        if length(e.state.occupations) == 0
            log(e, "ErrorCorrector: has an empty State in Results")
            # This results usually because of some server side issue where something crashed for no reason
            # Can also be because things converged before constraints were released (see process_Hubbard)
            if e in m[SimJob]
                set_flow!(m[SimJob][e].job, "" => true)
                should_rerun(m, e)
            end
        elseif e.constraining_steps == -1 && !(e in m[Intersection])
            log(e,
                "ErrorCorrector: has converged scf while constraints still applied, increasing Hubbard_conv_thr.")
            new_template = deepcopy(m[Template][e])
            new_template.calculation[:system][:Hubbard_conv_thr] *= 1.5
            should_rerun(m, e, new_template)
        end
    end

    # @error_capturing for e in @safe_entities_in(m, SimJob && !ShouldRerun && !Done)
    #     if !any(x-> x.run, e.job.calculations)
    #         log(e, "ErrorCorrector: Nothing happened during postprocessing stage, resubmitting")
    #         set_flow!(e.job, "" => true)
    #         should_rerun(m, e, BandsResults, Results, RelaxResults, FlatBands, Completed, Pulled)
    #     end
    # end

    for e in @safe_entities_in(m, Done && Error)
        if occursin("HTTP.Exceptions.StatusError(500, \"POST\", \"/rm/?", string(e.err))
            pop!(m[Error], e)
        end
    end
end

"""
    Stopper

Checks whether the [`StopCondition`](@ref) is met, and handles the automatic switching of
the mode of the [`Searcher`](@ref).
"""
struct Stopper <: System end
function Overseer.requested_components(::Stopper)
    return (Done, SimJob, IntersectionSearcher, Simulation, StopCondition)
end

function stop_check(maxgen::Int, m::AbstractLedger)
    stop_condition = singleton(m, StopCondition)
    if maxgen < stop_condition.n_generations
        return false, Int[], Int[]
    end
    n_unique = zeros(Int, maxgen)
    n_total  = zeros(Int, maxgen)
    for e in @safe_entities_in(m, Results && Generation && !Parent)
        if e in m[Unique] && e.generation != 0
            n_unique[e.generation] += 1
        end
        if e.generation != 0 && e in m[Trial] && m[Trial][e].origin != IntersectionMixed
            n_total[e.generation] += 1
        end
    end
    n_conseq = 0
    for i in 1:maxgen
        r = n_unique[i] / n_total[i]
        if r < stop_condition.unique_ratio
            n_conseq += 1
        elseif n_conseq >= stop_condition.n_generations
            break
        else
            n_conseq = 0
        end
    end
    @debugv 2 "Unique to trial ratio for the last $maxgen generations: $(n_unique./n_total)"
    return n_conseq >= stop_condition.n_generations, n_unique, n_total
end

function stop_check(m::AbstractLedger)
    return stop_check(maximum(x -> x.generation, m[Generation]; init = 0), m)
end

# Check if BaseCase was ran with the magnetizations of the minimum state
function check_basecase!(m::AbstractLedger)
    base_e = entity(m[BaseCase], length(m[BaseCase]))
    base_str = m[Template][base_e].structure

    magats = filter(ismagnetic, base_str.atoms)
    basecase_magmoms = map(x -> x.magnetization[3], magats)

    unique_es = collect(@entities_in(m, Unique && Results))

    if isempty(unique_es)
        return true
    end

    minid = findmin(x -> x.total_energy, unique_es)[2]

    minimum_magmoms = unique_es[minid].state.magmoms[1:length(magats)]
    minimum_zmag    = map(x -> abs(x) < 1e-2 ? 1e-5 : sign(x), minimum_magmoms)

    if any(!iszero, basecase_magmoms .- minimum_zmag) &&
       any(!iszero, basecase_magmoms .+ minimum_zmag)
        @debug "Running \"vanilla\" QE with starting magnetizations of global minimum."
        new_str = deepcopy(base_str)

        magats = filter(ismagnetic, new_str.atoms)
        for (mag, at) in zip(minimum_zmag, magats)
            at.magnetization = [0, 0, mag]
        end

        new_e = Entity(m, BaseCase(),
                       Template(new_str, deepcopy(m[Template][base_e].calculation)),
                       Generation(maximum_generation(m)))

        for c in Iterators.filter(x -> x isa PostProcessSettings, m[base_e])
            component = m[typeof(c)]
            if component isa PooledComponent
                component[new_e] = base_e
            else
                component[new_e] = deepcopy(c)
            end
        end
        return false
    else
        return true
    end
end

"""
    verify_groundstates!(m::AbstractLedger)

Makes sure that the global groundstates (including postprocessing) and
the ground state found from purely searching both hold the necessary `Components` like `FlatBands`
"""
function verify_groundstates!(m::AbstractLedger)
    gs_full   = ground_state(filter(x -> x.converged, @entities_in(m, Results)))
    gs_search = ground_state(filter(x -> x.converged, @entities_in(m, Results && !Parent)))

    for gs in (gs_full, gs_search)
        ldir = local_dir(m, gs)
        if gs ∉ m[FlatBands] && ispath(ldir)
            j = load(local_server(), Job(ldir))
            o = outputdata(j; calcs = [j.calculations[1].name])[j.calculations[1].name]
            m[gs] = FlatBands(flatbands(o))
        end
    end
end

function Overseer.update(::Stopper, m::AbstractLedger)

    # Handle cleanup stopping
    if m.mode == :cleanup
        entities_to_clean = @entities_in(m, SimJob && !Submit && !Error)
        if length(entities_to_clean) == 0
            @debug "All cleanup has finished"
            m.stop     = true
            m.finished = true
        end
        return
    end
    if m.mode == :manual
        return
    end

    # Handle search stopping
    if isempty(m[StopCondition])
        @warn "No StopCondition was set."
    end

    maxgen = maximum_generation(m)
    stop_condition_met, n_unique, n_total = stop_check(maxgen, m)

    search_entities = @entities_in(m, Trial && !Intersection)
    search_done     = all(x -> x in m[Done] || x in m[Error], search_entities)

    if (!stop_condition_met && search_done)
        set_mode!(m, :search)
        prepare(m)
        return
    end

    if mode(m) == :postprocess
        all_entities = @entities_in(m, (Trial || BaseCase) && !(Done || Error))
        if length(all_entities) == 0 && all(x->x.cleaned, @entities_in(m, Done && !Error))
            # Check if stop condition is still met
            if stop_condition_met && check_basecase!(m)
                verify_groundstates!(m)
                @debug "All postprocessing has finished, stopping search."
                @debug "Found $(sum(n_unique)) Unique states after $(sum(n_total)) Trials."

                try
                    # During testing this might fail
                    final_report(m)
                catch
                    nothing
                end

                m.stop = true
                m.finished = true
            else
                @debug "Stop condition not met after postprocessing at Generation($maxgen). Continuing search."
                set_mode!(m, :search)
                prepare(m)
            end
        end
    elseif stop_condition_met
        set_mode!(m, :postprocess)
        prepare(m)
        @debug "Stop condition was met at Generation($maxgen), switching to pure postprocessing mode."
    end
    @debugv 2 "Generation($maxgen): $(sum(n_unique)) unique states after $(sum(n_total)) trials."
end
