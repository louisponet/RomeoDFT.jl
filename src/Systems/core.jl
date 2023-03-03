function add_calc!(job, settings)
    calc = gencalc(job, settings)
    calc.run = true
    cid = findfirst(x -> x.name == calc.name, job.calculations)
    if cid !== nothing
        job.calculations[cid] = calc
    else
        push!(job, calc)
    end
end

"""
    JobCreator

Takes care of creating jobs from [`SCFSettings`](@ref), [`NSCFSettings`](@ref), [`BandsSettings`](@ref), [`ProjwfcSettings`](@ref), and new [`Trials`](@ref Trial) added to a [`Simulation`](@ref). It creates a [`SimJob`](@ref) component for each and a [`Submit`](@ref) component to signal that it
should be submitted.
"""
struct JobCreator <: System end
function Overseer.requested_components(::JobCreator)
    return (SimJob, Submit, Template, SCFSettings, NSCFSettings, BandsSettings,
            BandsResults,
            ProjwfcSettings, Simulation, Trial, Generation, Results, WannierResults,
            WannierSettings, TimingInfo, BaseCase)
end

function Overseer.update(::JobCreator, m::AbstractLedger)
    @debugv 2 "[START] JobCreator"
    sinfo = m[SimJob]
    # Firefly
    sublock = ReentrantLock()
    max_new = sum(x -> Server(x.server).max_concurrent_jobs, m[ServerInfo].data) - length(@entities_in(m, SimJob && !Error))
    tot_new = 0
    @sync for e in @entities_in(m[Simulation] && m[Trial] && m[Generation] && !sinfo &&
                                !m[Results])
        if tot_new > max_new
            break
        elseif e.generation == e.current_generation
            tot_new += 1
            Threads.@spawn begin
                simn = simname(m)
                loc_dir = local_dir(m, e)
                if ispath(loc_dir)
                    rm(loc_dir; recursive = true)
                end

                if Hybrid in m && e in m[Hybrid]
                    c1 = deepcopy(e.template_calculation)
                    c2 = deepcopy(e.template_calculation)
                    suppress() do
                        c1[:system][:Hubbard_conv_thr] = 1e-5
                        c1[:system][:Hubbard_maxstep] = 100000
                        c1[:restart_mode] = "from_scratch"
                        set_name!(c1, "scf_1")

                        c2[:restart_mode] = "restart"
                    end
                        
                    calcs = Calculation[c1, c2]
                else
                    calcs = Calculation[deepcopy(e.template_calculation)]
                end
                suppress() do
                    set_name!(calcs[end], "scf")
                end
                j = Job("$(simn)_gen_$(e.current_generation)_id_$(Entity(e).id)",
                        deepcopy(e.template_structure),
                        calcs;
                        dir = loc_dir, server = local_server().name,
                        environment = "default")
                set_flow!(j, "" => true)
                suppress() do
                    j.calculations[1][:system][:Hubbard_occupations] = generate_Hubbard_occupations(e[Trial].state, j.structure)
                    delete!(j.calculations[1], :disk_io)
                end
                lock(sublock)
                m[e] = SimJob(loc_dir, "", j)
                m[e] = Submit()
                unlock(sublock)
            end
        end
    end

    # initial job
    @sync for e in @entities_in(m, Template && (BaseCase || SCFSettings || Trial) && !SimJob && !Results)
        if tot_new > max_new && e ∉ m[BaseCase]
            continue
        end
        tot_new += 1
        Threads.@spawn begin
            simn = simname(m)
            loc_dir = local_dir(m, e)
            tc = deepcopy(e.calculation)
            # Some setup here that's required
            delete!(tc, :disk_io)
             
            suppress() do
                set_name!(tc, "scf")
                if e in m[Trial]
                    tc[:system][:Hubbard_occupations] = generate_Hubbard_occupations(m[Trial][e].state, e.structure)
                end
                if e in m[SCFSettings]
                    for (f, v) in e.replacement_flags
                        if v isa Dict
                            for (f2, v2) in v
                                tc[f][f2] = v2
                            end
                        else
                            tc[f] = v
                        end
                    end
                elseif RelaxSettings in m && e in m[RelaxSettings]
                    relset = m[RelaxSettings][e]
                    tc[:forc_conv_thr] = relset.force_convergence_threshold
                    tc[:etot_conv_thr] = relset.energy_convergence_threshold
                    tc[:ion_dynamics] = relset.ion_dynamics 
                    tc[:cell_dynamics] = relset.cell_dynamics
                    if !relset.symmetry
                        tc[:nosym] = true
                    end
                    if relset.variable_cell 
                        set_name!(tc, "vcrelax")
                        tc[:calculation] = "vc-relax"
                    else
                        set_name!(tc, "relax")
                        tc[:calculation] = "relax"
                    end
                end
                if e in m[BaseCase]
                    delete_Hubbard!(tc)
                end
            end
            if Hybrid in m && e in m[Hybrid]
                c1 = deepcopy(tc) 
                c2 = deepcopy(tc)
                suppress() do
                    c1[:system][:Hubbard_conv_thr] = 1e-9
                    c1[:system][:Hubbard_maxstep] = 100000
                    c1[:restart_mode] = "from_scratch"
                    set_name!(c1, "scf_1")

                    c2[:restart_mode] = "restart"
                end
                calcs = Calculation[c1, c2]
            else
                calcs = Calculation[tc]
            end
            job = Job("$(simn)_id_$(Entity(e).id)",
                      deepcopy(e.structure),
                      calcs;
                      dir = loc_dir, server = local_server().name, environment = "default")
            set_flow!(job, "" => true)

            #making sure that the scf calculations are always named scf

            lock(sublock)
            m[e] = SimJob(loc_dir, "", job)
            m[e] = Submit()
            unlock(sublock)
        end
    end

    # BandSubmission
    bsettings = m[BandsSettings]
    nsettings = m[NSCFSettings]
    psettings = m[ProjwfcSettings]
    @sync for e in @entities_in(sinfo && !m[Simulation])
        Threads.@spawn begin
            if e in bsettings && !any(x -> x.name == "bands", e.job.calculations)
                suppress() do
                    add_calc!(e.job, bsettings[e])
                end
                lock(sublock)
                m[e] = Submit()
                unlock(sublock)
            end
            if e in nsettings && !any(x -> x.name == "nscf", e.job.calculations)
                suppress() do
                    add_calc!(e.job, nsettings[e])
                end
                lock(sublock)
                m[e] = Submit()
                unlock(sublock)
            end
            if e in psettings && !any(x -> x.name == "projwfc", e.job.calculations)
                projwfc = deepcopy(psettings[e])
                if e in m[Results]
                    projwfc.Emin = m[Results][e].fermi - projwfc.Emin
                    projwfc.Emax = m[Results][e].fermi + projwfc.Emax
                end
                suppress() do
                    add_calc!(e.job, psettings[e])
                end
                lock(sublock)
                m[e] = Submit()
                unlock(sublock)
            end
        end
    end

    @sync for e in @entities_in(sinfo && m[WannierSettings] && m[BandsResults] &&
                                !m[WannierResults] && !m[Submit] && !m[Error])
        Threads.@spawn if !any(x -> occursin("wan", x.name), e.job.calculations) &&
                          state(e.job) == RemoteHPC.Completed &&
                          any(x -> x.name == "projwfc", e.job.calculations) &&
                          ispath(e.job, splitdir(e.job["projwfc"].outfile)[end])
            for (el, projs) in e.projections
                for a in e.job.structure[element(el)]
                    a.projections = projs
                end
            end

            wans = suppress() do
                try
                    Calculations.gencalc_wan(e.job, e.dos_ratio)
                catch err
                    m[e] = Error("Error generating wannier inputs.")
                    return nothing
                end
            end
            if wans !== nothing
                for w in wans
                    id = findfirst(x -> x.name == w.name, e.job.calculations)
                    if id !== nothing
                        e.job.calculations[id] = w
                    else
                        push!(e.job, w)
                    end
                end

                Jobs.set_flow!(e.job, "" => false, "wan" => true)
                suppress() do
                    e.job[:dis_num_iter] = 1000
                    return e.job[:num_iter] = 3000
                end
                e.created = true
                lock(sublock)
                m[e] = Submit()
                unlock(sublock)
            end
        end
    end
    @debugv 2 "[STOP] JobCreator"
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

function Overseer.update(::JobSubmitter, m::AbstractLedger)
    @debugv 2 "[START] JobSubmitter"
    # We find the server to submit the next batch to by finding the one
    # with the pool with the least entities
    if !isalive(local_server())
        @warn "Local server not alive, not submitting"
        @debugv 2 "[STOP] JobSubmitter"
        return
    end
    sinfo = m[ServerInfo]
    jinfo = m[SimJob]
    sub = m[Submit]
    lck = ReentrantLock()
    function lock_(f)
        lock(lck)
        try
            f()
        finally
            unlock(lck)
        end
    end
    
    
    if !isempty(m[Simulation]) || !isempty(m[Template])
        if isempty(m[Simulation])
            # This can by definition only happen if we're restarting from previous finished run, so a template structure should exist
            str = m[Template][1].structure
            calc = m[Template][1].calculation
        else
            str = m[Simulation][1].template_structure
            calc = m[Simulation][1].template_calculation
        end
        
        for s in [map(x -> x.server, m[ServerInfo].data); local_server().name]
            server = Server(s)
            pseudodir = joinpath(server, m, "pseudos")
            if !isalive(server) || ispath(server, pseudodir)
                continue
            end
            mkpath(server, pseudodir)
            for a in str.atoms
                ppath = joinpath(pseudodir, "$(a.element.symbol).UPF")
                if !ispath(server, ppath)
                    pseudo = a.pseudo
                    pserver = Server(pseudo.server)
                    if !isempty(pseudo.path) && isalive(pserver)
                        write(server, ppath, read(pserver, pseudo.path))
                    elseif !isempty(pseudo.pseudo)
                        write(server, ppath, pseudo.pseudo)
                    else
                        error("Some pseudos could not be saved, set them first")
                    end
                end
            end
        end
    end

    # TODO: For now only 1 server
    server_info = sinfo[1]
    server_entity = entity(sinfo, 1)
    server = Server(server_info.server)
    
    if !isalive(server)
        @debugv 2 "[STOP] JobSubmitter"
        return
    end

    # First we check whether all previously submitted simjobs are in a pending or running state
    ready_to_submit = true
    @sync for e in @entities_in(m, SimJob && !Submit && !Done)
        if e.job.server == server.name
            Threads.@spawn begin
                if state(e.job) ∉ (RemoteHPC.Running, RemoteHPC.Pending, RemoteHPC.Completed, RemoteHPC.Failed, RemoteHPC.Completing)
                    ready_to_submit = false
                end
            end
        end
    end
    if !ready_to_submit
        return
    end
                
    to_resubmit_due_to_server = Entity[]
    entities_to_submit = sort(collect(@entities_in(sub)), by = x->x.id)[1:min(length(sub), 20)]
    @sync for e in entities_to_submit
        if !(e in m[SimJob])
            continue
        end
        Threads.@spawn begin
            lock_() do
                sinfo[e] = server_entity
            end
            j = jinfo[e].job
            j.server = server_info.server
            if !(state(server, jinfo[e].remote_dir) in (RemoteHPC.Running, RemoteHPC.Pending, RemoteHPC.Submitted))
                try
                    # So we don't have to always do abspath(server, remote_dir)
                    jinfo[e].remote_dir = joinpath(server, m, e)

                    for c in filter(x -> x.run, j.calculations)
                        if eltype(c) == QE
                            if c.exec.exec == "pw.x"
                                c.exec = load(server, Exec(server_info.pw_exec))
                            elseif c.exec.exec == "projwfc.x"
                                c.exec = load(server, Exec(server_info.pw_exec))
                                c.exec.exec = "projwfc.x"
                            end
                        elseif eltype(c) == Wannier90
                            c.exec = load(server, Exec("wannier90"))
                        end
                    end
                    firstid = findfirst(x -> x.run, j.calculations)

                    cid = findfirst(x -> occursin("wan", x.name) && x.run, j.calculations)
                    j.environment = server_info.environment

                    set_server_pseudos!(j, local_server(), m)
                    local_save(j, jinfo[e].local_dir)
                    j.dir = jinfo[e].remote_dir
                

                    j.server = server_info.server
                    set_server_pseudos!(j, server, m)

                    suppress() do
                        return submit(j; fillexecs = false, versioncheck=false, priority = e in m[NSCFSettings] || e in m[BaseCase] ? server_info.priority + 1 : server_info.priority)
                    end
                    curt = Dates.datetime2unix(now())
                    if e ∉ m[TimingInfo]
                        lock_() do 
                            m[e] = TimingInfo(curt, curt, 0.0, 0.0, 0.0, 0.0, "", 0.0, 0.0)
                        end
                    end
                catch err
                    lock_() do
                        m[e] = Error("Submission error of with:\n$err\n$(stacktrace(catch_backtrace()))")
                    end
                end
            end
        end
    end
    for e in entities_to_submit
        pop!(sub, e)
    end
    @debugv 2 "[STOP] JobSubmitter"
end

"""
    JobMonitor

Monitors the progress of running entities, i.e. with [`SimJob`](@ref) and [`TimingInfo`](@ref) components.
It updates the timing information in [`TimingInfo`](@ref) and potentially aborts and resubmits stalled jobs.
"""
struct JobMonitor <: System end

isparseable(s::RemoteHPC.JobState) = s ∈ (RemoteHPC.Completed, RemoteHPC.Failed, RemoteHPC.Completing, RemoteHPC.Timeout)
# isparseable(s::RemoteHPC.JobState) = true

function Overseer.update(::JobMonitor, m::AbstractLedger)
    @debugv 2 "[START] JobMonitor"
    lck = ReentrantLock()
    function lock_(f)
        lock(lck)
        try
            f()
        finally
            unlock(lck)
        end
    end
    @sync for e in @entities_in(m[SimJob] && m[TimingInfo] && !m[Done] && !m[Submit])
        Threads.@spawn if any(x -> x.run, e.job.calculations) 
            server = Server(e.job.server)
            curt   = Dates.datetime2unix(now())
            prevt  = e.cur_time
            dt     = curt - prevt
            s      = state(e.job)
            if s == RemoteHPC.Pending
                e.pending += dt
            elseif s == RemoteHPC.Running
                e.running += dt
                cur_running = Client.last_running_calculation(e.job).name
                cur_filesize = e.current_filesize 
                e.current_filesize = Float64(filesize(server, joinpath(e.remote_dir, e.job[cur_running].outfile)))
                if cur_running != e.current_running || e.current_filesize != cur_filesize
                    e.current_runtime = 0.0
                else
                    e.current_runtime += dt
                end
                e.current_running = cur_running
                success_id = findfirst(x->x.converged, m[Results])
                run_check_time = success_id !== nothing ? m[TimingInfo][entity(m[Results], success_id)].running : 1800
                run_check_time = run_check_time == 0 ? 1800 : run_check_time
                # If filesize didn't change for 30min we abort
                if e.current_runtime > run_check_time &&  cur_filesize == e.current_filesize
                    if length(DFC.versions(e.job)) < 3
                        abort(e.job)
                        for c in e.job.calculations
                            if c.name == cur_running
                                break
                            end
                            c.run = false
                        end
                        @debug "Aborted and resubmitted job with entity $(e.e).\nDuring: $(Client.last_running_calculation(e.job).name)."
                        lock_() do 
                            m[e] = Submit()
                        end
                    else
                        lock_() do
                            m[e] = Error("Job stalled.")
                        end
                    end
                end
            elseif s == RemoteHPC.Failed && ispath(server, e.remote_dir) &&
                   filesize(server,
                            joinpath(e.remote_dir,
                                     Client.last_running_calculation(e.job).outfile)) == 0.0
                if length(DFC.versions(e.job)) < 3
                    lock_() do 
                        m[e] = Submit()
                    end
                else
                    lock_() do
                        m[e] = Error("Job state Failed.")
                    end
                end
            elseif s == RemoteHPC.Unknown
                if ispath(server, e.remote_dir)
                    try
                        res = outputdata(e.job)
                        for (k, v) in res
                            e.job[k].run = !(get(v, :finished, false) || length(v) > 2)
                        end
                        @debug "Something happened with $(e.e)... Resubmitting to run: $(join(map(x->x.name, filter(y -> y.run, e.job.calculations)), " - "))"
                        lock_() do
                            m[e] = Submit()
                        end
                    catch err
                        lock_() do
                            m[e] = Error(err, "Something went wrong while pulling outputdata in JobMonitor.")
                        end
                    end
                else
                    lock_() do
                        m[e] = Submit()
                    end
                end
            elseif s in (RemoteHPC.NodeFail, RemoteHPC.Cancelled)
                lock_() do
                    m[e] = Submit()
                end
            end
            e.cur_time = curt
        end
    end
    @debugv 2 "[STOP] JobMonitor"
end

"""
    Cleaner

Removes temporary files of entities with a finished [`SimJob`](@ref).
"""
struct Cleaner <: System end
Overseer.requested_components(::Cleaner) = (Done, SimJob)

function Overseer.update(::Cleaner, m::AbstractLedger)
    @debugv 2 "[START] Cleaner"
    @sync for e in @entities_in(m, Done && SimJob)
        if !e.cleaned
            Threads.@spawn begin
                s = Server(e.job.server)
                jdir = e.remote_dir
                if e.remote_dir != e.local_dir && ispath(s, e.remote_dir)
                    @debugv 2 "Cleaning $jdir on Server $(s.name)"
                    rm(s, jdir)
                end
                m[Done][e] = Done(true)
            end
        end
    end
    for e in @safe_entities_in(m, SimJob && Done)
        if e.cleaned
            pop!(m[SimJob], e)
        end
    end
    # for e in @safe_entities_in(m, SimJob && BandsResults && Results)
    #     pop!(jinfo, e)
    # end
    # for e in @safe_entities_in(m, SimJob && !BandsSettings && Results)
    #     pop!(jinfo, e)
    # end
    @debugv 2 "[STOP] Cleaner"
end

struct OutputPuller <: System end
Overseer.requested_components(::OutputPuller) = (Done, SimJob)
function Overseer.update(::OutputPuller, m::AbstractLedger)
    @debugv 2 "[START] OutputPuller"
    jinfo = m[SimJob]
    sinfo = m[ServerInfo]
    tinfo = m[TimingInfo]
    lck = ReentrantLock()
    @sync for e in @entities_in(m, SimJob && !Submit && ServerInfo && TimingInfo && !Error)
        Threads.@spawn begin
            server = Server(e.job.server)
            if isparseable(DFC.state(e.job)) && ispath(server, abspath(server, e.remote_dir))
                apath = abspath(server, e.remote_dir)
                files = readdir(server, apath)
                running_calcs = filter(x -> x.run, e.job.calculations)
                curt = Dates.datetime2unix(now())
                for c in running_calcs

                    rem_path = joinpath(apath, c.outfile)
                    if ispath(server, rem_path)
                        loc_path = joinpath(e.local_dir, c.outfile)
                        if filesize(server, rem_path) < 100e6
                            write(loc_path,
                                  read(server, rem_path))
                        else
                            RemoteHPC.pull(server, rem_path, loc_path)
                        end
                        fs = filesize(loc_path)
                        if fs == 0
                            @debugv 2 "$(c.name) output file was empty, resubmitting."
                            c.run=true
                            lock(lck)
                            try
                                m[e] = Submit()
                            finally
                                unlock(lck)
                            end
                        elseif c.name == "projwfc"
                            for f in filter(x->occursin("pdos", x), files)
                                write(joinpath(e.local_dir, f), read(server, joinpath(apath, f)))
                            end
                            lock(lck)
                            try
                                m[e] = Done(false)
                            finally
                                unlock(lck)
                            end
                        else
                            c.run = false
                        end
                    end
                end
                for f in filter(x->occursin("slurm", x), files)
                    write(joinpath(e.local_dir, f), read(server, joinpath(apath, f)))
                end
                e.postprocessing += Dates.datetime2unix(now()) - curt
            end
        end
    end
    @debugv 2 "[STOP] OutputPuller"
end

struct ErrorCorrector <: System end

function Overseer.update(::ErrorCorrector, m::AbstractLedger)
    @debugv 2 "[START] ErrorCorrector"
    for e in @safe_entities_in(m, Results)
        if length(e.state.occupations) == 0
            @debugv 1 "ErrorCorrector: $(e.e) has an empty State in Results"
            # This results usually because of some server side issue where something crashed for no reason
            # Can also be because things converged before constraints were released (see process_Hubbard)
            e in m[Error] && pop!(m[Error], e)
            e in m[Done] && pop!(m[Done], e)
            e in m[BandsResults] && pop!(m[BandsResults], e)
            if e in m[SimJob]
                if e.niterations == e.constraining_steps
                    @debugv 1 "ErrorCorrector: $(e.e) had converged scf while constraints still applied, increasing Hubbard_conv_thr."
                    pop!(m[ServerInfo], e)
                    pop!(m[TimingInfo], e)
                    m[Template][e].calculation[:system][:Hubbard_conv_thr] = 1.5 * m[SimJob][e].job.calculations[1][:system][:Hubbard_conv_thr]
                    pop!(m[SimJob], e)
                else
                    set_flow!(m[SimJob][e].job, "" => true)
                    m[Submit][e] = Submit()
                end
            else
                m[Template][e] = m[Template][end]
            end
            pop!(m[Results], e)
        end
    end
    for e in @entities_in(m, SimJob && !Results && !Submit)
        if !any(x-> x.run, e.job.calculations) || e ∉ m[TimingInfo]
            e in m[Error] && pop!(m[Error], e)
            @debugv 1 "ErrorCorrector: resubmitting $(e.e)"
            set_flow!(e.job, "" => true)
            m[e] = Submit()
        end
    end

    # To add BaseCase if there was already a simulation running
    # if !(BaseCase in m) || (isempty(m[BaseCase]) && !isempty(m[Simulation]) && maximum(x -> x.current_generation >= 0, m[Simulation]))
    #     sime = entity(m[Simulation], 1)
    #     e = Entity(m, BaseCase(), Template(m[Simulation][sime].template_structure, m[Simulation][sime].template_calculation))
    #     if Hybrid in m && sime in m[Hybrid]
    #         m[e] = Hybrid()
    #     end
    # end
    @debugv 2 "[STOP] ErrorCorrector"
end

struct Stopper <: System end
Overseer.requested_components(::Stopper) = (Done, SimJob, IntersectionSearcher, Simulation, StopCondition)

function stop_check(maxgen::Int, m::AbstractLedger)
    stop_condition = singleton(m, StopCondition)
    if maxgen < stop_condition.n_generations
        return false, Int[], Int[]
    end
    n_unique = zeros(Int, maxgen)
    n_total = zeros(Int, maxgen)
    for e in @entities_in(m, Results && Generation)
        if e in m[Unique] && e.generation != 0
            n_unique[e.generation] += 1
        end
        if e.generation != 0
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

stop_check(m::AbstractLedger) = stop_check(maximum(x->x.generation, m[Generation], init=0), m)

# Check if BaseCase was ran with the magnetizations of the minimum state
function check_basecase(m::AbstractLedger)
    base_e         = entity(m[BaseCase], length(m[BaseCase]))
    str            = m[Template][base_e].structure
    magats         = filter(x -> sum(x.magnetization) != 0 && x.dftu.U != 0, str.atoms)
    curmags        = map(x -> x.magnetization[3], magats)
    es             = collect(@entities_in(m, Unique && Results))
    min, minid     = findmin(x -> x.total_energy, es)
    magnetizations = map(x -> abs(x) < 1e-2 ? 1e-5 : sign(x), es[minid].state.magmoms)[1:length(magats)]
    if any(!iszero, curmags .- magnetizations) && any(!iszero, curmags .+ magnetizations)
        @debug "Postprocessing finished, running \"vanilla\" QE with starting magnetizations of global minimum."
        str = deepcopy(str)
        magats = filter(x -> sum(x.magnetization) != 0 && x.dftu.U != 0, str.atoms)
        for (mag, at) in zip(magnetizations, magats)
            at.magnetization = [0, 0, mag]
        end
        Entity(m, BaseCase(), Template(str, deepcopy(m[Template][base_e].calculation)), Generation(maximum(x->x.generation, m[Generation], init=0)))
        return false
    else
        return true
    end
end

function Overseer.update(::Stopper, m::AbstractLedger)

    # Handle cleanup stopping
    if m.mode == :cleanup
        if length(@entities_in(m, SimJob && !Submit && !Error)) == 0
            @debug "All cleanup has finished"
            m.stop = true
            m.finished = true
        end
        return
    end

    # Handle search stopping
    if isempty(m[StopCondition])
        @warn "No StopCondition was set."
    end

    maxgen =  maximum(x->x.generation, m[Generation], init=0)
    stop, n_unique, n_total = stop_check(maxgen, m)
    if isempty(n_unique)
        set_mode!(m, :search)
        prepare(m)
        return
    end
        
    if m.ledger.stages[1].name == :core
        if length(@entities_in(m[SimJob] && !m[Error])) == 0 && all(x->x.cleaned, m[Done]) && all(x -> x in m[Done] || x in m[Error], @entities_in(m[Trial])) 
            # Check if stop condition is still met
            if stop && check_basecase(m)
                @debug "All postprocessing has finished, stopping search."
                @debug "Found $(sum(n_unique)) Unique states after $(sum(n_total)) Trials."
                m.stop = true
                m.finished = true
            else
                @debug "Stop condition not met after postprocessing at Generation($maxgen). Continuing search."
                set_mode!(m, :search)
                prepare(m)
            end
        end
    elseif stop
        set_mode!(m, :postprocess)
        prepare(m)
        @debug "Stop condition was met at Generation($maxgen), switching to pure postprocessing mode."
    end
    @debugv 2 "Generation($maxgen): $(sum(n_unique)) unique states after $(sum(n_total)) trials." 
end
