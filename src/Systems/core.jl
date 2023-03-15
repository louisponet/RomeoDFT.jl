"""
    JobCreator

Takes care of creating jobs from [`SCFSettings`](@ref), [`NSCFSettings`](@ref), [`BandsSettings`](@ref), [`ProjwfcSettings`](@ref), and new [`Trials`](@ref Trial) added to a [`Simulation`](@ref). It creates a [`SimJob`](@ref) component for each and a [`Submit`](@ref) component to signal that it
should be submitted.
"""
struct JobCreator <: System end
function Overseer.requested_components(::JobCreator)
    return (SimJob, Submit, Template, SCFSettings, NSCFSettings, BandsSettings,
            BandsResults, Hybrid,
            ProjwfcSettings, Simulation, Trial, Generation, Results, WannierResults,
            WannierSettings, TimingInfo, BaseCase)
end

function Overseer.update(::JobCreator, m::AbstractLedger)
    sinfo = m[SimJob]
    # Firefly
    max_new = sum(x -> Server(x.server).max_concurrent_jobs, m[ServerInfo].data) - length(@entities_in(m, SimJob && !Error))
    tot_new = 0
    @error_capturing_threaded for e in @safe_entities_in(m, Simulation && Trial && Generation && !SimJob && !Results)
        if tot_new > max_new
            break
        elseif e.generation == e.current_generation
            tot_new += 1
            simn = simname(m)
            loc_dir = local_dir(m, e)
            if ispath(loc_dir)
                rm(loc_dir; recursive = true)
            end
            
            if e in m[Hybrid]
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
            m[e] = SimJob(loc_dir, "", j)
            set_state!(m, e, Submit())
        end
    end

    # initial job
    @error_capturing_threaded for e in @safe_entities_in(m, Template && (BaseCase || SCFSettings || Trial) && !SimJob && !Done)
        if tot_new > max_new && e ∉ m[BaseCase]
            continue
        end
        tot_new += 1
        simn    = simname(m)
        loc_dir = local_dir(m, e)
        tc      = deepcopy(e.calculation)
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
            end
            if e in m[BaseCase]
                delete_Hubbard!(tc)
            end
        end
        if e in m[Hybrid]
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

        m[e] = SimJob(loc_dir, "", job)
        set_state!(m, e, Submit())
    end

    # @sync for e in @entities_in(sinfo && m[WannierSettings] && m[BandsResults] &&
    #                             !m[WannierResults] && !m[Submit] && !m[Error])
    #     Threads.@spawn if !any(x -> occursin("wan", x.name), e.job.calculations) &&
    #                       state(e.job) == RemoteHPC.Completed &&
    #                       any(x -> x.name == "projwfc", e.job.calculations) &&
    #                       ispath(e.job, splitdir(e.job["projwfc"].outfile)[end])
    #         for (el, projs) in e.projections
    #             for a in e.job.structure[element(el)]
    #                 a.projections = projs
    #             end
    #         end

    #         wans = suppress() do
    #             try
    #                 Calculations.gencalc_wan(e.job, e.dos_ratio)
    #             catch err
    #                 m[e] = Error("Error generating wannier inputs.")
    #                 return nothing
    #             end
    #         end
    #         if wans !== nothing
    #             for w in wans
    #                 id = findfirst(x -> x.name == w.name, e.job.calculations)
    #                 if id !== nothing
    #                     e.job.calculations[id] = w
    #                 else
    #                     push!(e.job, w)
    #                 end
    #             end

    #             Jobs.set_flow!(e.job, "" => false, "wan" => true)
    #             suppress() do
    #                 e.job[:dis_num_iter] = 1000
    #                 return e.job[:num_iter] = 3000
    #             end
    #             e.created = true
    #             lock(sublock)
    #             m[e] = Submit()
    #             unlock(sublock)
    #         end
    #     end
    # end
end

function ensure_pseudos_uploaded!(m)
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
    # We find the server to submit the next batch to by finding the one
    # with the pool with the least entities
    if !isalive(local_server())
        @warn "Local server not alive, not submitting"
        return
    end
    sinfo = m[ServerInfo]
    jinfo = m[SimJob]
    sub = m[Submit]
    ensure_pseudos_uploaded!(m) 

    # TODO: For now only 1 server
    server_info = sinfo[1]
    server_entity = entity(sinfo, 1)
    server = Server(server_info.server)
    
    # First we check whether server is alive and
    # all previously submitted simjobs are in a pending or running state
    if !isalive(server) || !isempty(@entities_in(m, Submitted && !Error))
        return
    end

    pvec = sortperm(sub.indices.packed, rev=true)
    permute!(sub, pvec)
    @error_capturing_threaded for e in @safe_entities_in(m, Submit && SimJob)
        sinfo[e] = server_entity
        j = jinfo[e].job
        j.server = server_info.server
        already_submitted = state(server, jinfo[e].remote_dir) in (RemoteHPC.Running, RemoteHPC.Pending, RemoteHPC.Submitted)
        
        if already_submitted
            @info "Entity($(e.e)) was already submitted somehow..."
            curt = Dates.datetime2unix(now())
            if e ∉ m[TimingInfo]
                m[e] = TimingInfo(curt, curt, 0.0, 0.0, 0.0, 0.0, "", 0.0, 0.0)
            end
            set_state!(m, e, Running())
            continue
        end
        
        # So we don't have to always do abspath(server, remote_dir)
        jinfo[e].remote_dir = joinpath(server, m, e)

        for c in filter(x -> x.run, j.calculations)
            if eltype(c) == QE
                ex = load(server, Exec(server_info.pw_exec))
                ex.path = joinpath(dirname(ex.path), exec(c.exec))
                c.exec = ex
            elseif eltype(c) == Wannier90
                c.exec = load(server, Exec("wannier90"))
            end
        end

        j.environment = server_info.environment

        set_server_pseudos!(j, local_server(), m)
        local_save(j, jinfo[e].local_dir)
        j.dir = jinfo[e].remote_dir
    

        j.server = server_info.server
        set_server_pseudos!(j, server, m)

        suppress() do
            priority = e in m[NSCFSettings] || e in m[BaseCase] ? server_info.priority + 1 : server_info.priority
            submit(j; fillexecs = false, versioncheck=false, priority = priority)
        end
        curt = Dates.datetime2unix(now())
        if e ∉ m[TimingInfo]
            m[e] = TimingInfo(curt, curt, 0.0, 0.0, 0.0, 0.0, "", 0.0, 0.0)
        end
        set_state!(m, e, Submitted())
    end
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
    @error_capturing_threaded for e in @safe_entities_in(m, SimJob && TimingInfo && !Completed && !Submit && !Pulled)
        if !any(x -> x.run, e.job.calculations)
            continue
        end
        server = Server(e.job.server)
        curt   = Dates.datetime2unix(now())
        prevt  = e.cur_time
        dt     = curt - prevt
        s      = state(e.job)
        set_state!(m, e, s)
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
                abort(e.job)
                for c in e.job.calculations
                    if c.name == cur_running
                        break
                    end
                    c.run = false
                end
                log(e, "Aborted and resubmitted job during: $(Client.last_running_calculation(e.job).name).")
                should_rerun(m, e)
            end
        elseif s == RemoteHPC.Failed && ispath(server, e.remote_dir) &&
               filesize(server,
                        joinpath(e.remote_dir,
                                 Client.last_running_calculation(e.job).outfile)) == 0.0
            should_rerun(m, e)
        elseif s in (RemoteHPC.NodeFail, RemoteHPC.Cancelled)
            should_rerun(m, e)
        elseif isparseable(s)
            set_state!(m, e, Completed())
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
        s = Server(e.job.server)
        jdir = e.remote_dir
        if e.remote_dir != e.local_dir && ispath(s, e.remote_dir)
            @debugv 2 "Cleaning $jdir on Server $(s.name)"
            rm(s, jdir)
        end
        e[Done] = Done(true)
    end
end

struct OutputPuller <: System end
Overseer.requested_components(::OutputPuller) = (Done, SimJob, Completed, Pulled)
function Overseer.update(::OutputPuller, m::AbstractLedger)
    @error_capturing_threaded for e in @safe_entities_in(m, SimJob && Completed)
        server = Server(e.job.server)
           
        should_resubmit = false
        
        curt          = Dates.datetime2unix(now())
        running_calcs = filter(x -> x.run, e.job.calculations)
        pulled_files  = RemoteHPC.pull(e.job, e.local_dir, calcs=running_calcs)
        for c in running_calcs
            ofile = joinpath(e.local_dir, c.outfile)
            if !ispath(ofile) || filesize(ofile) == 0
                @debugv 2 "$(c.name) output file is missing or empty, resubmitting."
                c.run = true
                should_resubmit = true
                ispath(ofile) && rm(ofile)
            else
                c.run=false
            end
        end
        if should_resubmit
            should_rerun(m, e)
        else
            set_state!(m, e, Pulled())
        end
        m[e].postprocessing += Dates.datetime2unix(now()) - curt
    end
end

struct SimJobRemover <: System end

function Overseer.update(::SimJobRemover, m::AbstractLedger)
    for e in @safe_entities_in(m, Done && SimJob)
        if e.cleaned
            pop!(m[SimJob], e)
        end
    end
end

@component struct ShouldRerun
    data_to_pop::Set{DataType}
end
ShouldRerun(args::DataType...) = ShouldRerun(Set(args))
Base.push!(s::ShouldRerun, d::DataType) = push!(s.data_to_pop, d)

@component struct Rerun
    count::Int
end

struct Rerunner <: System end

Overseer.requested_components(::Rerunner) = (ShouldRerun, Rerun)

function Overseer.update(::Rerunner, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, ShouldRerun)
        from_scratch = false
        for d in e.data_to_pop
            trypop!(m[d], e)
            if d == SimJob
                trypop!(m[ServerInfo], e)
                trypop!(m[Unique], e)
                trypop!(m[Child], e)
                from_scratch = true
            end
        end
        trypop!(m[Done], e)
        if !from_scratch
            set_state!(m, e, Submit())
        end
        pop!(m[ShouldRerun], e) 
        
        m[e] = Rerun(e in m[Rerun] ? m[Rerun][e].count + 1 : 1)
    end
end

struct ErrorCorrector <: System end

function Overseer.update(::ErrorCorrector, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, Results && !ShouldRerun)
        if length(e.state.occupations) == 0
            @debugv 1 "ErrorCorrector: $(e.e) has an empty State in Results"
            # This results usually because of some server side issue where something crashed for no reason
            # Can also be because things converged before constraints were released (see process_Hubbard)
            if e in m[SimJob]
                if e.niterations == e.constraining_steps != 0
                    @debugv 1 "ErrorCorrector: $(e.e) had converged scf while constraints still applied, increasing Hubbard_conv_thr."
                    m[Template][e].calculation[:system][:Hubbard_conv_thr] = 1.5 * m[SimJob][e].job.calculations[1][:system][:Hubbard_conv_thr]
                    should_rerun(m, e, SimJob)
                else
                    set_flow!(m[SimJob][e].job, "" => true)
                end
            end
            should_rerun(m, e, BandsResults, Results, RelaxResults, FlatBands, Completed, Pulled)
        end
    end
    
    @error_capturing for e in @safe_entities_in(m, SimJob && !ShouldRerun && !Done && !Submit)
        if !any(x-> x.run, e.job.calculations) || e ∉ m[TimingInfo]
            @debugv 1 "ErrorCorrector: resubmitting $(e.e)"
            set_flow!(e.job, "" => true)
            should_rerun(m, e, BandsResults, Results, RelaxResults, FlatBands, Completed, Pulled)
        end
    end

    for e in @safe_entities_in(m, Done && Error)
        if occursin("HTTP.Exceptions.StatusError(500, \"POST\", \"/rm/?", string(e.err))
            pop!(m[Error], e)
        end
    end
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
    for e in @safe_entities_in(m, Results && Generation && !Parent)
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
        new_e = Entity(m, BaseCase(), Template(str, deepcopy(m[Template][base_e].calculation)), Generation(maximum(x->x.generation, m[Generation], init=0)))
        if base_e in m[RelaxSettings]
            m[RelaxSettings][new_e] = base_e
        end
        if base_e in m[HPSettings]
            m[HPSettings][new_e] = base_e
        end
        return false
    else
        return true
    end
end

function verify_groundstates!(m::AbstractLedger)
    gs_full = ground_state(filter(x->x.converged, @entities_in(m, Results)))
    gs_search = ground_state(filter(x->x.converged, @entities_in(m, Results && !Parent)))

    for gs in (gs_full, gs_search)
        if gs ∉ m[FlatBands]
            j = load(local_server(), Job(local_dir(m, gs)))
            o = outputdata(j, calcs = [j.calculations[1].name])[j.calculations[1].name]
            m[gs] = FlatBands(flatbands(o))
        end
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
    if isempty(n_unique) || (!stop && all(x -> x in m[Done] || x in m[Error], @entities_in(m, Trial && (RandomSearcher || Intersection))))
        set_mode!(m, :search)
        prepare(m)
        return
    end
        
    if m.ledger.stages[1].name == :core
        if length(@entities_in(m[SimJob] && !m[Error])) == 0 && all(x->x.cleaned, @entities_in(m, Done && !Error)) && all(x -> x in m[Done] || x in m[Error], @entities_in(m[Trial])) 
            # Check if stop condition is still met
            if stop && check_basecase(m)
                verify_groundstates!(m)
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

