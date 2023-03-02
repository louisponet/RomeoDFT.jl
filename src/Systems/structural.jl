struct Relaxor <: System end

Overseer.requested_components(::Relaxor) = (Template, RelaxSettings, RelaxResults, Parent, RelaxChild)

function Overseer.update(::Relaxor, m::AbstractLedger)
    if isempty(m[RelaxSettings])
        return
    end

    relset = entity(m[RelaxSettings], 1)
   
    for e in @entities_in(m, Unique && Results && !RelaxChild)
        enew = Entity(m, Parent(e.e), m[Generation][e], Trial(e.state, PostProcess))
        m[Template][enew] = e
        m[RelaxSettings][enew] = relset
        m[e] = RelaxChild(enew)
    end
end

struct RelaxProcessor <: System end

Overseer.requested_components(::RelaxProcessor) = (Template, RelaxSettings, RelaxResults, Parent, RelaxChild)

function Overseer.update(::RelaxProcessor, m::AbstractLedger)
    if isempty(m[RelaxSettings])
        return
    end
    lck = ReentrantLock()
    function lock_(f)
        lock(lck)
        try
            f()
        finally
            unlock(lck)
        end
    end
    to_pop = Entity[]
    @sync for e in @entities_in(m, SimJob && TimingInfo && RelaxSettings && !RelaxResults)
        cname = e.variable_cell ? "vcrelax" : "relax"
        if ispath(joinpath(e.local_dir, "$cname.out"))
            # Threads.@spawn begin
                curt = now()
                try
                    j = local_load(Job(e.local_dir))
                    o = hubbard_outputdata(j; calcs = [cname])
                    if !isempty(o)
                        e.job[cname].run = false
                        res = o[cname]
                        results = results_from_output(res)
                        if haskey(res, :total_force)
                            if res[:converged]
                                relres = RelaxResults(res[:n_scf], res[:total_force][end], res[:final_structure])
                                
                                lock_() do
                                    m[Results][e] = results
                                    m[RelaxResults][e] = relres
                                    m[Done][e] = Done(false)
                                end

                                if haskey(res, :bands) && haskey(res, :total_magnetization)
                                    bands = flatbands(res)
                                    lock_() do
                                        m[FlatBands][e] = FlatBands(bands)
                                    end
                                end
                            else
                                str = res[:final_structure]
                                if any(isnan, str.cell)
                                    m[RelaxSettings].data[1].symmetry = false
                                    lock_() do
                                        push!(to_pop, e.e)
                                    end
                                else
                                    orig_str = deepcopy(m[Template][e].structure)
                                    Structures.update_geometry!(orig_str, str)
                                    @debug "Updating structure of $(e.e) and rerunning" 
                                    lock_() do
                                        m[Template][e] = Template(orig_str, m[Template][e].calculation)
                                        push!(to_pop, e.e)
                                    end
                                end
                            end
                        else
                            relres = RelaxResults(0, Inf, e.job.structure)
                            lock_() do
                                m[Results][e] = results
                                m[RelaxResults][e] = relres
                                m[Done][e] = Done(false)
                            end
                        end
                    end
                catch err
                    lock_() do 
                        @debug "Error for $(e.e) $(stacktrace(catch_backtrace()))" 
                        m[e] = Error(err, stacktrace(catch_backtrace()))
                    end
                end
                e.postprocessing += Dates.datetime2unix(now()) - Dates.datetime2unix(curt)
            # end
        end
    end
    for e in to_pop
        pop!(m[SimJob], e)
        pop!(m[TimingInfo], e)
    end
end
