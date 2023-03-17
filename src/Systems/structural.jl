"""
    Relaxor

Creates the inputs for a relaxation calculation from [`RelaxSettings`](@ref).
"""
struct Relaxor <: System end

Overseer.requested_components(::Relaxor) = (RelaxSettings, RelaxResults) 

function Overseer.update(::Relaxor, m::AbstractLedger)
    if isempty(m[RelaxSettings])
        return
    end
    
    @error_capturing for e in @safe_entities_in(m, SimJob && RelaxSettings)
        # We assume that a job has been created with an SCF calculation instead of a vcrelax,
        # so we highjack that one
        cname = e.variable_cell ? "vcrelax" : "relax"
        if any(x -> x.name == cname, e.job.calculations)
            continue
        end
        tc = e.job.calculations[end]
        suppress() do
            tc[:forc_conv_thr] = e.force_convergence_threshold
            tc[:etot_conv_thr] = e.energy_convergence_threshold
            tc[:ion_dynamics]  = e.ion_dynamics 
            tc[:cell_dynamics] = e.cell_dynamics
            if !e.symmetry
                tc[:nosym] = true
            end
            if ispath(joinpath(local_dir(m, e), "$(tc.name).out"))
                rm(joinpath(local_dir(m, e), "$(tc.name).out"))
            end
            if e.variable_cell 
                set_name!(tc, "vcrelax")
                tc[:calculation] = "vc-relax"
            else
                set_name!(tc, "relax")
                tc[:calculation] = "relax"
            end
        end
    end
end

"""
    RelaxProcessor

Processes the results of a relaxation and puts them in [`Results`](@ref) and
[`RelaxResults`](@ref).
"""
struct RelaxProcessor <: System end

Overseer.requested_components(::RelaxProcessor) = (Template, RelaxSettings, RelaxResults, Parent)

function Overseer.update(::RelaxProcessor, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, Pulled && SimJob && RelaxSettings)
        cname = e.variable_cell ? "vcrelax" : "relax"
        if !ispath(joinpath(e.local_dir, "$cname.out"))
            continue
        end
        curt = now()
        j = local_load(Job(e.local_dir))
        o = hubbard_outputdata(j; calcs = [cname])
        if !isempty(o)
            e.job[cname].run = false
            res = o[cname]
            results = results_from_output(res, e in m[BaseCase])
            m[e] = results
            if e in m[Trial] && !isempty(results.state.occupations)
                m[e] = Trial(results.state, PostProcess)
            end
            
            if haskey(res, :bands) && haskey(res, :total_magnetization)
                bands = flatbands(res)
                m[e] = FlatBands(bands)
            end
                
            if haskey(res, :total_force) && haskey(res, :final_structure)
                str = res[:final_structure]
                
                if any(isnan, str.cell)
                    # We copy here otherwise we overwrite all relaxsettings
                    t = deepcopy(m[RelaxSettings][e])
                    t.symmetry = false
                    m[e] = t
                    log(e, "Found NaNs in final structure, rerunning with no symmetry")
                    diff = typemax(Float64)
                else
                    # Relatively arbitrary check for convergence
                    orig_str = m[Template][e].structure
                    diff = abs(Structures.ustrip(Structures.volume(orig_str)) - Structures.ustrip(Structures.volume(str)))
                    
                    for (a1, a2) in zip(orig_str.atoms, str.atoms)
                        diff += Structures.ustrip(norm(a1.position_cart .- a2.position_cart))
                    end
                    log(e, "Diff of relaxed structure with original: $diff") 
                    
                    # Update structure for possibly rerunning
                    Structures.update_geometry!(m[Template][e].structure, str)
                    Structures.update_geometry!(e.job.structure, str)
                    if !isempty(m[BaseCase]) && e.e == entity(m[BaseCase], 1)
                        search_e = entity(m[Unique], 1)
                        Structures.update_geometry!(m[Template][search_e].structure,
                                                             res[:final_structure])
                    end
                end

                relres = RelaxResults(res[:n_scf], res[:total_force][end], str)
                m[e] = relres
                
                if diff > e.force_convergence_threshold && (!res[:converged] || !res[:finished])
                    should_rerun(m, e, SimJob)
                    continue
                elseif res[:converged] && !res[:finished]
                    scf_id = findfirst(x->x.name == "scf", e.job.calculations)
                    if scf_id === nothing
                        tc = deepcopy(e.job[cname])
                        suppress() do
                            tc[:calculation] = "scf"
                            Calculations.set_name!(tc, "scf")
                        end
                        cid = findfirst(x->x.name == cname, e.job.calculations)
                        insert!(e.job, cid + 1, tc)
                        scf_id = cid + 1
                    end
                    for i = scf_id:length(e.job.calculations)
                        e.job.calculations[i].run = true
                    end
                        
                    should_rerun(m, e) 
                end
                
                if e.job.calculations[end].name == cname
                    m[e] = Done(false)
                end
            else
                if res[:finished] && !res[:converged]
                    calc = m[Template][e].calculation
                    if calc[:electron_maxstep] < 1000
                        suppress() do 
                            calc[:electron_maxstep] *= 2
                        end
                        log(e, "vcrelax did not converge, increased electron maxstep to $(calc[:electron_maxstep])")
                        should_rerun(m, e, SimJob)
                    elseif calc[:mixing_beta] > 0.01
                        suppress() do
                            calc[:mixing_beta] *= 0.5
                        end
                        log(e, "vcrelax did not converge, lowered mixing beta to $(calc[:mixing_beta])")
                        should_rerun(m, e, SimJob)
                    elseif e in m[Rerun]
                        if e.job.calculations[end].name == cname
                            m[e] = Done(false)
                        end
                    else
                        m[e] = Error(e, "Something went wrong with vcrelax calculation")
                    end
                else
                    m[e] = Error(e, "Something went wrong with vcrelax calculation")
                end
            end
        end
        m[TimingInfo][e].postprocessing += Dates.datetime2unix(now()) - Dates.datetime2unix(curt)
    end
end
