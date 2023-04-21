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

            opath = joinpath(local_dir(m, e), "$(tc.name).out")
            if ispath(opath)
                rm(opath)
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

function process_structure_update!(m::Searcher, e::AbstractEntity, final_structure)
    if any(isnan, final_structure.cell)
            # We copy here otherwise we overwrite all relaxsettings
        e[RelaxSettings].symmetry = false
        log(e, "Found NaNs in final structure, rerunning with no symmetry")
        return typemax(Float64)
        
    else
        # Relatively arbitrary check for convergence
        orig_str = m[Template][e].structure
        orig_V   = Structures.volume(orig_str)
        new_V    = Structures.volume(final_structure)

        diff = Structures.ustrip(abs(orig_V - new_V))
        
        for (a1, a2) in zip(orig_str.atoms, final_structure.atoms)
            diff += Structures.ustrip(norm(a1.position_cart .- a2.position_cart))
        end
        log(e, "Diff of relaxed structure with original: $diff") 
        
        # Update structure for possibly rerunning
        new_template = deepcopy(m[Template][e])
        Structures.update_geometry!(new_template.structure, final_structure)
        
        if !isempty(m[BaseCase]) && Entity(oldest_parent(m, e)) == entity(m[BaseCase], 1)
            search_e = entity(m[StopCondition], 1)
            Structures.update_geometry!(m[Template][search_e].structure, final_structure)
        end
        
        return diff, new_template
    end
end

function Overseer.update(::RelaxProcessor, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, Pulled && SimJob && RelaxSettings)
        
        cname = e.variable_cell ? "vcrelax" : "relax"
        if !ispath(joinpath(e.local_dir, "$cname.out"))
            continue
        end
        
        curt = now()
        
        j = local_load(Job(e.local_dir))
        o = hubbard_outputdata(j; calcs = [cname])
        
        if isempty(o)
            continue
        end
        
        e.job[cname].run = false
        
        res = o[cname]

        if !any(x->x.name == "scf", j.calculations)
            results = results_from_output(res, oldest_parent(m, e) in m[BaseCase])
            
            m[e] = results
            
            if haskey(res, :bands) && haskey(res, :total_magnetization)
                bands = flatbands(res)
                m[e] = FlatBands(bands)
                
                fermi = res[:fermi]
                
                up   = res[:bands].up
                down = res[:bands].down

                n_up_conduction   = count(x->maximum(x.eigvals) > fermi, up)
                n_down_conduction = count(x->maximum(x.eigvals) > fermi, down)

                if n_up_conduction == 0 || n_down_conduction == 0
                    log(e, "There were no conduction bands for a spin channel, increasing nbnd and rerunning.")
                    
                    suppress() do
                        new_template = deepcopy(m[Template][e])
                        c = new_template.calculation
                        if haskey(c, :nbnd)
                            c[:nbnd] = ceil(Int, 1.2*c[:nbnd])
                        else
                            c[:nbnd] = ceil(Int, res[:n_KS_states] * 1.2)
                        end
                        should_rerun(m, e, new_template)
                    end
                end
            end
        end

        # We set the trial to the updated occupations in case we need to rerun
        if e in m[Trial] && !isempty(m[Results][e].state.occupations)
            m[e] = Trial(m[Results][e].state, PostProcess)
        end
        
        # Update structures
        if haskey(res, :total_force) && haskey(res, :final_structure)
            str = res[:final_structure]
            
            diff, new_template = process_structure_update!(m, e, str)
            
            relres = RelaxResults(res[:n_scf], res[:total_force][end], str, diff)
            m[e] = relres
            
            if diff > e.force_convergence_threshold && (!res[:converged] || !res[:finished])
                should_rerun(m, e, new_template)
                continue
                
            elseif res[:converged] && !res[:finished]
                # Sometimes something random goes wrong at the very end of the vcrelax calculation.
                # We run an scf with the final parameters in that case
                
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
                    for i = scf_id:length(e.job.calculations)
                        e.job.calculations[i].run = true
                    end
                        
                    log(e, "vcrelax went wrong even though it converged, could be symmetry issues.\nCreating an scf with the final structure for further steps.")
                    should_rerun(m, e)
                    continue
                end

            end
            m[e] = Done(false)
            
            
        elseif res[:finished] && !res[:converged]
            # Standard non converged vcrelax (at the first scf step)
            new_template = deepcopy(m[Template][e])
            calc = new_template.calculation
            
            if calc[:electron_maxstep] < 1000
                
                suppress() do 
                    calc[:electron_maxstep] *= 2
                end
                should_rerun(m, e, new_template)
                
                log(e, "vcrelax did not converge, increased electron maxstep to $(calc[:electron_maxstep])")
                
            elseif calc[:mixing_beta] > 0.01
                
                suppress() do
                    calc[:mixing_beta] *= 0.5
                end
                should_rerun(m, e, new_template)
                
                log(e, "vcrelax did not converge, lowered mixing beta to $(calc[:mixing_beta])")
                
            elseif e in m[Rerun]
                
                if e.job.calculations[end].name == cname
                    m[e] = Done(false)
                end
                
            else
                m[e] = Error(e, "vc-relax doesn't converge even after changing input parameters.")
            end
        else
            m[e] = Error(e, "Something went wrong with vcrelax calculation")
        end
        
        m[TimingInfo][e].postprocessing += Dates.datetime2unix(now()) - Dates.datetime2unix(curt)
    end
end
