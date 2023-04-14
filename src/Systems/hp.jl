"""
    HPCreator

Creates the input for a HP calculation from [`HPSettings`](@ref).
"""
struct HPCreator <: System end

Overseer.requested_components(::HPCreator) = (HPSettings, HPResults)

function Overseer.update(::HPCreator, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, SimJob && HPSettings)
        if any(x->x.name == "hp", e.job.calculations)
            continue
        end
        # Inject HP into all entities before they are submitted
        hp_calc = Calculation{QE}(name  = "hp", exec=Exec(path="hp.x"),
                                  flags = Dict(:inputhp => Dict{Symbol, Any}(:nq1 => e.nq[1],
                                                                             :nq2 => e.nq[2],
                                                                             :nq3 => e.nq[3],
                                                                             :conv_thr_chi => e.conv_thr_chi,
                                                                             :find_atpert  => e.find_atpert)))
        push!(e.job, hp_calc)
    end
end

function update_hubbard_u!(str, hub_ats)
    ihub = 1
    totdiff = 0.0
    while ihub <= length(hub_ats)
        atsym = hub_ats[ihub][1]
        for a in str[element(atsym)]
            hub = hub_ats[ihub]
            a.name = hub[2]
            totdiff += abs(a.dftu.U - hub[3])
            a.dftu.U = hub[3]
            ihub += 1
        end
    end
    return totdiff/length(hub_ats)
end

"""
    HPProcessor

Processes the outputs of a HP calculation and puts them in [`HPResults`](@ref).
"""
struct HPProcessor <: System end

function getfirst_in_outputs(output, name)
    for (c, outdata) in output
        if haskey(outdata, name)
            return outdata[name]
        end
    end
end

function setup_scf_for_hp!(m, e, o, insulating_from_hp=false)
    if any(x->x.name == "scf_for_U", e.job.calculations)
        m[e] = Error(e, "Already ran with scf_for_U and still errors.")
        return false
    end
    bands = getfirst_in_outputs(o, :bands)
    fermi = getfirst_in_outputs(o, :fermi)

    scf_calc = deepcopy(m[Template][e].calculation)
    
    suppress() do
        set_name!(scf_calc, "scf_for_U")
        scf_calc.run = true
        delete_Hubbard!(scf_calc)
    end
    
    insulating = insulating_from_hp || (bands !== nothing && fermi !== nothing && DFControl.bandgap(bands, fermi) > 0.1)
    if insulating
        
        totmags = getfirst_in_outputs(o, :total_magnetization)
        totmags === nothing && return false
        totmag = totmags[end]
        
        n_ks = getfirst_in_outputs(o, :n_KS_states)
        n_ks === nothing && return false
        
        suppress() do 
            scf_calc[:tot_magnetization] = round(Int, totmag)
            scf_calc[:occupations] = "fixed"
            scf_calc[:nbnd] = n_ks
            scf_calc[:startingpot] = "file"
            scf_calc[:startingwfc] = "file"
            delete!(scf_calc, :degauss)
            delete!(scf_calc, :smearing)
        end
        
        for a in e.job.structure.atoms
            a.magnetization = [0,0,0]
        end
    end
    
    insert!(e.job.calculations, length(e.job.calculations), scf_calc)
    e.job["hp"].run = true
    
    return true
end

function Overseer.update(::HPProcessor, m::AbstractLedger)
    @error_capturing for e in @safe_entities_in(m, Pulled && SimJob && HPSettings && !ShouldRerun)
        if !ispath(joinpath(e.local_dir, "hp.out"))
            continue
        end
        j = local_load(Job(e.local_dir))
        o = outputdata(j)
        if isempty(o)
            continue
        end
        if get(o["hp"], :error, false)
            
            if o["hp"][:fermi_dos] < 0.1
                # We assume then it's insulating
                if setup_scf_for_hp!(m, e, o, true)
                    log(e, "HP: Fermi level shift 0. Creating insulating, 2 step HP job")
                    # Nothing to pop since we just want to rerun starting from the new scf
                    should_rerun(m, e)
                else
                    m[e] = Error(e, "Couldn't generate insulating scf for HP")
                end
                continue
            end
                
            # Probably too low cutoffs
            if m[Template][e].calculation[:ecutwfc] >= 100
                m[e] = Error(e, "Cutoff larger than 100 and still no HP results")
                continue
            end
            
            calc = m[Template][e].calculation
            tcut = calc[:ecutwfc] *= 1.2
            log(e, "HP: Fermi level shift too big. Increasing cutoff to $tcut")
            
            if haskey(calc, :ecutrho)
                calc[:ecutrho] *= 1.2
            end
            
            should_rerun(m, e, SimJob)
            
        elseif !haskey(o["hp"], :Hubbard_U)
            
            if setup_scf_for_hp!(m, e, o)
                log(e, "HP failed, creating a scf to run before it and try again")
                should_rerun(m, e)
                
            elseif e ∉ m[Rerun]
                m[e] = Error(e, "Issue with Hubbard output")
                if e.job.calculations[end].name == "hp"
                    m[e] = Done(false)
                end
            end
            
            continue
        else
            hub_ats = o["hp"][:Hubbard_U]
            m[e] = HPResults(hub_ats)
            
            diff = update_hubbard_u!(m[Template][e.e].structure, hub_ats)
            if !isempty(m[BaseCase]) && Entity(oldest_parent(m, e)) == entity(m[BaseCase], 1)
                update_hubbard_u!(m[Template][entity(m[StopCondition],1)].structure, hub_ats)
            end
            
            if diff > e.U_conv_thr
                log(e, "HP U diff with original structure: $diff, rerunning with updated U.")
                
                files = readdir(joinpath(e.local_dir))
                for f in files
                    if occursin("Hubbard_parameters", f)
                        rm(joinpath(e.local_dir, f))
                    end
                end
                should_rerun(m, e, SimJob)
                
            else
                log(e, "HP U diff with original structure: $diff, thr = $(e.U_conv_thr) --> Converged.")
                if e.job.calculations[end].name == "hp"
                    m[e] = Done(false)
                end
            end
        end
    end
end
