function rand_trial(s::Structure, norb, nelec)
    mag_ats = filter(x -> x.dftu.U != 0.0 || sum(x.magnetization) != 0.0, s.atoms)
    nat = length(mag_ats)
    occs = map(1:nat) do ia
        n = norb[ia]
        nangles = div(n * (n - 1), 2)
        rand_angles() = Angles([π * (rand() - 0.5) for i in 1:nangles], 1.0, (n, n))
        ox_state_offset = rand([-1, 0, 1])
        diagvec = zeros(2n)
        while sum(diagvec) < min(nelec + ox_state_offset, 2n)
            diagvec[rand(1:2n)] = 1.0
        end
        if n == 1
            return DFWannier.ColinMatrix(diagm(0 => diagvec[1:1]), diagm(0 => diagvec[2:2]))
        else
            diag_occupations = DFWannier.MagneticVector(diagvec)
            rotmat = DFWannier.ColinMatrix(Matrix(rand_angles()), Matrix(rand_angles()))
            return Matrix(Eigen(diag_occupations, rotmat))
        end
    end

    return Trial(State(occs), RandomMixed)
end

function mix_eulerangles(flies, trial, generation, norb, nelec)
    new_occs = Vector{Trial}(undef, length(flies))
    handled = falses(length(flies))
    rand_trials = 0
    for (i1, e1) in enumerate(flies)
        if handled[i1]
            continue
        end

        new_angles = deepcopy(e1[Results].state.angles)
        for (i2, e2) in enumerate(flies)
            if i2 == i1 || (handled[i2] && new_occs[i2].origin == RandomMixed)
                continue
            end
            if !e2.converged
                @debugv 2 "$(e2.e) did not converge"
                new_occs[i2] = rand_trial(e1.template_structure, norb, nelec)
                rand_trials += 1
                handled[i2] = true
                continue
            end

            s1 = e1[Results].state
            s2 = e2[Results].state
            d = Euclidean()(s1, s2)
            if d < 1e-2 
                new_occs[i2] = rand_trial(e1.template_structure, norb, nelec)
                rand_trials += 1
                handled[i2] = true
            else
                factor = dft_energy(e2) * e2.β * exp(-e2.γ * d^2)

                for i in 1:length(new_angles)
                    a1 = s1.angles[i]
                    a2 = s2.angles[i]
                    for j in 1:length(a1.θs)
                        angle1, angle2 = a1.θs[j], a2.θs[j]
                        new_angles[i].θs[j] += rem2pi(factor * (angle2 - angle1) +
                                                      e1.α * 2π * rand(), RoundNearest)
                    end
                end
            end
        end
        occ = map(1:length(new_angles)) do i
            rotmat = Matrix(new_angles[i])
            return Matrix(Eigen(e1[Results].state.eigvals[i], rotmat))
        end
        new_occs[i1] = Trial(State(occ), EulerAngleMixed)
        handled[i1] = true
    end
    @debugv 2 "$rand_trials random trials..."
    return new_occs
end

function mix_linear(flies, trial, generation, norb, nelec)
    new_occs = Vector{Trial}(undef, length(flies))
    handled = falses(length(flies))
    rand_trials = 0
    
    # to settle attraction factor
    min_energy = minimum(x->x.total_energy, filter(x->x.converged, flies))
    max_energy = maximum(x->x.total_energy, filter(x->x.converged, flies))
    attraction_factor(energy, β) = (energy - max_energy)/(min_energy - max_energy) * β 

    for (i1, e1) in enumerate(flies)
        if handled[i1]
            continue
        end

        new_occ = deepcopy(e1[Results].state.occupations)
        
        for (i2, e2) in enumerate(flies)
            if i2 == i1 || (handled[i2] && new_occs[i2].origin == RandomMixed)
                continue
            end
            if !e2.converged
                @debugv 2 "$(e2.e) did not converge"
                new_occs[i2] = rand_trial(e1.template_structure, norb, nelec)
                rand_trials += 1
                handled[i2] = true
                continue
            end

            s1 = e1[Results].state
            s2 = e2[Results].state
            d = Euclidean()(s1, s2)
            if d < 1e-2
                new_occs[i2] = rand_trial(e1.template_structure, norb, nelec)
                rand_trials += 1
                handled[i2] = true
                continue
            end
            
            factor = attraction_factor(e2.total_energy, e2.β) * exp(-e2.γ * d^2)
            @assert 0 <= factor <= 1 "Attraction factor error $factor"
            @debugv 2 "Attracting fly $i1 to fly $i2 with strength $factor"
            new_occ .+= factor .* (s2.occupations .- s1.occupations)
        end

        #Add some randomness without changing the trace
        for (io, o) in enumerate(new_occ)
            totocc = tr(o)
            vals, vecs = eigen(o)
            ts = 0.0
            for i = 1:length(vals) - 1
                vals[i] += e1.α * randn()/2
                ts += vals[i]
            end
            vals[end] = totocc - ts
            new_occ[io] = Matrix(Eigen(vals, vecs))
        end

        new_occs[i1] = Trial(State(new_occ), LinearMixed)
        handled[i1] = true
    end
    @debugv 2 "$rand_trials random trials..."
    return new_occs
end

function mix_random(flies, trial, generation, norb, nelec)
    e1 = flies[1]
    nflies = length(flies)
    return [rand_trial(e1.template_structure, norb, nelec) for i = 1:nflies]
end

"""
    FireFly

This controls the FireFly [`Simulation`](@ref). It gathers the results from the [`Trials`](@ref Trial) of the previous [`Generation`](@ref) and
creates new [`Trials`](@ref Trial) for the next [`Generation`](@ref).
This uses a variation on the update rule
```math
n_i^{t+1} = n_i^{t} + \\sum_j \\beta e^{-\\gamma r_{ij}^2} (n_i^{t} - n_j^{t}) + \\alpha \\varepsilon
```
"""
struct FireFly <: System end
function Overseer.requested_components(::FireFly)
    return (Simulation, Results, Generation, Trial, SimJob, TimingInfo, ServerInfo, BaseCase)
end

function Overseer.update(::FireFly, m::AbstractLedger)
    sim_comp = m[Simulation]
    res_comp = m[Results]

    # First make sure the base case calculation is finished
    if isempty(m[BaseCase]) || !(entity(m[BaseCase], 1) in m[Results])
        return
    end

    # Check if all intersections are done
    if any(x-> !(x in m[Results] || x in m[Error]), @entities_in(m, Intersection))
        return
    end

    base_state = m[Results][entity(m[BaseCase],1)].state
    if isempty(base_state.occupations)
        error("Something went wrong with the basecase calculation")
    end

    nelec = round(Int, sum(base_state.totoccs) / length(base_state.totoccs))
    norb = size.(base_state.occupations, 1)
    
    for e in @entities_in(sim_comp && !m[Generation])
        # Random Initialization
        m[e] = rand_trial(e.template_structure, norb, nelec)
        m[e] = Generation(e.current_generation)
    end
    trial_comp = m[Trial]
    gen_comp = m[Generation]
    simn = simname(m)

    for (sim, es) in pools(sim_comp)
        if !all(x -> x in res_comp, es)
            continue
        end
        curgen = maximum(x->gen_comp[x].generation, es)

        flies = map(y -> m[y], filter(e -> gen_comp[e].generation == curgen, es))

        issue = false
        
        # Check on whether each fly has actually an output that's usable
        for e in es
            if isempty(m[Results][e].state.occupations)
                pop!(m[Results], e)
                set_flow!(m[SimJob][e].job, "" => true)
                m[e] = Submit()
                issue = true
            elseif e in m[Error]
                pop!(m[Error], e)
                issue = true
            end
        end
        issue && continue
        prevtrials = @entities_in(trial_comp && !sim_comp)
        if !isempty(prevtrials)
            dists = zeros(length(flies))
            Threads.@threads for i in 1:length(flies)
                f = flies[i]
                if f[Results].converged
                    dists[i] = minimum(y -> Euclidean()(y.state, f[Results].state), prevtrials)
                end
            end
        else
            dists = ones(length(flies))
        end
        conv_flies = filter(x->x.converged, flies)

        finished = !any(x -> x > 1e-2, dists) # If no new states have been found

        last_postprocess_gen = isempty(prevtrials) ? 0 : maximum(x->x.generation, prevtrials)
        finished = finished || (isempty(conv_flies) && curgen - last_postprocess_gen > 2)  # If no converged flies for 2 generations
                
        if finished
            @debug "$simn - No new states at generation ($(sim.current_generation + 1))."
            for e in flies
                m[e] = Archived("", false)
            end
            if sim.current_best != Entity(0)
                m[sim.current_best] = Archived("", false)
            end
            continue
        elseif !isempty(prevtrials)
            @debug "$simn - Found $(length(filter(x->x > 1e-2, dists))) states at Generation $curgen."
        end

        # Potentially update the best entity
        conv_flies = filter(x->x.converged, flies)
        if isempty(conv_flies)
            @debug "None of the flies converged, trying a new batch of randomized flies."
            new_occ = map(e -> rand_trial(e.template_structure, norb, nelec), flies)
        else
            mine, id = findmin(x -> res_comp[x].total_energy, conv_flies)
            if sim.current_best == Entity(0) || !(sim.current_best in res_comp) || mine < res_comp[sim.current_best].total_energy
                sim.current_best = conv_flies[id]
                @debug "New global minimum found by $(conv_flies[id].e) at Generation($curgen): $mine Ry."
                if BaseCase in m
                    base_e = entity(m[BaseCase], 1)
                    str = m[Template][base_e].structure
                    magats = filter(x -> sum(x.magnetization) != 0 && x.dftu.U != 0, str.atoms)
                    curmags = map(x -> x.magnetization[3], magats)
                    magnetizations = map(x -> abs(x) < 1e-2 ? 1e-5 : sign(x), conv_flies[id][Results].state.magmoms)[1:length(magats)]
                    
                    if any(!iszero, curmags .- magnetizations) && any(!iszero, curmags .+ magnetizations)
                        
                        @debug "Updating \"vanilla\" case to magnetizations: $(join(string.(magnetizations), "; "))"
                        
                        base_e in res_comp && pop!(res_comp, base_e)
                        if base_e in m[SimJob]
                            sj = m[SimJob][base_e]
                            if state(sj.job) in (RemoteHPC.Submitted, RemoteHPC.Pending, RemoteHPC.Running)
                                abort(sj.job)
                            end
                            pop!(m[SimJob], base_e)
                        end

                        for (mag, at) in zip(magnetizations, magats)
                            at.magnetization = [0, 0, mag]
                        end
                    end
                end
            end
            
            if sim.mixing_mode ∈ (EulerAngleMixing, UnknownMixing)
                new_occ = mix_eulerangles(flies, m[Trial], gen_comp, norb, nelec)
            elseif sim.mixing_mode == LinearMixing
                new_occ = mix_linear(flies, m[Trial], gen_comp, norb, nelec)
            elseif sim.mixing_mode == RandomMixing
                new_occ = mix_random(flies, m[Trial], gen_comp, norb, nelec)
            else
                error("Sim mode not found, this should never happen.")
            end
            
            dists = zeros(length(flies))
            Threads.@threads for i in 1:length(flies)
                f = new_occ[i]
                dists[i] = Euclidean()(flies[i][Results].state, f.state)
            end
            @debugv 2 "$simn - Distances to prev result: $(join(dists, " - "))"
        end
        
        
        for (t, f) in zip(new_occ, flies)
            e = Entity(m, Generation(sim.current_generation + 1), t)
            sim_comp[e] = f
            if Hybrid in m && f in m[Hybrid]
                m[Hybrid][e] = Hybrid()
            end
        end
        @debug "$simn - New generation ($(sim.current_generation + 1))"
        sim.current_generation += 1
        sim.α *= 0.8
    end
end

"""
    PostFireflyExplorer

Takes states newly found by the [`FireFly`](@ref) simulation that are unique and generates the [`SCFSettings`](@ref), [`BandsSettings`](@ref),
[`NSCFSettings`](@ref) and [`ProjwfcSettings`](@ref) for the post processing workflow.
After this, the flies themselves are [`Archived`](@ref).
"""
struct PostFireflyExplorer <: System end
Overseer.requested_components(::PostFireflyExplorer) = (Archived,)

function Overseer.update(::PostFireflyExplorer, m::AbstractLedger)
    results = m[Results]
    generation = m[Generation]
    simulation = m[Simulation]
    nsettings = m[NSCFSettings]
    curgen = isempty(generation) ? 0 : (isempty(simulation) ? 0 : maximum(x->x.current_generation, simulation))
    new_states = 0
    for e in @entities_in(results && generation && simulation && !m[Archived])
        if e.generation < curgen
            if e.current_best == Entity(e)
                continue
            end
            m[e] = Archived("", false)
        end
        if e.converged && !(e in nsettings)

            # test_states = (e.state, State([ColinMatrix(x[Down()], x[Up()]) for x in e.state.occupations]))
            test_states = (e.state,)
            for s in test_states
                exists = false
                for e2 in @entities_in(m[Trial] && nsettings)
                    s2 = e2.state
                    d = Euclidean()(s, s2)
                    if d < 1e-2
                        exists = true
                    end
                end
                if !exists
                    occs    = generate_Hubbard_occupations(s, e.template_structure)
                    kpts    = (e.template_calculation.data[1].data...,)
                    scf     = SCFSettings(Dict(:system => Dict(:Hubbard_occupations => occs,
                    :smearing => "gaussian"),
                    :control => Dict(:disk_io => "low")),
                    kpts)
                    nscf    = NSCFSettings(; kpoints = (kpts[1:3]...,))
                    bands   = BandsSettings()
                    projwfc = ProjwfcSettings(; Emin = e.fermi - 20, Emax = e.fermi + 10)
                    # mag_els = unique(map(x->x.element.symbol, filter(a-> norm(a.magnetization) != 0 || a.dftu.U != 0.0, e.template_structure.atoms)))
                    # projs = Dict()
                    # for el in mag_els
                    #     if el in (:U, :Gd)
                    #         projs[el] = [Projection("f")]
                    #     elseif el in (:N, )
                    #         projs[el] = [Projection("p")]
                    #     elseif el in (:K,)
                    #         projs[el] = [Projection("s")]
                    #     else
                    #         projs[el] = [Projection("d")]
                    #     end
                    # end
                    # norm_els = unique(map(x->x.element.symbol, filter(a-> norm(a.magnetization) == 0 || a.dftu.U == 0.0, e.template_structure.atoms)))
                    # for el in norm_els
                    #     projs[el] = [Projection(el in (:O, :Sb) ? "p" : "s")]
                    # end
                    # wan = WannierSettings(projections=projs, plot_wannier=true)
                    # Entity(m, Template(e.template_structure, e.template_calculation), e[Generation], scf, nscf, bands, projwfc, wan, Trial(s))
                    new_e = Entity(m, Template(e.template_structure, e.template_calculation),
                            e[Generation], scf, bands, nscf, projwfc, Trial(s, m[Trial][e].origin))
                    if Hybrid in m && e in m[Hybrid]
                        m[new_e] = Hybrid()
                    end
                    new_states += 1
                end
            end
        end
    end
    simn = simname(m)
    if new_states != 0
        @debugv 2 "$simn - Found $new_states new states to explore."
    end
end

"""
    Archiver

Takes care of archiving certain entities that have finished their job, it saves all their components to a `archive.jld2` file inside the
`local_dir` of their [`SimJob`](@ref).
"""
struct Archiver <: System end

function Overseer.update(::Archiver, m::AbstractLedger)
    for e in @entities_in(m[Archived])
        if !e.isarchived
            p = joinpath(m.rootdir, "job_backups", "$(e.e.id)")
            mkpath(p)
            archive_path = joinpath(p, "archive.jld2")
            comps_to_store = []
            for (n, c) in components(m)
                if n == Archived
                    continue
                end
                if e in c
                    push!(comps_to_store, pop!(c, e))
                end
            end
            JLD2.save(archive_path, "components", comps_to_store, "version",
                      DATABASE_VERSION)
            e.archive_path = archive_path
            e.isarchived = true
        end
    end
end

