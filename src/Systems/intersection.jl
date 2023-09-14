struct Intersector <: System end

function Overseer.requested_components(::Intersector)
    return (Trial, Generation, Results, Unique, IntersectionSearcher, Intersection)
end

function Overseer.prepare(i::Intersector, m::AbstractLedger)
    for c in Overseer.requested_components(i)
        Overseer.ensure_component!(m, c)
    end
end

function Overseer.update(::Intersector, m::AbstractLedger)
    lck = ReentrantLock()
    gencomp = m[Generation]
    unique_es = collect(@entities_in(m, Unique && Results))
    intersections_per_generation = Int[]
    tgen = 0
    for g in @entities_in(m[Generation])
        if g.generation == 0
            continue
        end
        curlen = length(intersections_per_generation)
        if g.generation > curlen
            if g.generation - 1 > curlen
                append!(intersections_per_generation, zeros(Int, g.generation - 1 - curlen))
            end
            if !(g in m[Unique]) && !(g in m[Done]) && !(g in m[Error])
                push!(intersections_per_generation, 1)
            end
        else
            if !(g in m[Unique]) && !(g in m[Done]) && !(g in m[Error])
                intersections_per_generation[g.generation] += 1
            end
        end
    end

    # We determine the minimum distance for each generation by looking at the mean distance between
    # uniques in that generation
    pairwise_distances = [Euclidean()(e1.state, e2.state) for e1 in unique_es, e2 in unique_es]
    generation_distances = map(1:length(intersections_per_generation)) do gen
        n = 0
        totdist = 0.0
        gen_es = Iterators.filter(x->m[Generation][x].generation <= gen, unique_es)
        @inbounds for (i, e) in enumerate(gen_es)
            for (j, e1) in enumerate(gen_es)
                if e1.e != e.e
                    n += 1
                    totdist += pairwise_distances[i, j]
                end
            end
        end
        return n < 2 ? typemax(Float64) : totdist/n
    end

    for (i, intpergen) in enumerate(intersections_per_generation)
        if intpergen == 0
            for e in unique_es
                if m[Generation][e] == i
                    pop!(m[IntersectionSearcher], e)
                end
            end
        end
    end


    search_e = entity(m[IntersectionSearcher], 1)
    potential_intersections = Dict{Int, Vector{Tuple{Float64, Trial, Intersection}}}()
    for (p, es) in pools(m[IntersectionSearcher])
        tot = 0
        # Iterate over all entities that haven't been added to the intersection searcher
        for e1 in Iterators.filter(x -> !(x in es), unique_es)
            g1 = gencomp[e1].generation
            intersection_gen = g1 + 1
            if g1 > length(generation_distances)
                continue
            end
            mindist = p.mindist * generation_distances[g1]
            if intersection_gen > length(intersections_per_generation)
                push!(intersections_per_generation, 0)
            end
            @inbounds if intersections_per_generation[intersection_gen] < p.max_intersections_per_generation
                
                if intersection_gen ∉ keys(potential_intersections)
                    potential_intersections[intersection_gen] = Tuple{Float64, Trial, Intersection}[]
                end
                intersections = potential_intersections[intersection_gen]
                if length(intersections) < 2 * p.max_intersections_per_generation
                    
                    others = filter(unique_es) do x
                        g2 = gencomp[x].generation
                        return x.e != e1.e && g2 <= g1
                    end
                        
                    ilock = ReentrantLock() 
                    Threads.@threads for e2 in others
                        tstate = 0.5*(e1.state + e2.state)

                        # Find closest prev unique
                        dist, minid = findmin(x -> Euclidean()(x.state, tstate), unique_es)
                        sink = unique_es[minid]
                        trial_to_sink = Euclidean()(sink.state, tstate)
                        if trial_to_sink > mindist
                            # from this unique, find all prev trials that led to it
                            funnel = filter(x -> x in m[Trial] && sink.e != x.e && !isempty(x.state.occupations) && Euclidean()(x.state, sink.state) < 1e-2, @entities_in(m, Results))

                            # We test that all prev trials in funnel are further than mindist
                            # and that the distance to all prev trials is larger than towards the sink
                            # Very crude only distance based determination of funnel
                            test = all(funnel) do t
                                fstate = m[Trial][t].state
                                funnel_to_sink = Euclidean()(sink.state, fstate)
                                trial_to_funnel = Euclidean()(tstate, fstate)
                                return trial_to_funnel > mindist && trial_to_funnel > funnel_to_sink
                            end
                        
                            if test
                                lock(ilock)
                                try
                                    push!(intersections, (dist, Trial(tstate, IntersectionMixed), Intersection(e1.e, e2.e)))
                                finally
                                    unlock(ilock)
                                end
                            end
                        end
                    end
                end
            end
            m[IntersectionSearcher][e1] = es[1]
        end

        @inbounds for (gen, intersections) in potential_intersections
            n_int = length(intersections)
            # First recalculate mindist with respect to the new additions
            if n_int > 1
                Threads.@threads for i in 1:n_int
                    int = intersections[i]
                    min_to_others = minimum(j->Euclidean()(intersections[j][2].state, int[2].state), 1:i-1, init=int[1])
                    min_to_prev_trials = minimum(x->Euclidean()(x.state, int[2].state), m[Trial])
                    intersections[i] = (min(int[1], min(min_to_others, min_to_prev_trials)), int[2], int[3])
                end
            end
            # Then insert n_new best intersections
            n_new = min(n_int, p.max_intersections_per_generation - intersections_per_generation[gen])
            new_intersections = partialsort!(intersections, 1:n_new, rev=true, by = x -> x[1])
            n_added = 0
            for (dist, trial, intersection) in new_intersections
                if dist > p.mindist * generation_distances[gen - 1]
                    add_search_entity!(m, search_e, Generation(gen), trial, intersection)
                    n_added += 1
                end
            end
            if n_added != 0
                @debugv 2 "Added $n_added intersections at Generation($gen) with smallest distance $(new_intersections[n_added][1])"
            end
            tot += n_added
        end
        if tot != 0
            @debugv 1 "Total intersections added for pool $(pool(m[IntersectionSearcher], es[1])): $tot"
        end
    end
end

struct RandomTrialGenerator <: System end
function Overseer.requested_components(::RandomTrialGenerator)
    return (RandomSearcher, Intersection, BaseCase)
end

function rand_trial(n_orb_per_at, n_elec_per_at)
    
    occs = map(zip(n_orb_per_at, n_elec_per_at)) do (norb, nelec)
        nangles = div(norb * (norb - 1), 2)
        
        rand_angles() = Angles([π * (rand() - 0.5) for i in 1:nangles], 1.0, (norb, norb))
        
        ox_state_offset = rand([-1, 0, 1])
        diagvec = zeros(2norb)
        
        while sum(diagvec) < min(nelec + ox_state_offset, 2norb)
            diagvec[rand(1:2norb)] = 1.0
        end
        
        if norb == 1
            return ColinMatrixType(diagm(0 => diagvec[1:1]), diagm(0 => diagvec[2:2]))
        else
            D = DFWannier.MagneticVector(diagvec)
            V = DFWannier.ColinMatrix(Matrix(rand_angles()), Matrix(rand_angles()))
            return Matrix(Eigen(D, V))
        end
    end

    return Trial(State(occs), RandomMixed)
end

function Overseer.update(::RandomTrialGenerator, m::AbstractLedger)
    if isempty(m[RandomSearcher]) || isempty(m[BaseCase])
        return
    end
    # First make sure the base case calculation is finished
    base_e = entity(m[BaseCase], 1)
    if !all_children_done(m, base_e) || youngest_child(m, base_e) ∉ m[Results] || youngest_child(m, base_e) in m[Error]
        return
    end
    
    base_state = m[Results][base_e].state
    if isempty(base_state.occupations)
        @error "Something went wrong with the basecase calculation"
        return
    end
    rand_search_comp = m[RandomSearcher]
    rand_search_e    = entity(rand_search_comp, 1)
    rand_search      = rand_search_comp[rand_search_e]

    # Wait until all intersections based on random generation have been finished
    random_search_entities = @entities_in(m, RandomSearcher && Trial)
    
    all_in_results = all(x -> x ∈ m[Results], random_search_entities)
    all_intersection_finished = all(x -> x ∈ m[Done] || x ∈ m[Error], @entities_in(m, Intersection))
    if !all_in_results || !all_intersection_finished 
        return
    end

    # Set Generation to the maximum one considering the intersections
    maxgen = maximum_generation(m)
    
    template_str = m[Template][rand_search_e].structure
    nelec    = round.(Int, base_state.totoccs)
    norb     = size.(base_state.occupations, 1)
    
    nsearchers = rand_search.nsearchers
    for i = 1:nsearchers
        
        e = add_search_entity!(m, rand_search_e,
                               rand_trial(norb, nelec),
                               Generation(maxgen + 1),
                               RandomSearcher(nsearchers))
                               
        if Hybrid in m && length(m[Hybrid]) != 0
            m[e] = Hybrid()
        end
        
    end
    @debug "New random search at Generation($(maxgen + 1))."
end
