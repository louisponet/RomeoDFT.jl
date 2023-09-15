struct Intersector <: System end

function Overseer.requested_components(::Intersector)
    return (Trial, Generation, Results, Unique, IntersectionSearcher, Intersection)
end

function Overseer.prepare(i::Intersector, m::AbstractLedger)
    for c in Overseer.requested_components(i)
        Overseer.ensure_component!(m, c)
    end
end

# TODO: Think about whether we need all the max number of intersection per generation etc
function Overseer.update(::Intersector, m::AbstractLedger)
    lck = ReentrantLock()
    gencomp = m[Generation]
    
    unique_es = collect(@entities_in(m, Unique && Results))

    # The reason for this weirdness is that like this we only have to iterate over
    # m[Generation] once (calling intersections_per_generation = zeros(Int, maximum_generation(m)) would
    # be already a full iteration)
    intersections_per_generation = Int[]
    tgen = 0
    for g in @entities_in(m[Generation])
        if g.generation == 0
            continue
        end
        cur_maxgen = length(intersections_per_generation)
        if g.generation > cur_maxgen
            if g.generation - 1 > cur_maxgen
                append!(intersections_per_generation, zeros(Int, g.generation - 1 - cur_maxgen))
            end

            # This are pending intersections i.e. calculations that haven't ran yet
            if !(g in m[Unique]) && !(g in m[Done]) && !(g in m[Error])
                push!(intersections_per_generation, 1)
            end
        else
            # This are pending intersections i.e. calculations that haven't ran yet
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
    intersections = Tuple{Float64, Trial, Intersection, Generation}[]
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
                
               
            others = filter(unique_es) do x
                g2 = gencomp[x].generation
                return x.e != e1.e && g2 <= g1
            end
                
            ilock = ReentrantLock() 
            Threads.@threads for e2 in others
                for α in (0.25, 0.5, 0.75)
                    tstate = α * (e1.state + e2.state)

                    # Find closest prev unique
                    dist, minid = findmin(x -> Euclidean()(x.state, tstate), unique_es)
                    sink = unique_es[minid]
                    trial_to_sink = Euclidean()(sink.state, tstate)
                    if trial_to_sink > mindist
                        # from this unique, find all prev trials that led to it
                        funnel = filter(x -> x.origin != IntersectionMixed && sink.e != x.e && !isempty(x[Results].state.occupations) && Euclidean()(x[Results].state, sink.state) < 1e-2, @entities_in(m, Results && Trial))

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
                                push!(intersections, (dist, Trial(tstate, IntersectionMixed), Intersection(e1.e, e2.e), Generation(intersection_gen)))
                            finally
                                unlock(ilock)
                            end
                        end
                    end
                end
            end
            m[IntersectionSearcher][e1] = es[1]
        end

        n_int = length(intersections)
        
        Threads.@threads for i in 1:n_int
            int = intersections[i]
            min_to_others = minimum(j->Euclidean()(intersections[j][2].state, int[2].state), 1:i-1, init=int[1])
            min_to_prev_trials = minimum(x->Euclidean()(x.state, int[2].state), m[Trial])
            intersections[i] = (min(int[1], min(min_to_others, min_to_prev_trials)), int[2], int[3], int[4])
        end
            
        # Then insert n_new best intersections
        info = m[SearcherInfo][1]
        n_new = min(length(intersections), max(0, info.max_concurrent_trials - (info.n_running_calcs + info.n_pending_calcs)))
        new_intersections = sort!(intersections, rev=true, by = x -> x[1])
        n_added = 0
        gen = 0
        for i = 1:n_new
            dist, trial, intersection, generation = new_intersections[i]
            if dist > p.mindist * generation_distances[generation.generation - 1]
                gen = max(gen, generation.generation)
                add_search_entity!(m, search_e, generation, trial, intersection)
                n_added += 1
            end
        end
        # Remove entities that weren't used for intersecting so they might happen in the future
        for i = n_new+1:length(new_intersections)
            dist, trial, intersection, generation = new_intersections[i]
            trypop!(m[IntersectionSearcher], intersection.parent1)
        end
            
        if n_added != 0
            @debugv 1 "Added $n_added intersections at Generation($gen) with smallest distance $(new_intersections[n_added][1])"
        end
    end
end

