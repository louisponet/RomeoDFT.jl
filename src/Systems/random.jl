struct RandomTrialGenerator <: System end
function Overseer.requested_components(::RandomTrialGenerator)
    return (RandomSearcher, Intersection, BaseCase)
end

function rand_trial(l::Searcher, n=1)
    base_e     = entity(l[BaseCase], 1)
    base_state = l[Results][base_e].state
    nelec      = round.(Int, base_state.totoccs)
    norb       = size.(base_state.occupations, 1)

    # Here we check whether different oxidation states have been tried
    # if yes we start randomly taking the total occupations for eigenvalues
    diff_ox_tried = !isempty(l[Trial]) && round(Int, maximum(x -> sum(x.state.totoccs), l[Trial])) != round(Int, minimum(x->sum(x.state.totoccs), l[Trial]))

    out = Trial[]
    for i = 1:n
        if !diff_ox_tried
            push!(out, rand_trial(norb, nelec))
        else
            eigvals = map(1:length(nelec)) do i
                e = Entity(l[Unique], rand(1:length(l[Unique])))
                
                while !(e in l[Results])
                    e = Entity(l[Unique], rand(1:length(l[Unique])))
                end
            
                l[Results][e].state.eigvals[i]
            end
            
            push!(out, rand_trial(eigvals))
        end
    end
    return out
end

function rand_trial(eigvals::Vector)
    n_orb = div(length(eigvals[1]), 2)
    nangles = div(n_orb * (n_orb - 1), 2)
    
    rand_angles() = Angles([π * (rand() - 0.5) for i in 1:nangles])
    
    occs = map(eigvals) do eig
        D = DFWannier.MagneticVector(eig)
        V = DFWannier.ColinMatrix(Matrix(rand_angles()), Matrix(rand_angles()))
        return Matrix(Eigen(D, V))
    end
    return Trial(State(occs), RandomMixed)
end
    
function rand_trial(n_orb_per_at::Vector, n_elec_per_at::Vector)
    
    occs = map(zip(n_orb_per_at, n_elec_per_at)) do (norb, nelec)
        nangles = div(norb * (norb - 1), 2)
        
        rand_angles() = Angles([π * (rand() - 0.5) for i in 1:nangles])
        
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
    if isempty(m[BaseCase]) || !all_children_done(m, base_e) || base_e ∉ m[Results]
        return
    end
    
    base_state = m[Results][base_e].state
    if isempty(base_state.occupations)
        @error "Something went wrong with the basecase calculation"
        return
    end
    rand_search_e    = entity(m[RandomSearcher], 1)

    maxgen = maximum_generation(m)
    n_new = max_new(m)
    for trial in rand_trial(m, n_new)
        e = add_search_entity!(m, rand_search_e,
                               trial,
                               Generation(maxgen))
                               
        if Hybrid in m && length(m[Hybrid]) != 0
            m[e] = Hybrid()
        end
        n_new += 1
        
    end
    if n_new != 0
        @debug "$n_new random trials at Generation($(maxgen))."
    end
end
