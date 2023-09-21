using Flux
using Zygote

const HUBTYPE = Vector{Vector{NamedTuple{(:id, :trace, :eigvals, :eigvecs, :occupations, :magmom), Tuple{Int64, NamedTuple{(:up, :down, :total), Tuple{Float64, Float64, Float64}}, NamedTuple{(:up, :down), Tuple{Vector{Float64}, Vector{Float64}}}, NamedTuple{(:up, :down), Tuple{Matrix{Float64}, Matrix{Float64}}}, NamedTuple{(:up, :down), Tuple{Matrix{Float64}, Matrix{Float64}}}, Float64}}}}
function prepare_data(l::Searcher)
    xs = Vector{Float64}[]
    ys = Vector{Float64}[]
    for e in @entities_in(l, Results)
        if e.constraining_steps > 0 && e.converged
            scf = joinpath(l, e, "scf.out")
            if ispath(scf)
                o = DFC.FileIO.qe_parse_pw_output(scf)
                if !haskey(o, :Hubbard)
                    continue
                end
                hub::HUBTYPE = o[:Hubbard]
                for i = e.constraining_steps:length(o[:Hubbard])-1
                    s1_ = State(o[:Hubbard][i])
                    s2_ = State(o[:Hubbard][end])
                    if any(x -> x < 0, s1_.eigvals[1]) || any(x -> x < 0, s2_.eigvals[1]) || any(x-> x < 0, diag(s1_.occupations[1])) || any(x-> x < 0, diag(s2_.occupations[1]))
                        continue
                    end
                    s1 = mat2features(s1_.occupations[1])
                    s2 = mat2features(s2_.occupations[1])
                       
                    if isempty(xs) || s1 != xs[end]
                        push!(xs, s1)
                        push!(ys, s2)
                    end
                end
                # if o[:converged]
                #     # duplicate the final steps
                #     for i = 1:100
                #         s1 = mat2features(State(o[:Hubbard][end]).occupations[1])
                #         s2 = mat2features(State(o[:Hubbard][end]).occupations[1])
                #         push!(xs, s1)
                #         push!(ys, s2)
                #     end
                # end
            end
        end
    end
    hcat(xs...), hcat(ys...)
end

function mlp_single_atom(n_features)
    # model = Chain(Dense(n_features, n_features, x -> leakyrelu(x, 0.2)),
    #               Dense(n_features, n_features, x -> leakyrelu(x, 0.2)),
    #               Dense(n_features, n_features, sigmoid))
    # second = div(n_features, 2)
    # third = div(n_features, 4)
    # model = Chain(Dense(n_features, second, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(second, third, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(third, second, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(second, n_features, x -> 2 * (sigmoid(x) - 0.5)),
    #               )
    model = Chain(Dense(rand(n_features, n_features).-0.5, true, x -> 2 * (sigmoid(x) - 0.5)),
                  Dense(rand(n_features, n_features).-0.5, true, x -> 2 * (sigmoid(x) - 0.5)),
                  # Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5)),
                  Dense(rand(n_features, n_features).-0.5, true, x -> 2 * (sigmoid(x) - 0.5)),
                  )
    return model
end

function mat2features(m::RomeoDFT.ColinMatrixType)
    n = size(m, 1)
    n_half = sum(1:n)
    n_features =n_half *2

    out = zeros(n_features)
    c = 1
    for i = 1:n
        for j = i:n
            out[c] = m.data[i, j] 
            out[c+n_half] = m.data[i, j+n] 
            c += 1
        end
    end
    out
end

function features2mat(v::AbstractVector)
    n = 0
    m = 0
    while 2m < length(v)
        m += n
        n += 1
    end
    n -= 1

    data = zeros(n, 2n)

    c = 1
    for i = 1:n
        for j = i:n
            data[i, j] = data[j,i] = v[c]
            data[i, j+n] = data[j, i+n] = v[c+m]
            c+=1
        end
    end
    return RomeoDFT.ColinMatrixType(data)
end

@pooled_component Base.@kwdef struct TrainerSettings
    n_iterations_per_training::Int
    n_points_per_training::Int
end

@pooled_component struct Model
    n_points::Int
    model_state
end

function train_model(l::Searcher)
    trainer_settings = l[TrainerSettings][1]
    X, y = prepare_data(l)
    ndat = size(X, 2)
    model = mlp_single_atom(size(X, 1))
    if !isempty(l[Model])
        Flux.loadmodel!(model, l[Model][end].model_state)
    end
    opt_state = Flux.setup(Adam(), model)
    
    @info "training on batch of size $ndat"
    train_set = [(X, y)]
    train_loss = []
    for i = 1:trainer_settings.n_iterations_per_training
        loss = 0
        suppress() do
            Flux.train!(model, train_set, opt_state) do m, x, y
                loss = Flux.Losses.mse(m(x), y)
            end
        end
        push!(train_loss, loss)
        if i % 100 == 0
            @debug i, loss
        end
    end
    
    return Model(ndat, Flux.state(model))
end

struct ModelTrainer <: System end
function Overseer.requested_components(::ModelTrainer)
    return (TrainerSettings, Intersection, Model)
end

function Overseer.update(::ModelTrainer, m::AbstractLedger)
    trainer_settings = m[TrainerSettings][1]
    
    prev_model = isempty(m[Model]) ? nothing : m[Model][end]
    
    n_points = sum(x -> x.converged ? x.niterations - x.constraining_steps : 0, m[Results], init=0)
    
    prev_points = prev_model === nothing ? 0 : prev_model.n_points

    if n_points - prev_points > trainer_settings.n_points_per_training
        
        model = train_model(m)
        
        Entity(m, m[Template][entity(m[SearcherInfo],1)], model, Generation(length(m[Model].c.data)+1))
    end
end

struct MLTrialGenerator <: System end
function Overseer.requested_components(::MLTrialGenerator)
    return (TrainerSettings, Intersection, Model)
end
    
function Overseer.update(::MLTrialGenerator, m::AbstractLedger)
    # if no model yet, skip
    model = isempty(m[Model]) ? nothing : m[Model][end]
    model === nothing && return
    
    model_e = last_entity(m[Model])
    
    U = getfirst(x->x.dftu.U != 0, m[Template][1].structure.atoms).dftu.U

    # check how many calcs are pending on server
    # if free room -> generate some more calcs with the current model
    info = m[SearcherInfo][1]
    max_new = max(0, info.max_concurrent_trials - (info.n_running_calcs + info.n_pending_calcs))
    max_new <= 0 && return 
 
    # find previous unique, intersection, other trials about to run
    # to check whether it's duplication
    unique_states = @entities_in(m, Unique && Results)
    pending_states = @entities_in(m, Trial && !Results)

    # Generate Trial with origin ModelOptimized + increment generation
    # Severely limit the amount of new intersections, maybe only 6 per generation or whatever
    #  new trial -> hopefully new unique -> 3 x whatever prev unique intersections -> train model whenever Y new unique states are found from those -> rinse repeat 

    str = m[Template][1].structure
    natoms = length(filter(ismagnetic, str.atoms))
    nshell = size(m[Results][1].state.occupations[1], 1)
    dist_thr = 0.5
    curgen = maximum_generation(m)
    
    base_e = entity(m[BaseCase], 1)
    base_state = m[Results][base_e].state
    nelec    = round.(Int, base_state.totoccs)
    norb     = size.(base_state.occupations, 1)
    
    flux_model = mlp_single_atom(30)
    Flux.loadmodel!(flux_model,model.model_state)
    # TODO multithreading
    lck = ReentrantLock()
    max_reached = false
    while max_new > 0 && !max_reached
        Threads.@threads for i = 1:Threads.nthreads()
            # random start
            min_e = Entity(0)
            s = nothing
            n_tries = 0 
            while  n_tries < 1000
                min_dist = Inf
                x0 = mat2features(rand_trial(norb, nelec).state.occupations[1])

                s = State([features2mat(flux_model(x0))])

                if any(x->x<0, diag(s.occupations[1])) || any(x -> x < 0, s.eigvals[1])
                    continue
                end
            
                for e in @entities_in(m, Unique && Results)
                    dist = Euclidean()(e.state, s)
                    if dist < min_dist
                        min_dist = dist
                        min_e = Entity(e)
                    end
                end

                for e in @entities_in(m, Trial)
                    dist = Euclidean()(e.state, s)
                    if dist < min_dist
                        min_dist = dist
                        min_e = Entity(e)
                    end
                end
                @show min_dist
                if min_dist > 1 && max_new > 0  
                    lock(lck) do 
                        trial = Trial(s, RomeoDFT.ModelOptimized) # TODO other tag?

                        # add new entity with optimized occ
                        new_e = add_search_entity!(m, model_e, trial, m[Generation][model_e])
                        m[Model][new_e] = model_e
                        info.n_pending_calcs += 1
                        max_new -= 1
                    end
                else 
                    n_tries += 1
                end
            end
            if n_tries == 1000
                max_reached = true
            end
        end
    end
end

## Save for later
# abstract type AbstractSettings end
# @component struct HighLevelSettings <: AbstractSettings
#     #
#     #
#     # 
# end

# @component struct HighLevelResults
# end

# result_comp(::Type{HighLevelSettings}) = HighLevelResults

# struct HighLevel <: System end

# function Overseer.update(::HighLevel, l::Searcher)
    
#     for e in @entities_in(l, HighLevelSettings)
        
#         if !(e in l[SCFSettings])
#             l[e] = SCFSettings(e.create_scf_settings)
#         elseif e in l[SCFSettings] && e in l[Results] && !(e in l[NSCFSettings])
#             # do something with e
#             l[e] = NSCFSettings(e.create_nscf_settings)
#         elseif e in l[NSCFSettings] && e in l[NSCFResults] && !(e in l[BandsSettings])
#             # do something with nscf output
#             l[e] = BandsSettings(e.create_bands_settings)
#         elseif e in l[BandsResults]
#             l[e] = HighLevelResults(x, y, z)
#         end
#     end
# end

# struct DoneChecker <: System end

# function Overseer.update(::DoneChecker, l::Searcher)
#     dones = Dict()
#     done_c = l[Done]
#     for c in components(l, AbstractSettings)
#         done_comp = l[result_comp(eltype(c))]

#         for e in @entities_in(c && !done_c)
#             dones[Entity(e)] = get!(dones, Entity(e), true) && e in done_comp
#         end
#     end

#     for (e, d) in dones
#         if d
#             done_c[e] = Done(false)
#         end
#     end
# end
