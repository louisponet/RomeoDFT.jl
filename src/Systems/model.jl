using Flux
using Zygote
const HUBTYPE = Vector{Vector{NamedTuple{(:id, :trace, :eigvals, :eigvecs, :occupations, :magmom), Tuple{Int64, NamedTuple{(:up, :down, :total), Tuple{Float64, Float64, Float64}}, NamedTuple{(:up, :down), Tuple{Vector{Float64}, Vector{Float64}}}, NamedTuple{(:up, :down), Tuple{Matrix{Float64}, Matrix{Float64}}}, NamedTuple{(:up, :down), Tuple{Matrix{Float64}, Matrix{Float64}}}, Float64}}}}

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
    # model = Chain(Dense(rand(n_features, n_features).-0.5, true, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(rand(n_features, n_features).-0.5, true, x -> 2 * (sigmoid(x) - 0.5)),
    #               # Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(rand(n_features, n_features).-0.5, true, x -> 2 * (sigmoid(x) - 0.5)),
    #               )
    model = Chain(Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5f0)),
                  Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5f0)),
                  # Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5)),
                  Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5f0)),
                  )
    # model = Chain(Dense(n_features, div(n_features,2), x -> 2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,2), div(n_features, 4), x -> 2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features, 4), div(n_features,2), x -> 2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,2), n_features, x -> 2 * (sigmoid(x) - 0.5f0)),
    #               )
    return model
end

function mlp_vae(n_features)
    # model = Chain(Dense(n_features, div(n_features, 2), x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,2), div(n_features, 4), x->2 * (sigmoid(x) - 0.5f0)),
    #               # Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(div(n_features,4), 2, x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(2, div(n_features,4), x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,4),div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,2),n_features, x->2 * (sigmoid(x) - 0.5f0)),
    #               )
    model = Chain(Dense(n_features, div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,2), div(n_features,4), x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,4), 5, sigmoid),
                  Dense(5, div(n_features,4), sigmoid),
                  Dense(div(n_features,4), div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
                  # Dense(n_features,n_features, x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,2), n_features, x->2 * (sigmoid(x) - 0.5f0)),
                  )
    return model
end
function mlp_encoder(n_features)
    # model = Chain(Dense(n_features, div(n_features, 2), x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,2), div(n_features, 4), x->2 * (sigmoid(x) - 0.5f0)),
    #               # Dense(n_features, n_features, x -> 2 * (sigmoid(x) - 0.5)),
    #               Dense(div(n_features,4), 2, x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(2, div(n_features,4), x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,4),div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
    #               Dense(div(n_features,2),n_features, x->2 * (sigmoid(x) - 0.5f0)),
    #               )
    model = Chain(Dense(n_features, div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,2), div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,2), 2, sigmoid),
                  )
    return model
end

function mlp_decoder(n_features)
    model = Chain(Dense(2, div(n_features,4), x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,4),div(n_features,2), x->2 * (sigmoid(x) - 0.5f0)),
                  Dense(div(n_features,2),n_features, x->2 * (sigmoid(x) - 0.5f0)),
                  )
    return model
end

function mlp_converge(n_features)
    Chain(Dense(n_features, div(n_features, 2), x -> leakyrelu(x, 0.2f0)),
          Dense(div(n_features, 2), div(n_features, 4), sigmoid),
          Dense(div(n_features, 4), 1, sigmoid))
end
    

function mat2features(m::RomeoDFT.ColinMatrixType)
    n = size(m, 1)
    n_half = sum(1:n)
    n_features =n_half *2

    out = zeros(Float32, n_features)
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

State(model, s::State) = State([features2mat(model(mat2features(s.occupations[1])))])

@component struct ModelData
    x::Matrix{Float32}
    y::Matrix{Float32}
end

struct ModelDataExtractor <: System end

Overseer.requested_components(::ModelDataExtractor) = (ModelData, Results)

function Overseer.update(::ModelDataExtractor, l::AbstractLedger)
    @error_capturing_threaded for e in @entities_in(l, Results && !ModelData)
        path = joinpath(l, e, "scf.out")
        if !ispath(path)
            continue
        end
        
        out = DFC.FileIO.qe_parse_pw_output(path)
        if !haskey(out, :Hubbard)
            continue
        end
        hubbard::HUBTYPE = out[:Hubbard]
        
        r = get(out, :Hubbard_iterations, 1):length(hubbard)-1
        last_state = State(hubbard[end])
        
        states = Vector{State}(undef, length(r))
        Threads.@threads for i in r
            states[i] = State(hubbard[i])
        end

        good_r = findall(x->minimum(x.eigvals[1]) > -0.01 && maximum(x.eigvals[1]) < 1.01, states)
        
        xs = Vector{Vector{Float32}}(undef, length(good_r))
        Threads.@threads for i in 1:length(good_r)
            xs[i] = mat2features(states[good_r[i]].occupations[1])
        end

        isunique = trues(length(xs))
        @inbounds for i in 1:length(xs)
            
            if isunique[i]
                xsi = xs[i]
                for j = i+1:length(xs)
                    if isunique[j]
                        xsj = xs[j]
                        
                        isunique[j] = sum(k -> abs(xsi[k] - xsj[k]), 1:length(xsi)) > 0.1
                    end
                end
            end
        end

        x = xs[isunique]
        
        if e.converged
            ys = fill(mat2features(last_state.occupations[1]), length(x))
        else
            # not the most performant
            ys = fill(0.0f0, length(x))
        end
        #TODO can we just 1 ys ?
        l[ModelData][e] = ModelData(reduce(hcat, x), reduce(hcat, ys))
    end
end
const MAX_DATA = 100000

# Actually maybe it's better to always use the full non converged samples since they are less 
function prepare_data(l::Searcher)
    xs          = Matrix{Float32}[]
    xs_converge = Matrix{Float32}[]
    ys          = Matrix{Float32}[]
    ys_converge = Matrix{Float32}[]
    for e in @entities_in(l, ModelData)
        if size(e.y, 1) == 1
            push!(xs_converge, e.x) 
            push!(ys_converge, e.y)
        else
            push!(xs_converge, e.x)
            push!(ys_converge, fill(1.0f0, size(e.x, 2))')
            push!(xs, e.x)
            push!(ys, e.y)
        end
    end
    reduce(hcat, xs), reduce(hcat, ys), reduce(hcat, xs_converge), reduce(hcat, ys_converge)
end

@pooled_component Base.@kwdef struct TrainerSettings
    n_iterations_per_training::Int
    # This determines how much additional data w.r.t. previous there needs to be, like 1.2= 20% more data
    new_train_data_ratio::Float64
end

@pooled_component struct Model
    n_points::Int
    model_state
    model_state_converge
end

# Things to try:
#   - boosting by learning on data that's predicted worst
#   - only update model if new one is better than old one

function train_model(l::Searcher; model = mlp_vae(30), max_time = Minute(1))
    trainer_settings = l[TrainerSettings][1]
    X, y, x_conv, y_conv = prepare_data(l)
    ndat = size(X, 2)
    
    model_conv = mlp_converge(size(X, 1))
    
    opt_state      = Flux.setup(Adam(), model)
    opt_state_conv = Flux.setup(Adam(), model_conv)

    unique_features = map(x -> mat2features(x.state.occupations[1]), @entities_in(l, Unique && Results))
    unique_x = Vector{Float32}[]
    unique_y = Vector{Float32}[]
    for i in 1:length(unique_features)
        for j in 1:length(unique_features)
            if i == j
                continue
            end
            push!(unique_x, unique_features[i])
            push!(unique_y, unique_features[j])
        end
    end

    u_x = reduce(hcat, unique_x)
    u_y = reduce(hcat, unique_y)
    
    train_conv = Flux.eachobs((x_conv, y_conv), batchsize=500)
    train = Flux.eachobs((X, y), batchsize=500)

    loss        = Inf
    loss_uneq = Inf
    d_loss      = Inf
    d_loss_uneq = Inf
    loss_conv   = Inf
    d_loss_conv = Inf
    time = now()


    rand = reduce(hcat, tmap(x->mat2features(x.state.occupations[1]), rand_trial(l, 1000)))
    
    i1 = 1
    while loss > 1e-5 && (d_loss > 1e-6 || d_loss_uneq > 1e-6) && now() - time < max_time
        losses_grads = tmap(enumerate(train)) do (i, (x, y))
            val, grads = Flux.withgradient(model) do m
                zx = m[1:3](x)
                zy = m[1:3](y)
                Flux.Losses.mse(m(x), x)/10 + Flux.Losses.mse(zx, zy)
            end
            val, grads[1]
        end

        tl = 0.0
        for (v, g) in losses_grads
            tl += v
            Flux.update!(opt_state, model, g)
        end
        tl /= length(train)

        d_loss = abs(loss - tl)
        loss = tl
        if i1 % 100 == 0
            @debug "i1 $i1\nloss_equal $loss"
        end
        
        Flux.train!(model, [(rand, rand)], opt_state) do m, x, y
            Flux.Losses.mse(m(x), y)
        end
        
        if i1 %100 ==0
            tl = Flux.Losses.mse(model(rand), rand)
            @debug "loss_rand $tl"
        
        end
        if i1 % 10 == 0
            tl, grads = Flux.withgradient(model) do m
                m_x = m[1:3](u_x)
                m_y = m[1:3](u_y)
                
                sum(1:size(m_y,2)) do i
                    tot = 0f0
                    for d in 1:size(m_x, 1)
                        tot += min((m_x[d, i] - m_y[d, i])^2, (m_x[d, i] - m_y[d, i] + 1)^2)
                    end
                    1/sqrt(tot)
                end/4000
            end
            d_loss_uneq = abs(loss_uneq - tl)
            loss_uneq = tl

            Flux.update!(opt_state, model, grads[1])
            if i1 %100 ==0
                @debug "loss_unequal $loss_uneq"
            end
        end
        i1 += 1
    end
    
    # time = now()
    # i2 = 1
    # while loss_conv > 1e-5 && d_loss_conv > 1e-6 && now() - time < Minute(1)
    #     losses_grads = tmap(enumerate(train_conv)) do (i, (x, y))
    #         val, grads = Flux.withgradient(model_conv) do m
    #             Flux.Losses.mse(m(x), y)
    #         end
    #         val, grads[1]
    #     end

    #     tl = 0.0
    #     for (v, g) in losses_grads
    #         tl += v
    #         Flux.update!(opt_state_conv, model_conv, g)
    #     end
    #     tl /= length(train_conv)

    #     d_loss_conv = abs(loss_conv - tl)
    #     loss_conv = tl
    #     if i2 % 100 == 0
    #         @debug "i2 $i2, loss_conv $loss_conv"
    #     end
    #     i2 += 1
    # end
    return model

    # return Model(now() - time > Minute(2) && !isempty(l[Model]) ? l[Model][end].n_points : n_points, Flux.state(model), Flux.state(model_conv))
end

struct ModelTrainer <: System end
function Overseer.requested_components(::ModelTrainer)
    return (TrainerSettings, Intersection, Model)
end

function Overseer.update(::ModelTrainer, m::AbstractLedger)
    # No training till unique larger than 10
    if length(m[Unique]) - 1 < 10
        return
    end
    
    trainer_settings = m[TrainerSettings][1]
    
    prev_model = isempty(m[Model]) ? nothing : (m_ = mlp_vae(30); Flux.loadmodel!(m_, m[Model][end].model_state); m_)
    
    n_points = sum(x -> x.converged ? x.niterations - x.constraining_steps : 0, m[Results], init=0)
    
    prev_points = prev_model === nothing ? 0 : m[Model][end].n_points
    curgen = maximum_generation(m)
    if prev_model === nothing || sum(x -> x.converged && x.generation == curgen ? x.niterations : 0, @entities_in(m, Results &&  Generation)) > 1000  
        
        # model = prev_model === nothing ? train_model(m) : train_model(m, model=prev_model)
        model = train_model(m)
        
        Entity(m, m[Template][1], Model(length(m[Unique]), Flux.state(model), Flux.state(mlp_converge(30))), Generation(length(m[Model].c.data)+2))
    end
end

@component Base.@kwdef struct MLTrialSettings
    n_tries::Int = 1000
    minimum_distance::Float64 = 1.0
end

struct MLTrialGenerator <: System end
function Overseer.requested_components(::MLTrialGenerator)
    return (MLTrialSettings, Model, SearcherInfo)
end

function model(m::AbstractLedger)
    # if no model yet, skip
    model = isempty(m[Model]) ? nothing : m[Model][end]
    
    model === nothing && return nothing, nothing
    
    flux_model = mlp_vae(30)
    Flux.loadmodel!(flux_model, model.model_state)

    model_conv = mlp_converge(30)
    Flux.loadmodel!(model_conv, model.model_state_converge)

    return flux_model, model_conv
end

converge_chance(model, s::State) = model(mat2features(s.occupations[1]))[1]

function decoder(m::AbstractLedger)
    flux_model, model_conv = model(m)
    flux_model === nothing && return
    n_layers = length(flux_model) 
    return flux_model[div(n_layers, 2)+1:end]
end    
function encoder(m::AbstractLedger)
    flux_model, model_conv = model(m)
    flux_model === nothing && return
    n_layers = length(flux_model) 
    return flux_model[1:div(n_layers, 2)]
end    

function potential_trials(m::AbstractLedger; en = encoder(m), de = decoder(m))
    en === nothing && return
    unique_states = @entities_in(m, Unique && Results)
    pending_states = @entities_in(m, Trial && !Results)

    max_tries = m[MLTrialSettings][1].n_tries
    dist_thr  = m[MLTrialSettings][1].minimum_distance

    latent_space = map(x->en(mat2features(x.state.occupations[1])), unique_states)
    trial_latent_space = Vector{Float64}[]
    for e in @entities_in(m, ModelData)
        append!(trial_latent_space, tmap(i -> en(view(e.x, :, i)), 1:size(e.x, 2)))
    end
    time = now()
    
    d_latent = length(trial_latent_space[1])
    n_latent = length(latent_space)
    
    lck = ReentrantLock()
    potential_intersections = Tuple{Float64, Vector{Float32}}[]    
    Threads.@threads for i in 1:n_latent
        z1 = latent_space[i]
        for j = i+1:n_latent
            z2 = latent_space[j]

            mid = (z1 .+ z2) ./ 2
            for i = 1:d_latent
                mid[i] = mid[i] > 1 ? mid[i] - 1 : mid[i]
            end
            # res = optimize([0f0,0f0], [1f0,1f0],mid, Optim.Fminbox(ParticleSwarm(n_particles=10)), Optim.Options(show_trace=true)) do x
            res = optimize(mid,ParticleSwarm(lower=fill(0f0, d_latent), upper= fill(1f0, d_latent), n_particles=10), Optim.Options(iterations=200,show_trace=false)) do x
                s = State([features2mat(de(x))])
                if minimum(s.eigvals[1]) < -0.01 || maximum(s.eigvals[1])>1.01
                    return typemax(Float32)
                end
                out = sum(1:length(latent_space)) do i
                    tot = 0f0
                    for d in 1:d_latent
                        tot += (latent_space[i][d] - x[d])^2
                    end
                    1/sqrt(tot)
                end

                t = sum(1:length(trial_latent_space)) do i
                    tot = 0f0
                    for d in 1:d_latent
                        tot += (trial_latent_space[i][d] - x[d])^2
                    end
                    1/sqrt(tot)
                end
                t2 = sum(1:length(potential_intersections), init=0.0f0) do i
                    tot = 0f0
                    for d in 1:d_latent
                        tot += (potential_intersections[i][2][d] - x[d])^2
                    end
                    1/sqrt(tot)
                end
                
                return out + t + 100*t2
            end
            coord = res.minimizer
            for i = 1:d_latent
                coord[i] = coord[i] > 1 ? coord[i] - 1 : coord[i]
            end
            
            lock(lck) do
                push!(potential_intersections, (res.minimum, coord))
            end
        end
    end
    return unique(x -> x[1], potential_intersections)
end

function Overseer.update(::MLTrialGenerator, m::AbstractLedger)
    if max_new(m) <= 0
        return
    end
    potential = potential_trials(m)
    potential === nothing && return
    de = decoder(m)
    model_e = last_entity(m[Model])
    
    n_new = 0
    intersections = sort(potential, by = x->x[1], rev=true)
    while max_new(m) > 0 && !isempty(intersections)
        dist, mid = pop!(intersections)

        s = State([features2mat(de(mid))])
        
        @debug "New ml trial: min dist = $dist"
        # add s to latent space
        trial = Trial(s, RomeoDFT.ModelOptimized) # TODO other tag?
        new_e = add_search_entity!(m, model_e, trial, m[Generation][model_e])
        m[Model][new_e] = model_e
        n_new += 1
    end
    if n_new != 0
        @debug "$n_new new ML trials at Generation($(m[Generation][model_e].generation))" 
    end
end

## This is super bad
struct MLIntersector <: System end
function Overseer.requested_components(::MLIntersector)
    return (TrainerSettings, Intersection, Model)
end
    
function Overseer.update(::MLIntersector, m::AbstractLedger)
    # if no model yet, skip
    flux_model = model(m)
    flux_model === nothing && return
    model_e = last_entity(m[Model])
    
    unique_states = collect(@entities_in(m, Unique && Results))
    pending_states = collect(@entities_in(m, Trial && !Results))

    max_dist = 0
    n_new = max_new(m)
    n_new <= 0 && return
    n_tries = 0
    max_tries = m[MLTrialSettings][1].n_tries
    dist_thr  = m[MLTrialSettings][1].minimum_distance

    lck = ReentrantLock()
    # Trial get the ones that have dist lower than dist_thr
    trial_intersections = Tuple{Float64, Trial, Intersection}[]
    # dist > dist_thr
    good_intersections  = Tuple{Float64, Trial, Intersection}[]
    
    # if length(good) > max_new -> stop and ship it
    Threads.@threads for e1 in unique_states
        if n_new <= 0
            break
        end
        for e2 in unique_states
            if n_new <= 0
                break
            end
            if e1 == e2 
                continue
            end
            for α in (0.25, 0.5, 0.75)
                if n_new <= 0
                    break
                end
                tstate = α * e1.state + (1-α) * e2.state
                tstate = State([features2mat(flux_model(mat2features(tstate.occupations[1])))])
                
                dist, minid   = isempty(unique_states) ? (Inf, 0) : findmin(x -> Euclidean()(x.state, tstate), unique_states)
                dist2, minid2 = isempty(pending_states) ? (Inf, 0) : findmin(x -> Euclidean()(x.state, tstate), pending_states)
                lock(lck) do
                    dist3, mid3   = isempty(good_intersections) ? (Inf, 0) : findmin(x -> Euclidean()(x[2].state, tstate), good_intersections)
                    dist = min(dist, dist2, dist3)
                end


                if dist > dist_thr
                    lock(lck) do
                        push!(good_intersections, (dist, Trial(tstate, ModelOptimized), Intersection(Entity(e1), Entity(e2))))
                        n_new -= 1
                    end
                else
                    lock(lck) do
                        push!(trial_intersections, (dist, Trial(tstate, ModelOptimized), Intersection(Entity(e1), Entity(e2))))
                    end
                end
            end
        end
    end
    while n_new > 0 && !isempty(trial_intersections)
        sort!(trial_intersections, by = x -> min(x[1], minimum(y->Euclidean()(x[2].state, y[2].state), good_intersections, init=Inf)))
        push!(good_intersections, pop!(trial_intersections))
        n_new-=1
    end
    for (_, trial, intersection) in good_intersections 
        new_e = add_search_entity!(m, model_e, trial, m[Generation][model_e], intersection)
        m[Model][new_e] = model_e
    end
    @debug "$(length(good_intersections)) new ML trials at Generation($(m[Generation][model_e].generation))" 
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
