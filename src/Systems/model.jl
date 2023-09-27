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
        out = DFC.FileIO.qe_parse_pw_output(path)
        if !haskey(out, :Hubbard)
            continue
        end
        hubbard::HUBTYPE = out[:Hubbard]
        
        r = get(out, :Hubbard_iterations, 1):length(hubbard)-1
        last_state = State(hubbard[end])
        
        xs = Vector{Vector{Float32}}(undef, length(r))
        Threads.@threads for i in r
            xs[i] = mat2features(State(hubbard[i]).occupations[1])
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

const MAX_DATA = 5000

function train_model(l::Searcher, n_points)
    trainer_settings = l[TrainerSettings][1]
    X, y, x_conv, y_conv = prepare_data(l)
    ndat = size(X, 2)
    
    model      = mlp_single_atom(size(X, 1))
    model_conv = mlp_converge(size(X, 1))
    
    opt_state      = Flux.setup(Adam(), model)
    opt_state_conv = Flux.setup(Adam(), model_conv)

    train_ratio      = clamp(MAX_DATA / size(X, 2), 0.0, 1.0)
    train_ratio_conv = clamp(MAX_DATA / size(x_conv, 2), 0.0, 1.0)

    train, test           = Flux.splitobs((X, y),           at = train_ratio, shuffle=true)

    r_non_conv = findall(iszero, y_conv[1, :])
    r_conv = union(r_non_conv, findall(!iszero, y_conv[1, :])[1:min(2*length(r_non_conv), size(y_conv,2)-length(r_non_conv))])
    train_conv = (x_conv[:, r_conv], y_conv[:, r_conv])
    
    @debug length(findall(iszero, train_conv[end])), length(train_conv[end])
    
    if size(X, 2)>1000 && !isempty(test[1])&& !isempty(l[Model])
        Flux.loadmodel!(model,      l[Model][end].model_state)
        Flux.loadmodel!(model_conv, l[Model][end].model_state_converge)
    end
    
    loss      = Inf
    d_loss = Inf
    loss_conv = Inf
    d_loss_conv = Inf
    time = now()
    Threads.@sync begin
        i1 = 1
        i2 = 1
        Threads.@spawn while loss > 1e-3 && d_loss > 1e-6 && now() - time < Minute(1)
            suppress() do
                Flux.train!(model, [train], opt_state) do m, x, y
                    t_loss = Flux.Losses.mse(m(x), y)
                    d_loss = abs(loss - t_loss)
                    loss = t_loss
                end
            end
            if i1 % 100 == 0
                @debug "i1 $i1, loss $loss"
            end
            i1 += 1
        end
        Threads.@spawn while loss_conv > 1e-3 && d_loss_conv > 1e-6 && now() - time < Minute(1)
            suppress() do
                Flux.train!(model_conv, [train_conv], opt_state_conv) do m, x, y
                    t_loss = Flux.Losses.mse(m(x), y)
                    d_loss_conv = abs(loss_conv - t_loss)
                    loss_conv = t_loss
                end
            end
            if i2 % 100 == 0
                @debug "i2 $i2, loss_conv $loss_conv"
            end
            i2 += 1
        end
    end
    return Model(!isempty(test[1]) ? l[Model][end].n_points : n_points, Flux.state(model), Flux.state(model_conv))
end

struct ModelTrainer <: System end
function Overseer.requested_components(::ModelTrainer)
    return (TrainerSettings, Intersection, Model)
end

function Overseer.update(::ModelTrainer, m::AbstractLedger)
    if length(m[Results]) < 10
        return
    end
    n_unique, n_total = unique_evolution(m)
    ilast = length(m[Model].c.data) + 1
    
    if length(m[Results]) < 10
        return
    elseif n_total[ilast] != 0 && n_unique[ilast] / n_total[ilast] > 0.4
        return
    end
    trainer_settings = m[TrainerSettings][1]
    
    prev_model = isempty(m[Model]) ? nothing : m[Model][end]
    
    n_points = sum(x -> x.converged ? x.niterations - x.constraining_steps : 0, m[Results], init=0)
    
    prev_points = prev_model === nothing ? 0 : prev_model.n_points

    if n_points > prev_points * trainer_settings.new_train_data_ratio
        
        model = train_model(m, n_points)
        
        Entity(m, m[Template][1], model, Generation(length(m[Model].c.data)+2))
    end
end

@component Base.@kwdef struct MLTrialSettings
    n_tries::Int = 10000
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
    
    flux_model = mlp_single_atom(30)
    Flux.loadmodel!(flux_model, model.model_state)

    model_conv = mlp_converge(30)
    Flux.loadmodel!(model_conv, model.model_state_converge)

    return flux_model, model_conv
end

converge_chance(model, s::State) = model(mat2features(s.occupations[1]))[1]

function Overseer.update(::MLTrialGenerator, m::AbstractLedger)
    # if no model yet, skip
    flux_model, model_conv = model(m)
    
    flux_model === nothing && return
    model_e = last_entity(m[Model])
    
    unique_states = @entities_in(m, Unique && Results)
    pending_states = @entities_in(m, Trial && !Results)

    n_new = 0
    max_tries = m[MLTrialSettings][1].n_tries
    dist_thr  = m[MLTrialSettings][1].minimum_distance
    
    while max_new(m) > 0
        max_dist = 0
        new_s = nothing
        conv_chance = 0f0
        for _ = 1:max_tries
            
            s = State(flux_model, rand_trial(m)[1].state)
            
            c_chance = converge_chance(model_conv, s)
            
            if c_chance < 0.6f0 || any(x->x<0, diag(s.occupations[1])) || any(x -> x < 0, s.eigvals[1])
                continue
            end
        
            min_dist = Inf
            for e in @entities_in(m, Unique && Results)
                dist = Euclidean()(e.state, s)
                if dist < min_dist
                    min_dist = dist
                    # min_e = Entity(e)
                end
            end

            for e in @entities_in(m, Trial)
                dist = Euclidean()(e.state, s)
                if dist < min_dist
                    min_dist = dist
                    # min_e = Entity(e)
                end
            end
            
            if min_dist > dist_thr
                max_dist = min_dist
                new_s = s
                conv_chance = c_chance
                break
            else
                if min_dist > max_dist
                    max_dist = min_dist
                    new_s = s
                    conv_chance = c_chance
                end
            end
        end
        @debug "New ml trial: max dist = $max_dist, conv chance = $conv_chance"
        trial = Trial(new_s, RomeoDFT.ModelOptimized) # TODO other tag?
        # add new entity with optimized occ
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
