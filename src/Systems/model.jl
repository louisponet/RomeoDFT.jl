using MCMCChains, Turing, Distributions, StatsBase, Combinatorics
using Optim
using Optim: minimizer, LineSearches 

function E_hund_exchange(diag_occs, J)
    # -J_H Σ (1/2 + S_α * S_β)
    s = 0.
    n = size(diag_occs, 1)
    combs = combinations(1:n, 2)
    @inbounds for i in 1:n
        o = diag_occs[i]
        for (p,q) in combs
            s += o[p, 1] * o[q, 1] + o[p, 2] * o[q, 2]
        end
    end
    return -J * s
end

function E_U(n, U)
    s = 0.0
    d = size(n, 1)
    @inbounds for m1 = 1:d
        s += n[m1,m1]
        s += n[m1,m1+d]
        for m2 = 1:d
            s -= n[m1, m2] *  n[m2, m1]
            s -= n[m1, m2+d] *  n[m2, m1+d]
        end
    end
    return 0.5 * U * s
end

function on_site_energy(Jh, β, C, constant_shift, data, T)    
    magmoms_squared_sum = data.magmoms_squared_sum
    Eu = data.E_Us
    diag_occs = data.diag_occs
    diag_occs_sum = data.diag_occs_sum
    diag_occs_squared_sum = data.diag_occs_squared_sum

    natoms = Int(length(diag_occs_squared_sum[1]) / length(C))

    etot = zeros(T, length(magmoms_squared_sum))
    @inbounds Threads.@threads for i = 1:length(magmoms_squared_sum)
        m = magmoms_squared_sum[i]
        diag_occ_square = diag_occs_squared_sum[i]
        diag_occ_sum = diag_occs_sum[i]
        diag_occ = diag_occs[i]
        Eu_ = Eu[i]  # +U energy
        Eh = E_hund_exchange(diag_occ, Jh)         # Hund's exchange
        Ecf = sum(y->y[1] * y[2], zip(diag_occ_square, repeat(C, natoms)))
        
        E_colombic = - β * sum(diag_occ_sum)       # e-N interactions
        etot[i] = Eh + Ecf + Eu_ + E_colombic
    end
    return etot .+ constant_shift * length(C)
end

function prepare_data(l::Searcher)

    base_energy = l[Results][entity(l[BaseCase],1)].total_energy
    
    nmagats = length(filter(ismagnetic, l[Template][1].structure.atoms))
    
    U = getfirst(x->x.dftu.U != 0, l[Template][1].structure.atoms).dftu.U
    
    conv_res = @entities_in(l, Results && (Unique || Intersection))
    
    states = map(x->x.state, conv_res)
    energies = map(x->(x.total_energy - base_energy) * 13.6, conv_res)
    occupations = map(x->x.occupations, states)
    return (; energies, prepare_data(occupations, U)...)
end

# TODO duplicate code
function prepare_data(occs, U)
    magmoms = tmap(occs) do occupation
        map(occupation) do occ
            tr(view(occ, Up())) - tr(view(occ, Down())) 
        end
    end
    magmoms_squared_sum   = tmap(magmom -> sum(m -> m^2, magmom), magmoms)
    E_Us                  = tmap(occ -> sum(m -> E_U(m, U), occ), occs)
    diag_occs             = tmap(occ -> map(x -> hcat(diag(view(x, Up())), diag(view(x, Down()))), occ), occs)
    diag_occs_sum         = tmap(occ -> sum(x -> (diag(view(x, Up())) .+ diag(view(x, Down()))), occ), occs)
    diag_occs_squared_sum = tmap(x-> x.^2, diag_occs_sum)
    return (; magmoms_squared_sum, diag_occs_sum, E_Us, diag_occs_squared_sum, diag_occs)
end

@model function total_energy(data, ::Type{T}=Float64) where {T}
    # physics
    Jh ~ truncated(Normal(1.0, 0.5); lower = 0)
    # α ~ truncated(Normal(0, 2); lower = 0)
    β ~ truncated(Normal(4, 0.5); lower = 0)
    n = size(data.diag_occs[1][1], 1)
    C ~ MvNormal(fill(0.7,n), 0.5*I)
    # etot = func(γ, α, β, [C1, C2, C3, C4, C5], data, T)
    constant_shift ~ Uniform(0, 10)
    # constant_shift =0 
    etot =on_site_energy(Jh,β, C, constant_shift, data, T)
    # etot = func(γ, β, [C1, C2, C3], data, T)
    # noise
    σ ~ truncated(Normal(0.5, 0.5); lower = 0)
    return y ~ MvNormal(etot, σ * I)
end

@pooled_component Base.@kwdef struct TrainerSettings
    n_samples_per_training::Int
    n_points_per_training::Int
    sampler = NUTS()
end

@pooled_component struct Model
    α::Float64
    Jh::Float64
    C::Vector{Float64}
    constant_shift::Float64
    n_points::Int
    chain::Union{Nothing, Chains}
end

function train_model(m::Searcher)
    trainer_settings = m[TrainerSettings][1]
    data_train = prepare_data(m)
    @info "Training with $(length(data_train.energies)) data points"
    # chain = suppress() do
    chain=    sample(total_energy(data_train) | (;y=data_train.energies), trainer_settings.sampler, trainer_settings.n_samples_per_training)
    # end
    parameters = mean(chain)

    C = map(1:5) do i
        parameters[Symbol("C[$i]"), :mean]
    end
    # return Model(parameters[:Jh, :mean], parameters[:β, :mean], C, 0, length(data_train.energies), chain)
    return Model(parameters[:Jh, :mean], parameters[:β, :mean], C, parameters[:constant_shift, :mean], length(data_train.energies), chain)
end

struct ModelTrainer <: System end
function Overseer.requested_components(::ModelTrainer)
    return (TrainerSettings, Intersection, Model)
end

function Overseer.update(::ModelTrainer, m::AbstractLedger)
    trainer_settings = m[TrainerSettings][1]
    
    prev_model = isempty(m[Model]) ? nothing : m[Model][end]
    
    n_points = length(m[Unique]) + length(@entities_in(m,!Unique && Intersection && Results))
    prev_points = prev_model === nothing ? 0 : prev_model.n_points

    if n_points - prev_points > trainer_settings.n_points_per_training

        prev_chain = prev_model === nothing ? nothing : prev_model.chain

        # I'm not sure what we were trying to do actually works...
        model = train_model(m)
        # or max ?
        Entity(m, m[Template][entity(m[SearcherInfo],1)], model, Generation(length(m[Model].c.data)+1))
    end
end

function occ2eigvals_angles(occ)
    eigs, vecs = eigen(occ)
    return [clamp.(eigs.data, 0, 1); Angles(vecs[Up()]).θs; Angles(vecs[Down()]).θs]
end


# flat vector in + leading dimension (e.g. 5 for d orbitals)
# [dim up eigvals; dim down eigvals; (dim * (dim-1))/2 angles up; (dim * (dim-1))/2 angles down; etc next atoms;...]  
function eigvals_angles2occ(eigvals_angles::Vector, dim::Int)
    n_vals_per_at = 2 * dim + dim * (dim - 1)

    up_vals_r = 1:dim
    dn_vals_r = up_vals_r[end] + 1:2dim

    up_angles_r = dn_vals_r[end] + 1:dn_vals_r[end] + div(dim * (dim - 1), 2)
    dn_angles_r = up_angles_r[end] + 1:up_angles_r[end] + div(dim * (dim - 1), 2)
    
    return @views map(1:n_vals_per_at:length(eigvals_angles)) do start
    
        all_vals = eigvals_angles[start:start+n_vals_per_at-1]

        up_vals = all_vals[up_vals_r]
        dn_vals = all_vals[dn_vals_r]

        up_angles = Angles(all_vals[up_angles_r])
        dn_angles = Angles(all_vals[dn_angles_r])
        
        up_vecs = Matrix(up_angles)
        dn_vecs = Matrix(dn_angles)
        up_occ = Matrix(Eigen(up_vals, up_vecs))
        dn_occ = Matrix(Eigen(dn_vals, dn_vecs))
        occs = DFWannier.ColinMatrix(hcat(up_occ, dn_occ))
    end
end

function model_diag(model, func, U, eigvals_angles::Vector; use_penalty=false)

    occs = eigvals_angles2occ(eigvals_angles, 5) # for now we have one atom
    data = prepare_data([occs], U)
    
    @assert size(eigvals_angles, 1) == 30
    
    ys = func(model.α, model.Jh, model.C, model.constant_shift, data, Float64)[1]
    if use_penalty
        eigenvals = eigvals_angles[1:10]
        penalty =  1e2 * sum(@. (max(0, eigenvals-1) + max(0, -eigenvals)))
        ys += penalty
    end
    return ys
end

function optimize_cg(x0, model_diag, lower, upper)
    # conjugate gradient
    res = optimize(
        model_diag,
        lower, upper, x0,
        Fminbox(ConjugateGradient(; linesearch= LineSearches.BackTracking())),
        Optim.Options(iterations=100, outer_iterations=10, show_trace=true))
    x_min = minimizer(res);
    return eigvals_angles2occ(x_min, 5);
end
function optimize_nm(x0, model_diag, lower, upper)
    # conjugate gradient
    res = optimize(
        model_diag,
        lower, upper, x0,
        Fminbox(NelderMead()),
        Optim.Options(iterations=1000, outer_iterations=10, show_trace=false, show_every=100000))
    x_min = minimizer(res);
    return eigvals_angles2occ(x_min, 5);
end

function optimize_sa(x0, model_diag, lower, upper)
    # simulated annearling with bounds
    res = optimize(
        model_diag,
        lower, upper, x0,
        SAMIN(
            rt=0.95, # temperature reduction rate controls space covered, higher -> more coverage
            verbosity=1,
        ), 
        Optim.Options(iterations=10^6, show_trace=false, show_every=100000)
    )
    x_min = minimizer(res);
    return eigvals_angles2occ(x_min, 5);
end

function optim_limits(s::State)
    n = size(s.occupations[1], 1)
    lower = repeat(vcat(fill(-0.0, 2n), fill(-4π-1e-3, n*(n-1))), length(s.occupations))
    upper = repeat(vcat(fill(1.0, 2n), fill(4π +1e-3, n*(n-1))), length(s.occupations))
    return (; lower, upper)
end
    
function Optim.optimize(model, U, s::State, optimizer=Fminbox(NelderMead()), args...;
                        iterations = 1000,
                        outer_iterations = 5,
                        show_trace = false,
                        kwargs...)
    x = vcat(map(x->occ2eigvals_angles(x), s.occupations)...)
    lower, upper = optim_limits(s)

    res = optimize(x -> model_diag(model, on_site_energy, U, x), lower, upper, x, optimizer, Optim.Options(;iterations, outer_iterations, show_trace, kwargs...), args...)

    x_min = minimizer(res)
    return State(eigvals_angles2occ(x_min, size(s.occupations[1],1))), res
end
    
# @pooled_component Base.@kwdef struct OptimizeSettings
#     opts::Vector{Optimizer}
#     opt_options::Vector{Dict}
# end

struct ModelOptimizer <: System end
function Overseer.requested_components(::ModelOptimizer)
    return (TrainerSettings, Intersection, Model)
end
    
function Overseer.update(::ModelOptimizer, m::AbstractLedger)
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
    
    rand_search_comp = m[RandomSearcher]
    rand_search_e    = entity(rand_search_comp, 1)
    rand_search      = rand_search_comp[rand_search_e]
    base_e = entity(m[BaseCase], 1)
    base_state = m[Results][base_e].state
    template_str = m[Template][rand_search_e].structure
    nelec    = round.(Int, base_state.totoccs)
    norb     = size.(base_state.occupations, 1)
    

    # TODO multithreading
    while max_new > 0
        # random start
        x0 = rand_trial(norb, nelec).state
        
        s, res= optimize(model, U, x0)
        energy = res.minimum

        min_dist = Inf
        min_e = Entity(0)
        min_energy = Inf
        
        for e in @entities_in(m, Unique && Results)
            dist = Euclidean()(e.state, s)
            if dist < min_dist
                min_dist = dist
                min_e = e.e
                min_energy = 13.6 * (e.total_energy - m[Results][base_e].total_energy)
            end
        end

        for e in @entities_in(m, Trial)
            dist = Euclidean()(e.state, s)
            if dist < min_dist
                min_dist = dist
                min_e = e.e
            end
        end

        if min_dist > dist_thr && energy < min_energy / 2
            trial = Trial(s, RomeoDFT.ModelOptimized) # TODO other tag?

            # add new entity with optimized occ
            new_e = add_search_entity!(m, model_e, trial, m[Generation][model_e])
            m[Model][new_e] = model_e
            max_new -= 1
        end
    end

    # Use particle swarm to try and find the "Ground state"
    # out, res = optimize(model, 4.0, rand_trial(norb, nelec).state, ParticleSwarm(;lower=[fill(0, 10); fill(-4pi, 20)], upper=[fill(1, 10); fill(4pi, 20)], n_particles=20), iterations=2000)
    # if minimum(x->Euclidean()(out, x.state), @entities_in(m, Results)) > dist_thr && minimum(x->Euclidean()(out, x.state), @entities_in(m, Trial)) > dist_thr
    #     add_search_entity!(m, model_e, Trial(out, RomeoDFT.IntersectionMixed), m[Generation][model_e], Intersection(model_e, model_e))
    #     info.n_pending_calcs += 1
    # end

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

    
function test_model(name)
    l = load(Searcher(name))
    l[Entity(1)] = RomeoDFT.TrainerSettings(1000, 10, RomeoDFT.NUTS())
    l[Template][Entity(1)] = Entity(2)
    l[Entity(1)] = RomeoDFT.SearcherInfo(max_concurrent_trials=10)
    update(RomeoDFT.ModelTrainer(), l)
    get_metrics(l)
end


function model_evaluation(model, data)
    return RomeoDFT.on_site_energy(model.α, model.Jh, model.C, model.constant_shift, data, Float64)
end

function rmse(y, y_pred)
    return sqrt(mean((y_pred - y) .^ 2))
end
function get_metrics(l::Searcher)
    model = l[RomeoDFT.Model][end]
    data = RomeoDFT.prepare_data(l)
    loss_params = rmse(data.energies, model_evaluation(model, data))
    E_samples = RomeoDFT.predict(RomeoDFT.total_energy(data), model.chain);
    ys_mean = vec(mean(Array(RomeoDFT.group(E_samples, :y)); dims=1));
    ys_median = quantile(E_samples)[:, 3];
    loss_mean = rmse(data.energies, ys_mean)
    loss_median = rmse(data.energies, ys_median)
    std_parameters = sum(describe(model.chain)[1][:, :std])
    return (; loss_params, loss_mean, loss_median, std_parameters)
end

