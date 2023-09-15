using MCMCChains, Turing, Distributions, StatsBase, Combinatorics
using Optim
using Optim: minimizer, LineSearches 

function E_hund_exchange(diag_occs, J)
    # -J_H Σ (1/2 + S_α * S_β)
    up = map(d ->d[:, 1], diag_occs)
    down = map(d ->d[:, 2], diag_occs)
    s = 0.
    for i = eachindex(up)
        for (p,q) = combinations(1:size(diag_occs, 1), 2)
            s += up[i][p] * up[i][q] + down[i][p] * down[i][q]
        end
    end
    return -J * s
end

function E_U(n, U)
    s = 0.0
    for m1 = 1:size(n, 1)
        s += n[m1,m1]
        for m2 = 1:size(n, 2)
            s -= n[m1, m2] *  n[m2, m1]
        end
    end
    return 0.5 * U * s
end

function on_site_energy(γ, β, C, constant_shift, data, T)    
    magmoms_squared_sum = data.magmoms_squared_sum
    Eu = data.E_Us
    diag_occs = data.diag_occs
    diag_occs_sum = data.diag_occs_sum
    diag_occs_squared_sum = data.diag_occs_squared_sum

    natoms = Int(length(diag_occs_squared_sum[1]) / length(C))

    etot = zeros(T, length(magmoms_squared_sum))
    @inbounds for i = 1:length(magmoms_squared_sum)
        m = magmoms_squared_sum[i]
        diag_occ_square = diag_occs_squared_sum[i]
        diag_occ_sum = diag_occs_sum[i]
        diag_occ = diag_occs[i]
        Eu_ = Eu[i]  # +U energy
        Eh = E_hund_exchange(diag_occ, γ)         # Hund's exchange
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
    
    conv_res = @entities_in(l, Results && Unique)
    
    states = map(x->x.state, conv_res)
    energies = map(x->x.total_energy - base_energy, conv_res)
    energies .*= 13.6
    occupations = map(x->x.occupations, states)
    return (; energies, prepare_data(occupations, U)...)
end

# TODO duplicate code
function prepare_data(occs, U)
    magmoms = map(occs) do occupation
        map(occupation) do occ
            tr(view(occ, Up())) - tr(view(occ, Down())) 
        end
    end
    magmoms_squared_sum   = map(magmom -> sum(m -> m^2, magmom), magmoms)
    E_Us                  = map(occ -> sum(m -> E_U(view(m, Up()), U) + E_U(view(m, Down()), U), occ), occs)
    diag_occs             = map(occ -> map(x -> hcat(diag(view(x, Up())), diag(view(x, Down()))), occ), occs)
    diag_occs_sum         = map(occ -> sum(x -> (diag(view(x, Up())) .+ diag(view(x, Down()))), occ), occs)
    diag_occs_squared_sum = [d.^2 for d in diag_occs_sum]
    return (; magmoms_squared_sum, diag_occs_sum, E_Us, diag_occs_squared_sum, diag_occs)
end

@model function total_energy(data, ::Type{T}=Float64) where {T}
    # physics
    α ~ truncated(Normal(0, 2); lower = 0)
    # α ~ truncated(Normal(0, 2); lower = 0)
    Jh ~ truncated(Normal(0, 2); lower = 0)
    n = size(data.diag_occs[1][1], 1)
    C ~ MvNormal(zeros(T, n), fill(T(2), n))
    # etot = func(γ, α, β, [C1, C2, C3, C4, C5], data, T)
    constant_shift ~ Uniform(-10, 10)
    etot =on_site_energy(α, Jh, C, constant_shift, data, T)
    # etot = func(γ, β, [C1, C2, C3], data, T)
    # noise
    σ ~ truncated(Normal(0, 2); lower = 0)
    return y ~ MvNormal(etot, σ * I)
end

@pooled_component Base.@kwdef struct TrainerSettings
    n_samples_per_training::Int
    n_points_per_training::Int
    sampler = NUTS()
end

@component struct Model
    α::Float64
    Jh::Float64
    C::Vector{Float64}
    constant_shift::Float64
    n_points::Int
    chain::Union{Nothing, Chains}
end

struct ModelTrainer <: System end
function Overseer.requested_components(::ModelTrainer)
    return (TrainerSettings, Intersection, Model)
end

function Overseer.update(::ModelTrainer, m::AbstractLedger)
    trainer_settings = m[TrainerSettings][1]
    
    prev_model = isempty(m[Model]) ? nothing : m[Model][end]
    
    n_points = length(m[Unique]) + length(m[Intersection])
    prev_points = prev_model === nothing ? 0 : prev_model.n_points

    if n_points - prev_points > trainer_settings.n_points_per_training
        data_train = prepare_data(m)

        prev_chain = prev_model === nothing ? nothing : prev_model.chain
        
        if prev_chain !== nothing
    
            chain = sample(prev_model.chain, total_energy(data_train) | (;y=data_train.energies), trainer_settings.sampler, trainer_settings.n_samples_per_training)
        else
            chain = sample(total_energy(data_train) | (;y=data_train.energies), trainer_settings.sampler, trainer_settings.n_samples_per_training)
        end
        
        parameters = mean(chain)
        display(parameters)

        C = map(1:5) do i
            parameters[Symbol("C[$i]"), :mean]
        end
        # or max ?
        Entity(m, Model(parameters[:α, :mean], parameters[:Jh, :mean], C, parameters[:constant_shift, :mean], n_points, chain))
    end
end

function occ2x(occ)
    eigen_vals = Float64[]
    angles = Float64[]
    up=view(occ, Up())
    down = view(occ, Down())
    for m in [up, down]
        vals, vecs = eigen(m)
        agl = Angles(vecs) 
        append!(eigen_vals, vals)
        append!(angles, agl.flat)
    end
    return vcat(eigen_vals, angles)
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

        up_angles = Angles(all_vals[up_angles_r], 1.0,(5,5))
        dn_angles = Angles(all_vals[dn_angles_r], 1.0, (5,5))
        
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

function optimize_sa(x0, model_diag, lower, upper)
    # simulated annearling with bounds
    res = optimize(
        model_diag,
        lower, upper, x0,
        SAMIN(
            rt=0.95, # temperature reduction rate controls space covered, higher -> more coverage
            verbosity=1,
        ), 
        Optim.Options(iterations=10^6, show_trace=true)
    )
    x_min = minimizer(res);
    return eigvals_angles2occ(x_min);
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
    
    U = getfirst(x->x.dftu.U != 0, m[Template][1].structure.atoms).dftu.U

    # check how many calcs are pending on server
    # if free room -> generate some more calcs with the current model
    info = m[SearcherInfo][1]
    max_new = max(0, info.max_concurrent_trials - (info.n_running_calcs + info.n_pending_calcs))
    max_new <= 0 && return 
 
    # find previous unique, intersection, other trials about to run
    # to check whether it's duplication
    unique_states = @entities_in(m, Unique && Results)
    pending_states = @entities_in(m, SimJob && Submit && Trial)

    # Generate Trial with origin ModelOptimized + increment generation
    # Severely limit the amount of new intersections, maybe only 6 per generation or whatever
    #  new trial -> hopefully new unique -> 3 x whatever prev unique intersections -> train model whenever Y new unique states are found from those -> rinse repeat 
    lower = vcat(ones(10) * (-0.0), ones(20) * (-π-1e-3))
    upper = vcat(ones(10) * 1.0, ones(20) * (π +1e-3))

    str = m[Template][1].structure
    natoms = length(filter(ismagnetic, str.atoms))
    nshell = size(m[Results][1].state.occupations[1], 1)
    dist_thr = 1e-2
    curgen = maximum_generation(m)
    # TODO multithreading
    while max_new > 0
        # random start
        x0 = rand(30) # TODO slightly better random?
        occ = optimize_cg(x0, x -> model_diag(model, on_site_energy, U, x), lower, upper)
        state0 = State(occ)
        distances = map(e->Euclidean()(e[Results].state, state0), unique_states) ∪
                    map(e->Euclidean()(e[Trial].state, state0), pending_states)
                    
        if minimum(distances; init=Inf) < dist_thr
            continue
        end
        trial = Trial(state0, RomeoDFT.ModelOptimized) # TODO other tag?
        gen = Generation(curgen + 1) # TODO this is wrong

        # add new entity with optimized occ
        new_e = Entity(m, trial, gen) 
        max_new -= 1
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

    
