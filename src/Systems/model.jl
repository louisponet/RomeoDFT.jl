using MCMCChains, Turing, Distributions, StatsBase, Combinatorics

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
    diag_occupations = data.diag_occs
    diag_occupations_sum = data.diag_occupations_sum
    diag_occupations_squared_sum = data.diag_occupations_squared_sum

    natoms = Int(length(diag_occupations_squared_sum[1]) / length(C))

    etot = zeros(T, length(magmoms_squared_sum))
    @inbounds for i = 1:length(magmoms_squared_sum)
        m = magmoms_squared_sum[i]
        diag_occ_square = diag_occupations_squared_sum[i]
        diag_occ_sum = diag_occupations_sum[i]
        diag_occ = diag_occupations[i]
        Eu_ = Eu[i]  # +U energy
        Eh = E_hund_exchange(diag_occ, γ)         # Hund's exchange
        
        Ecf = sum(y->y[1] * y[2], zip(diag_occ_square, repeat(C, natoms)))
        
        E_colombic = - β * sum(diag_occ_sum)       # e-N interactions
        etot[i] = Eh + Ecf + Eu_ + E_colombic
    end
    return etot + constant_shift * length(C)
end

function prepare_data(l::Searcher)

    base_energy = l[Results][entity(l[BaseCase],1)].total_energy
    
    nmagats = length(filter(ismagnetic, l[Template][1].structure.atoms))
    
    U = DFC.getfirst(x->x.dftu.U != 0, l[Template][1].structure.atoms).dftu.U
    
    conv_res = @entities_in(l, Results && Unique)
    
    states = map(x->x.state, conv_res)
    energies = map(x->x.total_energy - base_energy, conv_res)
    energies .*= 13.6
    
    occupations = map(x->x.occupations, states)
    E_Us = map(occ-> sum(m -> E_U(view(m, Up()), U) + E_U(view(m, Down()), U), occ), occupations)
    magmoms_squared_sum = map(x->sum(y-> y^2, x.magmoms), states)
    diag_occs = map(occs->map(x->hcat(diag(view(x, Up())), diag(view(x, Down()))), occs), occupations)
    diag_occupations_sum = map(occ -> sum(x-> (diag(view(x, Up())) .+ diag(view(x, Down()))), occ), occupations)
    diag_occupations_squared_sum = [d.^2 for d in diag_occupations_sum]
    return (; magmoms_squared_sum, diag_occs, E_Us, energies, states, diag_occupations_sum, diag_occupations_squared_sum)
end

@model function total_energy(data, ::Type{T}=Float64) where {T}
    # physics
    α ~ truncated(Normal(0, 2); lower = 0)
    # α ~ truncated(Normal(0, 2); lower = 0)
    Jh ~ truncated(Normal(0, 2); lower = 0)
    n = size(data.diag_occs[1], 1)
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

function Overseer.update(::ModelTrainer, m::AbstractLedger)

    trainer_settings = m[TrainerSettings][1]
    
    prev_model = m[Model][end]
    
    n_points = length(m[Unique]) + length(m[Intersection])
    
    if n_points - prev_model.n_points > trainer_settings.n_points_per_training
        data_train = prepare_data(m)

        prev_chain = prev_model.chain
        if prev_chain !== nothing
    
            chain = sample(prev_model.chain, total_energy(data_train) | (;y=data_train.energies), trainer_settings.sampler, trainer_settings.n_samples_per_training)
        else
            chain = sample(total_energy(data_train) | (;y=data_train.energies), trainer_settings.sampler, trainer_settings.n_samples_per_training)
        end
        
        parameters = mean(chain)

        # or max ?
        Entity(m, Model(parameters[:α, :mean], parameters[:Jh, :mean], parameters[:C, :mean], parameters[:constant_shift, :mean], chain))
    end
end

struct ModelOptimizer <: System end
    
function Overseer.update(::ModelOptimizer, m::AbstractLedger)
    # check how many calcs are pending
    # if free room -> generate some more calcs with the current model
    # check Unique + Results and Trial whether it's duplication on either
    # Generate Trial with origin ModelOptimized + increment generation

    # Severely limit the amount of new intersections, maybe only 6 per generation or whatever
    #  new trial -> hopefully new unique -> 3 x whatever prev unique intersections -> train model whenever Y new unique states are found from those -> rinse repeat 
end

