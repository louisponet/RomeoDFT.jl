"""
    State

Represents the local state of a system. This is given by the occupation matrices of the local orbitals.
These are usually the valence shells for example to which the +U correction is applied in DFT + U calculations.
Can be initialized using the `hubbard` entry in the outputdata of an scf calculation. Using [`generate_Hubbard_occupations`](@ref) this can generate the `:Hubbard_occupations` scf input parameter which will be
used as the target during a constrained scf calculation. [`generate_starting_ns_eigenvalue`](@ref)
generates the `:starting_ns_eigenvalue` parameter instead.
"""
mutable struct State{T<:AbstractFloat,VT,MT}
    occupations::Vector{MT}
    magmoms::Vector{T}
    eigvals::Vector{VT}
    eigvecs::Vector{MT}
    totoccs::Vector{T}
    angles::Vector{Angles{T,2,VT}} # These are the angles w.r.t. eigenvectors
end

State() = State(ColinMatrixType[], Float64[], MagneticVectorType[], ColinMatrixType[], Float64[], AnglesType[])

State(hubs::Vector{<:NamedTuple}) = State(map(x -> DFW.ColinMatrix(x.occupations.up, x.occupations.down), hubs))

function State(occs::Vector{<:AbstractMatrix})
    magmoms = [tr(x[Up()]) - tr(x[Down()]) for x in occs]
    return State(occs, magmoms, MagneticVectorType[], ColinMatrixType[], tr.(occs), AnglesType[])
end

function State(occmat::Array{T, 4} where T)
    occs = map(1:size(occmat,4)) do i
        nl = findfirst(iszero, diag(view(occmat,:,:, 1, i))) - 1
        DFWannier.ColinMatrix(view(occmat,1:nl, 1:nl, 1, i), view(occmat,1:nl,1:nl,2,i))
    end
    return State(filter(x->sum(x) != 0, occs))
end

Base.convert(::Type{T}, s)    where T <: State = State(s.occupations)
Base.convert(::Type{T}, s::T) where T <: State = State(s.occupations)

@inline Base.length(s::State) = length(s.occupations)

@inline function (dist::Distances.Euclidean)(s1::State, s2::State)
    return sum(x -> evaluate(dist, x[1], x[2]), zip(s1.occupations, s2.occupations))
end

@inline function Base.getproperty(s::State, sym::Symbol)
    if sym == :eigvals
        t = getfield(s, :eigvals)
        !isempty(t) && return t
        eigs = eigen.(s.occupations)
        s.eigvecs = [x.vectors for x in eigs]
        return s.eigvals = [x.values for x in eigs]
    elseif sym == :eigvecs
        t = getfield(s, :eigvecs)
        !isempty(t) && return t
        eigs = eigen.(s.occupations)
        s.eigvals = [x.values for x in eigs]
        return s.eigvecs = [x.vectors for x in eigs]
    elseif sym == :angles
        t = getfield(s, :angles)
        !isempty(t) && return t
        angles = [Angles(x) for x in s.eigvecs]
    else
        return getfield(s, sym)
    end
end

Base.show(io::IO, s::State) = print(io, "State(nions: $(length(s.magmoms)))")

function Base.show(io::IO, ::MIME"text/plain", s::State)
    println(io, "State:")
    print(io, "Magmoms:")
    # print(io, "\t")
    for m in s.magmoms
        print(io, " $(round(m, digits=3))")
    end
    print(io, "\n")
    println(io, "RomeoDFT:")
    for (n, u) in zip(("Up", "Down"), (Up(), Down()))
        print(io, "\t$n:")
        for e in s.occupations
            print(io, " $(round(tr(e[u]), digits=3))")
        end
        print(io, "\n")
    end
    print(io, "\tTotal:")
    for e in s.totoccs
        print(io, " $(round(e, digits=3))")
    end
    return println(io)
end

# State algebra
Base.:(+)(s1::State, s2::State) = State([o1 + o2 for (o1, o2) in zip(s1.occupations, s2.occupations)])  
Base.:(-)(s1::State, s2::State) = State([o1 - o2 for (o1, o2) in zip(s1.occupations, s2.occupations)])  
Base.:(*)(f::Number, s1::State) = State([f * o for o in s1.occupations])  
Base.:(*)(s1::State, f::Number) = State([f * o for o in s1.occupations])  
Base.:(/)(s1::State, f::Number) = State([f / o for o in s1.occupations])  


## EulerAngles interface
function EulerAngles.Angles(m::DFWannier.ColinMatrix)
    up = Angles(m[Up()])
    down = Angles(m[Down()])
    return Angles(DFWannier.MagneticVector([up.θs; down.θs]), 1.0, size(m))
end

function EulerAngles.Angles(m::DFWannier.NonColinMatrix)
    t = Angles(m.data)
    return Angles(DFWannier.MagneticVector(t.θs), 1.0, size(m))
end

function Base.Matrix(a::Angles{T,2,<:DFWannier.MagneticVector} where {T})
    if a.n[1] != a.n[2]
        up = Matrix(Angles(a.θs[Up()], 1.0, (a.n[1], div(a.n[2], 2))))
        down = Matrix(Angles(a.θs[Down()], 1.0, (a.n[1], div.(a.n[2], 2))))
        return DFWannier.ColinMatrix(up, down)
    else
        return DFWannier.NonColinMatrix(Matrix(Angles(a.θs.data, a.r, a.n)))
    end
end

EulerAngles.Angles(s::State) = map(Angles, s.occupations)

struct EulerDist <: Distances.PreMetric end

@inline function (dist::EulerDist)(s1::State, s2::State)
    return sum(x -> evaluate(dist, x[1], x[2]), zip(s1.angles, s2.angles))
end

@inline (dist::EulerDist)(a1::Angles, a2::Angles) = evaluate(dist, a1, a2)

function Distances.evaluate(::EulerDist, a1::Angles, a2::Angles)
    s = 0.0
    for (θ1, θ2) in zip(a1.θs, a2.θs)
        co = cos(θ1) * cos(θ2)
        si = sin(θ1) * sin(θ2)
        s += co + si
    end
    return 1 - s / length(a1.θs)
end

distance_to_target(states::Vector{State}) = tmap(x -> Euclidean()(x, target), states)

function distance_to_target(states)
    target = State(states[1])
    return tmap(x -> Euclidean()(State(x), target), states)
end

"""
    generate_Hubbard_occupations(state)

Generates the `:Hubbard_occupations` entry for a constrained scf calculation.
"""
function generate_Hubbard_occupations(hubstate::State)
    nat = length(hubstate)
    occs = hubstate.occupations
    Hubbard_occupations = zeros(7, 7, 4, nat)
    for (ia, occ) in enumerate(occs)
        for (is, s) in enumerate((Up(), Down()))
            o = occ[s]
            dim = size(o, 1)
            for m1 in 1:dim, m2 in 1:dim
                t = o[m1, m2]
                Hubbard_occupations[m1, m2, is, ia] = t == 0.0 ? 1e-5 : t
            end
        end
    end
    return Hubbard_occupations
end
function generate_Hubbard_occupations(hubstate::Vector{<:NamedTuple})
    return generate_Hubbard_occupations(State(hubstate))
end
function generate_Hubbard_occupations(s::State, str::Structure)
    occ = generate_Hubbard_occupations(s)
    nat = length(str.atoms)
    full_occ = zeros(7, 7, 4, nat)
    hub_at = 1
    for i = 1:nat
        if str.atoms[i].dftu.U != 0.0
            full_occ[:, :, :, i] .= occ[:, :, :, hub_at]
            hub_at += 1
        end
    end
    return full_occ
end

"""
    starting_ns_eigenvalue(state)

Generates the `:starting_ns_eigenvalue` entry for a constrained scf calculation.
"""
function generate_starting_ns_eigenvalue(hubstate::State)
    nat = length(hubstate)
    occs = hubstate.eigvals
    starting_ns = zeros(7, 4, nat)
    dim = length(occs[1][Up()])
    for (ia, occ) in enumerate(occs)
        for (is, s) in enumerate((Up(), Down()))
            for m1 in 1:dim
                t = occ[s][m1]
                starting_ns[m1, is, ia] = t == 0.0 ? 1e-5 : t
            end
        end
    end
    return starting_ns
end
function generate_starting_ns_eigenvalue(hubstate::Vector{<:NamedTuple})
    return generate_starting_ns_eigenvalue(State(hubstate))
end
