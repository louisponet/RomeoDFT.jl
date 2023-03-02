struct Binner{F<:Function,T} <: System
    components::T
    binfunc::F
end

function Overseer.update(b::Binner, m::AbstractLedger)
    Overseer.ensure_component!(m, Bin)
    comps = map(c -> m[c], b.components)
    for e in @entities_in(comps[1])
        if length(comps) < 2 || all(x -> e in x, comps[2:end])
            binned = false
            for bin in m[Bin]
                if any(x -> b.binfunc(map(c -> c[x], comps)..., map(c -> c[e], comps)...),
                       bin.entities)
                    push!(bin.entities, e.e)
                    binned = true
                    break
                end
            end
            if !binned
                Entity(m, Bin([e]))
            end
        end
    end
end

function flatbands(all_out)
    first = all_out[:total_magnetization][end] > 0 ? :up : :down
    second = first == :up ? :down : :up
    bands = all_out[:bands]
    outbands = Vector{Float64}(undef,
                               2 * length(bands[first]) * length(bands[first][1].eigvals))
    icur = 1
    @inbounds for ib in 1:length(bands[:up])
        for s in (first, second)
            b = bands[s][ib]
            for e in b.eigvals
                outbands[icur] = e
                icur += 1
            end
        end
    end
    return outbands
end

function add_bands!(l)
    Overseer.ensure_component!(l, FlatBands)
    sj = l[SimJob]
    for e in @entities_in(l[SimJob])
        tp = joinpath(e.local_dir, "scf.out")
        if ispath(tp)
            l[FlatBands][e] = FlatBands(Float64[])
        end
    end
    @sync for e in @entities_in(l, SimJob && FlatBands)
        Threads.@spawn try
            o = DFC.FileIO.qe_parse_pw_output(joinpath(e.local_dir, "scf.out"))
            e.bands = flatbands(o)
        catch
            @warn "Something went wrong for $(e.e)"
        end
    end
end

function sssp_distance(bands1, bands2, fermi)
    function obj(Δ)
        n = 0
        d = 0.0
        for (b1, b2) in zip(bands1, bands2)
            if b1 <= fermi && b2 <= fermi
                d += (b1 - b2 + Δ)^2
                n += 1
            end
        end
        return sqrt(d / n)
    end
    return optimize(x -> obj(x[1]), [0.0]).minimum
end

"Calculates the order of eigenvectors of `occ1` which corresponds most to the order of eigenvectors in `occ2`."
function eigvec_order(occ1, occ2)
    eigvals1, eigvecs1 = eigen(occ1)
    eigvals2, eigvecs2 = eigen(occ2)
    dim = size(eigvecs1, 1)
    order = zeros(Int, dim)
    for c1 in 1:dim
        best = 0.0
        for c2 in 1:dim
            t = abs(dot(eigvecs1[:, c1], eigvecs2[:, c2]))
            if t > best && !(c2 in order)
                best = t
                order[c1] = c2
            end
        end
    end
    return (; eigvals1, eigvals2, eigvecs1, eigvecs2, order)
end

"""
write_xsf(filename::String, wfc::Wfc3D{T}) where T<:AbstractFloat
Writes the real part of the Wfc3D to a .xsf file that is readable by XCrysden or VESTA.
"""
function write_xsf(filename::String, l::AbstractLedger)
    gs = ground_state(l)
    structure = deepcopy(l[Template][gs].structure)
    if RelaxResults in l && gs in l[RelaxResults]
        Structures.update_geometry!(structure, l[RelaxResults][gs].final_structure)
    end
    open(filename,"w") do f
        write(f,"# Generated from PhD calculations\n")
        write(f, "CRYSTAL\n")
        c = Structures.ustrip.(structure.cell')
        write(f, "PRIMVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "CONVVEC\n")
        write(f, "$(c[1,1]) $(c[1,2]) $(c[1,3])\n")
        write(f, "$(c[2,1]) $(c[2,2]) $(c[2,3])\n")
        write(f, "$(c[3,1]) $(c[3,2]) $(c[3,3])\n")
        write(f, "PRIMCOORD\n")
        write(f, "$(length(structure.atoms)) 1\n")
        hub_at = 1
        for at in structure.atoms
            n = at.element.symbol
            p = Structures.ustrip.(at.position_cart)
            if at.dftu.U != 0
                write(f, "$n $(p[1]) $(p[2]) $(p[3]) 0 0 $(gs[Results].state.magmoms[hub_at])\n")
                hub_at += 1
            else
                write(f, "$n $(p[1]) $(p[2]) $(p[3])\n")
            end
        end
    end
end

function relative_energies(es; include_hub_energy = true, eV = true)
    E_conv_fac = eV ? 13.6056980659 : 1.0
    nat = length(first(es)[Template].structure.atoms)
    energies = include_hub_energy ? map(x -> x.total_energy * E_conv_fac / nat, es) :
               map(x -> dft_energy(x) * E_conv_fac / nat, es)
    return energies .- minimum(energies)
end

function relative_energies(l::AbstractLedger; kwargs...)
    relative_energies(filter(x -> x.converged, @entities_in(l[Results] && l[FlatBands] && l[Template] && !l[Simulation])); kwargs...)
end

function min_energy_per_generation(l::AbstractLedger)
    minima = zeros(maximum(x->x.generation, l[Generation], init=0))
    curmin = typemax(Float64)
    curgen = 0
    for e in @entities_in(l, Generation && Results)
        if e.converged
            if e.generation > curgen
                if curgen > 0
                    minima[curgen:e.generation-1] .= curmin
                end
                curgen = e.generation
            end
            if e.total_energy < curmin
                curmin = e.total_energy
            end
        end
    end
    if curgen <= length(minima)
        minima[curgen:end] .= curmin
    end
    return minima
end

function DFControl.bandgap(bands::FlatBands, fermi::Float64)
    maxval = maximum(x -> x <= fermi ? x : -1000.0, bands.bands)
    mincon = minimum(x -> x > fermi ? x : +1000.0, bands.bands)
    return mincon - maxval
end
function DFControl.bandgap(e::Overseer.EntityState)
    @assert FlatBands in e ArgumentError("No FlatBands in Entity")
    @assert Results in e  ArgumentError("No Results in Entity")
    return bandgap(e[FlatBands], e.fermi)
end
