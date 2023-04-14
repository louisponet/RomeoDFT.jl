"""
write_xsf(filename::String, wfc::Wfc3D{T}) where T<:AbstractFloat
Writes the real part of the Wfc3D to a .xsf file that is readable by XCrysden or VESTA.
"""
function write_xsf(filename::String, gs::Overseer.EntityState)
    structure = deepcopy(gs[Template].structure)
    if RelaxResults in gs
        Structures.update_geometry!(structure, gs[RelaxResults].final_structure)
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

function write_xsf(filename::String, l::AbstractLedger)
    gs = ground_state(l)
    write_xsf(filename, gs)
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

function ground_state(es; by=x->x.total_energy)
    min_entity = Entity(0)
    min_energy = typemax(Float64)
    for e in es
        if e.converged
            test = by(e)
            if test < min_energy
                min_energy = test
                min_entity = e.e
            end
        end
    end
    return min_entity
end

function ground_state(l::AbstractLedger; kwargs...)
    return l[ground_state(@entities_in(l, Results && FlatBands); kwargs...)]
end

function unique_states(es; thr = 1e-4, bands = true)
    iu = ones(length(es))
    @inbounds for i in 1:length(es)-1
        e1 = es[i]
        if iu[i] == 1
            Threads.@threads for j in i+1:length(es)
                if iu[j] == 1
                    e2 = es[j]
                    dist = bands ? sssp_distance(e1.bands, e2.bands, e1.fermi) :
                           Euclidean()(e1.state, e2.state)
                    if dist < thr
                        iu[j] = 0
                    end
                end
            end
        end
    end
    return es[findall(isequal(1), iu)]
end

function unique_states(l::AbstractLedger; kwargs...)
    return unique_states(filter(x -> x.converged,
                                  @entities_in(l[Results] && l[FlatBands])); kwargs...)
end

