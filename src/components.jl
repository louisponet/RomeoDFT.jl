include("Components/Components.jl")

function version_convert(::Type{T}, x) where {T}
    return error("No version conversion function for type $T from $(typeof(x))")
end

function version_convert(v, type_mods, id)
    cur_type = type_mods[id]
    if id == length(type_mods)
        return convert(cur_type, v)
    else
        return version_convert(convert(cur_type, v), type_mods, id + 1)
    end
end


#############################################
#  CONSTRUCTORS                             #
#############################################

function Error(e::Overseer.AbstractEntity, exc::Exception, trace::StackTraces.StackTrace)
    s = IOBuffer()
    println(s, "Error for $(Entity(e)):")
    showerror(s, exc, trace; backtrace=true)
    errormsg = String(resize!(s.data, s.size))
    @error errormsg 
    return Error(exc, join(trace, "\n"))
end

#############################################
#  SHOW                                     #
#############################################
function Base.show(io::IO, sim::Simulation)
    println(io, "Simulation:")
    for f in fieldnames(Simulation)
        print(io, "\t$f: ")
        if f == :template_structure
            print(io, "structure with $(length(sim.template_structure.atoms)) atoms")
        elseif f == :template_calculation
            print(io, "calculation with ")
            for flag in (:conv_thr, :mixing_beta, :Hubbard_mixing_beta, :Hubbard_conv_thr)
                if haskey(sim.template_calculation, flag)
                    print(io, "$flag=$(sim.template_calculation[flag]) ")
                end
            end
        else
            print(io, getfield(sim, f))
        end
        print(io, "\n")
    end
end

function Base.show(io::IO, res::Results)
    print(io, "Results(state: ", res.state)
    for f in fieldnames(Results)
        f == :state && continue
        print(io, ", $f: $(getfield(res, f))")
    end
end

Base.show(io::IO, f::FlatBands) = print(io, "FlatBands(nkpt: $(length(f.bands)))")

function Base.show(io::IO, err::Error)
    (print(io, "Error: "); showerror(io, err.err); print(io, "\nTrace:\n"); println(io,
                                                                                    err.stack))
end

#############################################
#  gencalc                                  #
#############################################
gencalc(job, settings::NSCFSettings) = 
    Calculations.gencalc_nscf(job[1], settings.kpoints)

gencalc(job, settings::ProjwfcSettings) = 
    Calculations.gencalc_projwfc(job[1], settings.Emin, settings.Emax, settings.deltaE)
    
function gencalc(job, wsettings::WannierSettings)
    for (elsym, projs) in wsettings.projections
        for a in job.structure[elsym]
            a.projections = projs
        end
    end
    return Calculations.gencalc_wan(job, wsettings.dos_ratio; plot_wannier = wsettings.plot_wannier)
end

function gencalc(job, settings::BandsSettings)
    if settings.kpoints isa Int
        kpath = filter(x -> all(y -> y in (1 / 3, 0.5, 0.0), x[1:3]),
                       Structures.high_symmetry_kpath(job.structure,
                                                      settings.kpoints))
        kpath[end] = (kpath[end][1:3]..., 1.0)
        return Calculations.gencalc_bands(job["scf"], kpath)
    else
        return Calculations.gencalc_bands(job["scf"], kpoints)
    end
end

function add_calc!(job, settings)
    calc = gencalc(job, settings)
    calc.run = true
    cid = findfirst(x -> x.name == calc.name, job.calculations)
    if cid !== nothing
        job.calculations[cid] = calc
    else
        push!(job, calc)
    end
end

dft_energy(res::Results) = res.total_energy - res.Hubbard_energy

function trypop!(c::Overseer.AbstractComponent, e::AbstractEntity)
    if e in c
        pop!(c, e)
    end
end
