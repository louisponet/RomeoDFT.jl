# Should be ran after ResultsProcessor so as to not remove the SimJob before this system runs
struct ElectrideCreator <: System end

Overseer.requested_components(::ElectrideCreator) = (ElectrideSettings, PPSettings, Results)

function Overseer.update(::ElectrideCreator, l::AbstractLedger)

    @error_capturing for e in @entities_in(l, ElectrideSettings && Results && !PPSettings && !ShouldRerun && SimJob)
        
        j = local_load(Job(joinpath(l, e)))

        o = outputdata(j)["scf"]
        bands = o[:bands]

        metallic = bandgap(o[:bands], o[:n_electrons]) != 0.0

        if metallic
            emin = e.fermi - e.range_below
            emax = e.fermi + e.range_above
        else
            homo = maximum(x -> maximum(x.eigvals) > e.fermi ? -Inf : maximum(x.eigvals), Iterators.flatten(bands))
            lumo = minimum(x -> minimum(x.eigvals) < e.fermi ? Inf : minimum(x.eigvals), Iterators.flatten(bands))
            
            emin = homo - e.range_below
            emax = homo + e.range_above
        end

        l[e] = PPSettings(Dict(:emax           => emax,
                               :emin           => emin,
                               :filplot        => "density.filplot",
                               :plot_num       => 23,
                               :spin_component => 3,
                               :fileout        => "density.fileout",
                               :iflag          => 3,
                               :output_format  => 6))
        should_rerun(l, e)
    end
end


struct ElectrideProcessor <: System end

Overseer.requested_components(::ElectrideProcessor) = (ElectrideSettings,  Pulled)

function Overseer.update(::ElectrideProcessor, l::AbstractLedger)
    
    @error_capturing for e in @entities_in(l, ElectrideSettings && Pulled && SimJob)
        if e.job.calculations[end].name == "pp"
            l[e] = Done(false)
        end
    end
end
