using .Plots


function Overseer.update(::BandsPlotter, m::AbstractLedger)
    processed = 0
    simn = simname(m)
    @debugv 2 "$simn - [START] - BandsPlotter"
    @sync for e in @entities_in(m, SimJob && BandsSettings && TimingInfo && !BandsResults && !Error)
        if processed == 10
            break
        end
        curt = now()
        bp = joinpath(e.local_dir, "bands.out")
        if ispath(bp) && filesize(bp) != 0
            o = Dict()
            @tspawnat 1 try
                plot_path = joinpath(e.local_dir, "bands.png")
                p = plot()
                if e in m[ProjwfcSettings]
                    if ispath(joinpath(e.local_dir, "projwfc.out"))
                        psetting = m[ProjwfcSettings][e]
                        suppress() do
                            return p = plot(local_load(Job(e.local_dir)), e.ymin, e.ymax,
                                            psetting.dos_ratio)
                        end
                    else
                        return
                    end
                else
                    suppress() do
                        return p = plot(local_load(Job(e.local_dir)), bsetting.ymin,
                                        bsetting.ymax; outdat = o)
                    end
                end
                suppress() do
                    return savefig(p, plot_path)
                end
                processed += 1
                m[e] = BandsResults(plot_path)
            catch err
                m[e] = Error(err, stacktrace(catch_backtrace()))
            end
        end
        e.postprocessing += Dates.datetime2unix(now()) -
                                           Dates.datetime2unix(curt)
    end
    @debugv 2 "$simn - [STOP] - BandsPlotter"
end
function plot_states(tl;
                     include_hub_energy = true,
                     rel_energy = true,
                     fit = false,
                     y_property = "energies",
                     x_property = "magmoms",
                     z_property = "band_distances",
                     p = plot(),
                     kwargs...)
                     
    nat = isempty(tl[Template]) ? length(tl[Simulation][1].template_structure.atoms) : length(tl[Template][1].structure.atoms)
    es = filter(x -> x.converged, collect(@entities_in(tl[Results] && tl[FlatBands] && !tl[Simulation] && !tl[BaseCase])))
    energies = include_hub_energy ? map(x -> x.total_energy * E_conv_fac / nat, es) :
               map(x -> dft_energy(x) * E_conv_fac / nat, es)
    properties = Dict()

    properties["energies"] = rel_energy ? energies .- minimum(energies) : energies
    properties["abs_magmoms"] = map(x -> sum(y -> abs(y), x.state.magmoms), es)
    properties["magmoms"] = map(x -> sum(x.state.magmoms), es)

    m, mid = findmin(energies)
    occs = map(x -> x.state, es)
    properties["occupation_distances"] = map(occs) do x
        tot = 0.0
        for (o1, o2) in zip(x.occupations, occs[mid].occupations)
            d1 = Euclidean()(o1, o2)
            d2 = Euclidean()(o1, DFWannier.ColinMatrix(o2[Down()], o2[Up()]))
            tot += min(d1, d2)
        end
        return tot
    end
    minbands = tl[FlatBands][es[mid]].bands
    fermi = es[mid].fermi
    properties["band_distances"] = map(x -> sssp_distance(x.bands, minbands, fermi), es)

    labels = Dict("band_distances" => "\n"*"eta(n, n_GS)",
                  "occupation_distances" => "d_o(n, n_GS)",
                  "energies" => "(E - E_GS)/N (eV)",
                  "abs_magmoms" => nat > 2 ? "sum(|m_I|) (mu_b)" :
                                   "|m| (mu_b)",
                  "magmoms" => nat > 2 ? "sum(m_I) (mu_b)" :
                               "m (mu_b)")
    base_id = findfirst(x->x in tl[BaseCase], es)
    if z_property == "band_distances"
        scatter!(p, properties[x_property], properties[y_property];
                    zcolor = properties[z_property],
                    xguide = labels[x_property],
                    yguide = labels[y_property],
                    colorbar_title = labels[z_property], label = "",
                    size = (600, 300), kwargs...)
    elseif z_property == "none"
        scatter!(p,properties[x_property], properties[y_property];
                    xguide = labels[x_property],
                    yguide = labels[y_property],
                    label = "",
                    size = (600, 300), kwargs...)
    else
        scatter!(p, properties[x_property], properties[y_property];
                    zcolor = properties[z_property],
                    xguide = labels[x_property],
                    yguide = labels[y_property],
                    colorbar_title = labels[z_property], label = "",
                    size = (600, 300), kwargs...)
    end
    if base_id !== nothing
        scatter!(p, [properties[x_property][base_id]], [properties[y_property][base_id]], color = :red, label="vanilla QE")
    end
    return p
end

export plot_states
