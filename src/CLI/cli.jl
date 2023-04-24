using Comonicon

"""
Orchestrator commands
"""
module Orchestrator
    using ..RomeoDFT: start_orchestrator, client_stop_orchestrator, orchestrator_eval
    using ..RomeoDFT
    using Comonicon
    
    """
    start orchestrator.
    """
    @cast function start()
        return start_orchestrator()
    end
    
    """
    stop orchestrator.
    """
    @cast function stop()
        return client_stop_orchestrator()
    end
    """
    stop orchestrator.
    """
    @cast function status()
        return print(orchestrator_eval("status()"))
    end
    """
    pause orchestrator.
    """
    @cast function pause()
        return print(orchestrator_eval("pause_orchestrator()"))
    end
    """
    resume orchestrator.
    """
    @cast function resume()
        return print(orchestrator_eval("resume_orchestrator()"))
    end
end
@cast Orchestrator

"""
Searcher commands
"""
module searcher
    using ..RomeoDFT: searchers_dir, searcher_name, orchestrator_eval, setup_search, orchestrator_submit
    using ..RomeoDFT.RemoteHPC: InteractiveUtils
    using Comonicon
    """
    load a searcher.

    # Args
    
    - `dir`: directory where to find searcher
    """
    @cast function load(dir::String)
        if !isabspath(dir)
            dir = searchers_dir(dir)
        end
        if !ispath(dir)
            error("No such file or directory")
        end
        return print(orchestrator_eval("load_searcher(\"$dir\")"))
    end
    
    """
    unload a searcher.

    # Args
    
    - `dir`: directory where to find searcher
    """
    @cast function unload(dir::String)
        return print(orchestrator_eval("unload_searcher(\"$(searcher_name(dir))\")"))
    end

    """
    start a searcher.
    
    # Args

    - `name`: the name or directory of the searcher.

    # Options

    - `--verbosity, -v`: logging verbosity
    - `--sleep-time, -s`: sleep time between searcher updates
    """
    @cast function start(name::String; verbosity::Int = -1, sleep_time::Float64 = 30.0)
        return print(orchestrator_eval("start_searcher(\"$name\"; verbosity=$verbosity, sleep_time=$sleep_time)"))
    end
    
    """
    stop a searcher.
    
    # Args

    - `name`: the name or directory of the searcher.
    """
    @cast function stop(name::String)
        @info "Stopping searcher and removing pending jobs from queue..."
        return orchestrator_eval("stop_searcher(\"$name\")")
    end
    
    """
    list loaded searchers.
    """
    @cast function list()
        return print(orchestrator_eval("list_searchers()"))
    end
    
    """
    opens the log of a searcher in your favorite editor.
    
    # Args

    - `name`: the name or directory of the searcher.
    """
    @cast function log(name::String)
        InteractiveUtils.edit(joinpath(searchers_dir(name), "log.log"))
    end
    
    """
    create a new searcher.

    # Args

    - `name`: name or directory of the searcher
    - `scf_file`: template QE scf file for input parameters

    # Options
    
    - `--structure-file, -s`: CIF or QE input where to take the structure from. Defalts to `scf_file`
    - `--primitive`: set this to true to first try to find the primitive unit cell
    - `--supercell-a=<1>`: number specifying the number of unit cells along `a` direction
    - `--supercell-b=<1>`: number specifying the number of unit cells along `b` direction
    - `--supercell-c=<1>`: number specifying the number of unit cells along `c` direction
    - `--verbosity=<0>`: the logging verbosity, higher = more
    - `--sleep-time=<30>`: time in seconds between updates of the searcher
    - `--nrand=<10>`: how many random trials should be performed in a random search generation
    - `--unique-thr=<0.1>`: threshold that determines the uniqueness of electronic states (uses `sssp_distance`)
    - `--mindist-ratio=<0.25>`: minimum distance a trial should have as the ratio of the mean distance between current unique states
    - `--stopping-unique-ratio=<0.2>`: the ratio of unique states to trials below which the searching will stop
    - `--stopping-n-generations=<3>`: the search will stop if for this amount of successive generations the unique to trial ratio is lower than the `stopping-unique-ratio`
    - `--Hubbard-maxstep=<100>`: maximum constraining steps
    - `--Hubbard-mixing-beta=<0.4>`: mixing used to update the constraints during scf iterations
    - `--Hubbard-strength=<1.0>`: strength of the constraining potential
    - `--Hubbard-conv-thr=<0.1>`: threshold euclidean distance per atom between trial and current occupation matrices after which the constraints are released
    - `--electron-maxstep=<500>`: see QE documentation
    - `--relax-unique=<false>`: whether a structural relaxation should be ran on each unique state
    - `--relax-base=<false>`: whether a relaxation should be ran so the search happens for the relaxed structure
    - `--relax-force-convergence-threshold=<1e-3>`: `forc_conv_thr` in QE
    - `--relax-energy-convergence-threshold=<1e-4>`: `energy_conv_thr` in QE
    - `--relax-ion-dynamics=<"bfgs">`: `ion_dynamics` in QE
    - `--relax-cell-dynamics=<"bfgs">`: `cell_dynamics` in QE
    - `--relax-no-symmetry=<false>`: whether symmetry should be released during relaxation
    - `--relax-no-variable-cell=<false>`: whether the cell should be relaxed (`false`) or just the atomic positions (`true`)
    - `--hp-base=<false>`: whether to calculate U parameters self-consistently before starting the search
    - `--hp-unique=<false>`: whether to recalculate U self-consistently for each unique state
    - `--hp-nq1=<2>`: `nq1` in hp.x input
    - `--hp-nq2=<2>`: `nq2` in hp.x input
    - `--hp-nq3=<2>`: `nq3` in hp.x input
    - `--hp-conv-thr-chi=<1e-4>`: `conv_thr_chi` in hp.x input
    - `--hp-find-atpert=<2>`: `find_atpert` in hp.x input
    - `--hp-U-conv-thr=<0.1>`: threshold for self-consistency of U
    """
    @cast function create(name::String, scf_file::Comonicon.Arg.Path;
                          structure_file::Comonicon.Arg.Path = scf_file,
                          primitive::Bool = false,
                          supercell_a::Int = 1,
                          supercell_b::Int = 1,
                          supercell_c::Int = 1,
                          verbosity::Int = 0,
                          sleep_time::Float64 = 30.0,
                          nrand::Int = 10,
                          unique_thr::Float64 = 0.1,
                          mindist_ratio::Float64 = 0.25,
                          stopping_unique_ratio::Float64 = 0.2,
                          stopping_n_generations::Int =3,
                          Hubbard_maxstep::Int = 100,
                          Hubbard_mixing_beta::Float64 = 0.4,
                          Hubbard_strength::Float64 = 1.0,
                          Hubbard_conv_thr::Float64 = 0.1,
                          electron_maxstep::Int = 500,
                          relax_unique::Bool = false,
                          relax_base::Bool   = false,
                          relax_force_convergence_threshold::Float64  = 1e-3,
                          relax_energy_convergence_threshold::Float64 = 1e-4,
                          relax_ion_dynamics::String  = "bfgs",
                          relax_cell_dynamics::String = "bfgs",
                          relax_no_symmetry::Bool      = false,
                          relax_no_variable_cell::Bool = false,
                          hp_base::Bool= false,
                          hp_unique::Bool=false,
                          hp_nq1::Int = 2,
                          hp_nq2::Int = 2,
                          hp_nq3::Int = 2,
                          hp_conv_thr_chi::Float64 = 1e-6,
                          hp_find_atpert::Int = 1,
                          hp_U_conv_thr::Float64 = 0.1)
                          
        l = setup_search(name, abspath(scf_file.content), abspath(structure_file.content);
                      primitive              = primitive,
                      supercell              = [supercell_a, supercell_b, supercell_c],
                      sleep_time             = sleep_time,
                      nrand                  = nrand,
                      unique_thr             = unique_thr,
                      mindist_ratio          = mindist_ratio,
                      stopping_unique_ratio  = stopping_unique_ratio,
                      stopping_n_generations = stopping_n_generations,
                      Hubbard_maxstep        = Hubbard_maxstep,
                      Hubbard_mixing_beta    = Hubbard_mixing_beta,
                      Hubbard_strength       = Hubbard_strength,
                      Hubbard_conv_thr       = Hubbard_conv_thr,
                      electron_maxstep       = electron_maxstep,
                      relax_unique                 = relax_unique,
                      relax_base                   = relax_base,
                      force_convergence_threshold  = relax_force_convergence_threshold,
                      energy_convergence_threshold = relax_energy_convergence_threshold,
                      ion_dynamics                 = relax_ion_dynamics,
                      cell_dynamics            = relax_cell_dynamics,
                      symmetry                     = !relax_no_symmetry,
                      variable_cell                = !relax_no_variable_cell,
                      hp_base         = hp_base,
                      hp_unique       = hp_unique,
                      hp_nq           = (hp_nq1, hp_nq2, hp_nq3),
                      hp_conv_thr_chi = hp_conv_thr_chi,
                      hp_find_atpert  = hp_find_atpert,
                      hp_U_conv_thr   = hp_U_conv_thr)
        return orchestrator_submit(l; verbosity=verbosity)
    end

    """
    show status of a searcher.

    # Args

    - `name`: name of the searcher
    """
    @cast function status(name::String)
        print(orchestrator_eval("status(\"$name\")")) 
    end

    """
    plot states and save to a file.

    # Args

    - `name`: name of the searcher
    - `outfile`: path where to save the image

    # Options

    - `--exclude-hubbard-energy=<false>`: if false Hubbard energy is subtracted from the total energy
    - `--markersize=<3>`: size of markers
    - `--alpha=<1.0>`: alpha of markers
    """
    @cast function plot(name::String, outfile::Comonicon.Arg.Path; exclude_hubbard_energy::Bool = false, markersize::Int = 3, alpha::Float64=1.0)
        print(orchestrator_eval("plot_states(\"$name\", \"$(abspath(outfile.content))\"; include_hub_energy = $(!exclude_hubbard_energy), markersize=$markersize, alpha=$alpha)"))
    end

    """
    switch searcher mode.

    # Args

    - `name`: name of searcher
    - `mode`: new mode of searcher. Choose from `postprocess`, `search` or `cleanup`
    """
    @cast function mode(name::String, mode::String)
        print(orchestrator_eval("set_mode!(\"$name\", \"$mode\")"))
    end
end
@cast searcher

module server
    using Comonicon
    using ..RomeoDFT: RemoteHPC, orchestrator_eval, pretty_table, Server

    """
    list Servers.
    """
    @cast function list()
        servers = RemoteHPC.load(Server(""))
        
        data = Matrix{String}(undef, length(servers), 3)
        for (i, n) in enumerate(servers)
            data[i, 1] = n
            s = Server(n)
            data[i, 2] = string(RemoteHPC.isalive(s))
            try
                data[i, 3] = string(RemoteHPC.version(s))
            catch
                data[i, 3] = "UNKNOWN/UNREACHABLE"
            end
        end
        pretty_table(data, header = ["NAME", "ALIVE", "VERSION"])
    end
    
    """
    stop Server.

    # Args

    - `name`: name of server
    """
    @cast function stop(name::String)
        if RemoteHPC.exists(RemoteHPC.Server(name=name))
            RemoteHPC.kill(RemoteHPC.Server(name))
            print("Killed Server(\"$name\")")
        else
            @error "No such server"
        end
    end
    
    """
    start Server.

    # Args

    - `name`: name of server
    """
    @cast function start(name::String)
        if RemoteHPC.exists(RemoteHPC.Server(name=name))
            RemoteHPC.start(RemoteHPC.Server(name))
            print("Started Server(\"$name\")")
        else
            @error "No such server"
        end
    end
    
    """
    update Server to latest version of RemoteHPC.

    # Args

    - `name`: name of server
    """
    @cast function update(name::String)
        if RemoteHPC.exists(RemoteHPC.Server(name=name))
            #TODO: Print more stuff
            RemoteHPC.update_RemoteHPC(RemoteHPC.Server(name))
            print("Updated Server(\"$name\")")
        else
            @error "No such server"
        end
    end
end
@cast server


"""
configure interactively Servers, Environments, Executables and PseudoSets.
"""
@cast function configure()
    
    function ask_input(msg)
        t = ""
        return t
    end

    choices = ["searchers directory", "Storables (Servers, Execs, PseudoSets, etc...)"]
    
    id = request("What would you like to configure?", RadioMenu(choices))
    if id == -1
        return
    end

    if id == 1
        println("Current searchers directory: ", searchers_dir())
        id = request("Change?", RadioMenu(["yes", "no"]))
        if id == 2
            return
        end
        
        d = ""
        while isempty(d)
            println("Provide a directory where Searchers and their temporary files will be stored:")
            t = ""
            while isempty(t)
                t = readline()
            end
            t = abspath(t)
            if isdir(t)
                d = t
            else
                id = request("$t no such directory. Create?", RadioMenu(["yes", "no"]))
                if id == 1
                    mkdir(t)
                    d = t
                else
                    @error "$t no such directory. Try again"
                end
            end
        end
        open(config_path("config.toml"), "w") do f
            TOML.print(f, Dict("searchers_directory" => d))
        end
    elseif id == 2
        RemoteHPC.configure()
    end
end

@main
