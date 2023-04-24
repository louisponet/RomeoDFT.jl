var documenterSearchIndex = {"docs":
[{"location":"components/#Components","page":"Components","title":"Components","text":"","category":"section"},{"location":"orchestrator/#Orchestrator","page":"Orchestrator","title":"Orchestrator","text":"","category":"section"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"CurrentModule = RomeoDFT","category":"page"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"The Orchestrator can be thought of as a locally running daemon that handles running the various Searchers. There are three main commands that are useful to interact with it (from a standard shell):","category":"page"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"romeo orchestrator start: starts the Orchestrator and loads and starts the searchers that were previously running\nromeo orchestrator stop: stops and saves the currently loaded Searchers and gracefully stops the Orchestrator\nromeo orchestrator status: shows an overview of the currently loaded Searchers","category":"page"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"For more information use -h after these commands in the shell.","category":"page"},{"location":"orchestrator/#Advanced","page":"Orchestrator","title":"Advanced","text":"","category":"section"},{"location":"orchestrator/#RemoteREPL","page":"Orchestrator","title":"RemoteREPL","text":"","category":"section"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"Behind the scenes the Orchestrator runs in a RemoteREPL, which can be accessed to interact directly with the loaded Searchers. This can be achieved by starting a julia REPL and executing the following commands:","category":"page"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"using RomeoDFT\nconnect_orchestrator()","category":"page"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"You will see that now by pressing > you get a pink julia prompt, which means you're executing commands inside the Orchestrator REPL. The useful commands here are:","category":"page"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"start_searcher(\"<searcher_name>\"; verbosity = <number>): starts a Searcher with the provided name with a given logging verbosity (higher = more logged).\nstop_searcher(\"<searcher_name>\"): stops a Searcher with the name\nload_searcher: simply loads a Searcher from disk if no function is provided as the first argument. See load_searcher documentation for more info.","category":"page"},{"location":"orchestrator/#Submitting-Searcher","page":"Orchestrator","title":"Submitting Searcher","text":"","category":"section"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"In order to create and submit a Searcher manually you can use setup_search and orchestrator_submit in sequence.","category":"page"},{"location":"orchestrator/#Reference","page":"Orchestrator","title":"Reference","text":"","category":"section"},{"location":"orchestrator/","page":"Orchestrator","title":"Orchestrator","text":"start_searcher\nstop_searcher\nload_searcher\nconnect_orchestrator\norchestrator_submit","category":"page"},{"location":"orchestrator/#RomeoDFT.start_searcher","page":"Orchestrator","title":"RomeoDFT.start_searcher","text":"start_searcher(name; verbosity=0)\n\nStarts a Searcher and sets the logging verbosity. Higher verbosity means more will be printed to the log file.\n\n\n\n\n\n","category":"function"},{"location":"orchestrator/#RomeoDFT.stop_searcher","page":"Orchestrator","title":"RomeoDFT.stop_searcher","text":"stop_searcher(name)\n\nStops a Searcher.\n\n\n\n\n\n","category":"function"},{"location":"orchestrator/#RomeoDFT.load_searcher","page":"Orchestrator","title":"RomeoDFT.load_searcher","text":"load_searcher([f::Function,] name)\n\nLoads a Searcher to the Orchestrator, and executes f on it if it is provided. This can be useful to perform manual tasks or retrieve more granular information about the Searcher.\n\nExample\n\nload_searcher(\"MnO\") do s\n    # execute code on Searcher s which represents the MnO searcher\nend\n\n\n\n\n\n","category":"function"},{"location":"orchestrator/#RomeoDFT.orchestrator_submit","page":"Orchestrator","title":"RomeoDFT.orchestrator_submit","text":"orchestrator_submit(searcher)\n\nSave and submit searcher to the running Orchestrator.\n\n\n\n\n\n","category":"function"},{"location":"getting_started/#Installation","page":"Getting Started","title":"Installation","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Here we go through the steps required to install all the parts and software that are needed to perform global searches with Quantum Espresso and ROMEO. Afterwards continue with the Tutorial.","category":"page"},{"location":"getting_started/#Quantum-Espresso","page":"Getting Started","title":"Quantum Espresso","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Install the version of Quantum Espresso with constrained occupations on each cluster or computer you want to perform the calculations on. This can be done through:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"git clone https://gitlab.com/louisponet/q-e/ qe-occupations\ncd qe-occupations\ngit checkout broyden_constraints","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Followed by the standard commands that are required to install Quantum Espresso.","category":"page"},{"location":"getting_started/#Julia","page":"Getting Started","title":"Julia","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Download and install the latest version of Julia on the computer you are going to use to orchestrate the search, e.g. your personal workstation.","category":"page"},{"location":"getting_started/#RomeoDFT.jl","page":"Getting Started","title":"RomeoDFT.jl","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Next we install the main driver software which for now is not publicly available, so some steps are required.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Begin by starting a julia REPL session, and execute:","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"using Pkg\nPkg.add(\"RomeoDFT\")","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"or press the ] key to enter the julia pkg REPL mode and type add RomeoDFT. The package is now being intialized, sit back grab a coffee, this takes a long time (~5-10 min).","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"To update to a newer version you can do","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"using Pkg\nPkg.update(\"RomeoDFT\")","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"or up RomeoDFT in the pkg REPL mode.","category":"page"},{"location":"getting_started/#ROMEO-executable","page":"Getting Started","title":"ROMEO executable","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"After having performed the above steps, there should be the romeo executable in ~/.julia/bin. You can use it as, e.g. ~/.julia/bin/romeo -h.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"note: Note\nIt is strongly recommended to add ~/.julia/bin to your PATH environment variable by changing ~/.bashrc or similar.","category":"page"},{"location":"getting_started/#RemoteHPC-Configuration","page":"Getting Started","title":"RemoteHPC Configuration","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Next we configure the Servers, Execs and Environments that represent where, what and how executables can be ran on the remote clusters. In our case the qe-occupations/bin/pw.x executable.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"note: Note\nWe will be using the default editor to write the configuration files in the following steps. If the EDITOR environment variable is not set, this will be nano. Again it is advised to set it with export EDITOR=<vim vscode or any other editor> or in your ~/.bashrc.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"From inside a shell, execute","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"romeo configure","category":"page"},{"location":"getting_started/#Server","page":"Getting Started","title":"Server","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Use the arrow keys and enter to select Server. A Server represents a computer or cluster where executables can be ran. Your local computer is configured automatically so if you want to run everything locally, you can skip this step. Give the server a name that is used to load it later, and insert your username and domain that correspond to the ssh commands you use to connect to the cluster where you installed qe-occupations (e.g. ssh lponet@cscs.ch would mean username = lponet, domain = cscs.ch).","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"If you plan to run on multiple clusters, you can add them one after the other. For now, however, each Searcher runs on a single Server.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"If you configured a new Server, go out of the configuration at this point and run romeo server start <name of new server> then run romeo configure again.","category":"page"},{"location":"getting_started/#Exec","page":"Getting Started","title":"Exec","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Next select Exec from the menu. Read carefully the documentation that is printed, enter a name (for example qe-occupations), select the Server that corresponds to where you installed the patched QuantumESPRESSO. After pressing enter, an editor will open with a representation of the Exec with some fields already filled in. Complete the remaining fields, especially path which should be the absolute directory of the qe-occupations/bin/pw.x executable on the remote Server. Again, read the Documentation that is printed below the fields to help filling out the remaining information.","category":"page"},{"location":"getting_started/#Environment","page":"Getting Started","title":"Environment","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"An environment essentially represents the #SBATCH commands of a slurm job script, or similar. Carefully read the documentation and fill out the fields similar to above.","category":"page"},{"location":"getting_started/#Pseudo-Potentials","page":"Getting Started","title":"Pseudo Potentials","text":"","category":"section"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"Finally we set up a set of pseudopotentials. For our intention, we install them locally. First download and unpack your favorite set of pseudopotentials into a local directoyr, then supply the directory during the configuration of a PseudoSet. Behind the scenes it will read through the directory you specified and resolve element names for the files in that directory that have the .UPF extension.","category":"page"},{"location":"getting_started/","page":"Getting Started","title":"Getting Started","text":"You're now ready to start the Orchestrator and run your first Searcher.","category":"page"},{"location":"searcher/#Searchers","page":"Searchers","title":"Searchers","text":"","category":"section"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"CurrentModule = RomeoDFT","category":"page"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"A Searcher is the main object representing all the data related with a global search. Most interactions can be performed from the command line using the romeo searcher subset of commands. Each Searcher is assigned a label or name with which it can be identified.","category":"page"},{"location":"searcher/#Data","page":"Searchers","title":"Data","text":"","category":"section"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"The name also references the sub directory where all the data will be stored, which can be retrieved using dirname. This data consists of:","category":"page"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"$(DATABASE_VERSION)/ledger.jld2: stores all the search data that is kept in memory during the search. See load and save\njob_backups: a directory storing all the input and output data of the various calculations ran during the search\nlog.log: file that keeps a log of everything that happens during the search. The amount of information logged can be tuned using different levels passed to --verbosity <n> flag of romeo searcher create or romeo searcher start\nHTTP.log: file holding all HTTP related logs\nreport.out: a report generated when a Searcher is finished","category":"page"},{"location":"searcher/#CommandLine","page":"Searchers","title":"CommandLine","text":"","category":"section"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"Most interaction with a Searcher can be done through the command line using the romeo searcher subcommands (use -h to show more info on commands and flags). The most important ones are:","category":"page"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"romeo searcher create: interactive wizard to initialize and start a new Searcher, see setup_search\nromeo searcher start: starts a Searcher\nromeo searcher stop: stops a Searcher\nromeo searcher status: shows a status overview\nromeo searcher unload: unloads a Searcher from the orchestrator\nromeo searcher log: opens the log.log file in the editor specified by the \\$EDITOR environment variable\nromeo searcher plot: creates a standard State plot\nromeo searcher mode: sets the mode, see Modes for more information","category":"page"},{"location":"searcher/#Life-Cycle","page":"Searchers","title":"Life Cycle","text":"","category":"section"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"The lifecycle of a Searcher can be summarized by the following steps:","category":"page"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"Potential geometry and U parameter optimization when --hp-base or --relax-base are passed to romeo searcher create\nInitial unconstrained scf calculation to determine the number of electrons in the target manifold.\nRandom search with --nrand trials\nBased on the unique States found in 3. (determined by --unique-thr) perform midpoint intersections (n = (n1 + n2)/2) between each pair of unique states\nrepeat 4. for each newly found unique state\nif no new intersections can be created and the stopcondition is not yet reached go to step 3.\nif the stop condition determined by --stopping-unique-ratio and --stopping-n-generations is reached switch to postprocess mode\nif all jobs are finished and the stop condition is still satisfied, stop the Searcher and print results to report.out, otherwise switch back to 3.","category":"page"},{"location":"searcher/#Modes","page":"Searchers","title":"Modes","text":"","category":"section"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"There several execution modes that a Searcher can be in:","category":"page"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"search: search mode characterized by performing intersection and random search, and postprocessing on unique States.\npostprocess: when the stopcondition is satisfied, no new intersections or random trials will be generated, and only postprocessing actions will be taken\ncleanup: no new trials or postprocessing jobs will be created and submitted, only cleanup of the currently running jobs will be performed, after which the Searcher will finish\nmanual: only used when performing manual tasks and you don't want any searching or postprocessing to happen automatically","category":"page"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"These can be set using romeo searcher mode. ","category":"page"},{"location":"searcher/#Reference","page":"Searchers","title":"Reference","text":"","category":"section"},{"location":"searcher/","page":"Searchers","title":"Searchers","text":"Searcher\nsetup_search\nload(::Searcher)\nsave(::Searcher)\nsssp_distance","category":"page"},{"location":"searcher/#RomeoDFT.Searcher","page":"Searchers","title":"RomeoDFT.Searcher","text":"Searcher(name)\n\nRepresents a global search with label name. \n\n\n\n\n\n","category":"type"},{"location":"searcher/#RomeoDFT.setup_search","page":"Searchers","title":"RomeoDFT.setup_search","text":"setup_search(name, scf_file, structure_file=scf_file; kwargs...)\n\nCreates and saves a new Searcher with name taking the pw input parameters from the template scf_file and the structure from structure_file, which can be either a pw input or a .cif file. This is the backend method used for the romeo searcher create from the command line.\n\nStructure Kwargs\n\nprimitive=false: set this to true to first try to find the primitive unit cell\nsupercell=[1,1,1]: number specifying the number of unit cells along a, b and c direction\npseudoset=nothing: label (string) of pseudoset to use, if this is nothing it will trigger interactive setup\nU_values=nothing: Dict{Symbol, Float64} of U values to use for the constrained manifolds if it is nothing it will trigger interactive setup\n\nControl Kwargs\n\nverbosity=0: the logging verbosity, higher = more. See Data for more info\nsleep_time=30: time in seconds between update polls of the Searcher\nserver=nothing: label of Server on which to run everything\nexec=nothing: label of the pw.x executable on server to use for the search\nenvironment=nothing: label of the Environment to use for running all calculations\npriority=nothing: number signifying the priority of the Searcher\n\nSearch Kwargs\n\nnrand=10: how many random trials should be performed in a random search generation\nunique_thr=0.1: threshold that determines the uniqueness of electronic states (uses sssp_distance)\nmindist_ratio=0.25: minimum distance a trial should have to previous trials and unique states,                       defined as the ratio of the mean distance between current unique states\nstopping_unique_ratio=0.2: the ratio of new found unique states to trials in a generation below which the searching will stop\nstopping_n_generations=3: the search will stop if for this amount of successive generations the unique to trial ratio is lower than the stopping_unique_ratio\nHubbard_maxstep=100: maximum constraining steps\nHubbard_mixing_beta=0.4: mixing used to update the constraints during scf iterations\nHubbard_strength=1.0: strength of the constraining potential\nHubbard_conv_thr=0.1: threshold euclidean distance per atom between trial and current occupation matrices after which the constraints are released\nelectron_maxstep=500: see QE documentation\n\nPre/Post processing Kwargs:\n\nrelax_unique=false: whether a structural relaxation should be ran on each unique state\nrelax_base=false: whether a relaxation should be ran so the search happens for the relaxed structure\nrelax_force_convergence_threshold=1e-3: forc_conv_thr in QE\nrelax_energy_convergence_threshold=1e-4: energy_conv_thr in QE\nrelax_ion_dynamics=\"bfgs\": ion_dynamics in QE\nrelax_cell_dynamics=\"bfgs\": cell_dynamics in QE\nrelax_no_symmetry=false: whether symmetry should be released during relaxation\nrelax_no_variable_cell=false: whether the cell should be relaxed (false) or just the atomic positions (true)\nhp_base=false: whether to calculate U parameters self-consistently before starting the search\nhp_unique=false: whether to recalculate U self-consistently for each unique state\nhp_nq=(2,2,2): nq1, nq2, nq3 in hp.x input\nhp_conv_thr_chi=1e-4: conv_thr_chi in hp.x input\nhp_find_atpert=2: find_atpert in hp.x input\nhp_U_conv_thr=0.1: threshold for self-consistency of U\n\n\n\n\n\n","category":"function"},{"location":"searcher/#RemoteHPC.load-Tuple{Searcher}","page":"Searchers","title":"RemoteHPC.load","text":"load(s::Searcher; version)\n\nLoads the Searcher data from disk, from the ledger.jld2 that is associated with the latest Database version if version is unspecified. See Data for more info.\n\nExample\n\ns = load(Searcher(\"my_searcher_name\"))\n\n\n\n\n\n","category":"method"},{"location":"searcher/#RemoteHPC.save-Tuple{Searcher}","page":"Searchers","title":"RemoteHPC.save","text":"save(s::Searcher)\n\nSaves s to the standard location as specified in Data.\n\n\n\n\n\n","category":"method"},{"location":"searcher/#RomeoDFT.sssp_distance","page":"Searchers","title":"RomeoDFT.sssp_distance","text":"sssp_distance(bands1, bands2, fermi)\n\nCalculates\n\nsqrtfracsum_mk f_mk (varepsilon^1_mk - varepsilon^2_mk + Delta)^2sum_mkf_mk\n\nfor bands1 and bands2 flattened bandstructures. fermi is used to crudely determine f_mk and Delta is optimized to minimize the distance. \n\n\n\n\n\n","category":"function"},{"location":"states/#States","page":"States","title":"States","text":"","category":"section"},{"location":"states/","page":"States","title":"States","text":"CurrentModule = RomeoDFT","category":"page"},{"location":"states/","page":"States","title":"States","text":"A State represents a set of occupied local manifolds, one for each ion to which occupation constraints are applied, or that have a Hubbard U parameter. Usually no direct interaction with these is needed as the search is happening, but they are used for some postprocessing analysis.","category":"page"},{"location":"states/","page":"States","title":"States","text":"They hold:","category":"page"},{"location":"states/","page":"States","title":"States","text":"occupations: the occupation matrices for each ion\ntotoccs: total electron occupation of each ion\nmagmoms: magnetic moments defined as tr(up) - tr(down)\neigvals: the eigen values for each occupation matrix\neigvecs: the eigen vectors for each occupation matrix","category":"page"},{"location":"states/#Euclidean-distance","page":"States","title":"Euclidean distance","text":"","category":"section"},{"location":"states/","page":"States","title":"States","text":"The Euclidean distance metric is defined between states, which essentially does","category":"page"},{"location":"states/","page":"States","title":"States","text":"sqrtsum_I sum_alphabeta(n^I_1 alphabeta - n^I_2 alphabeta)^2","category":"page"},{"location":"states/","page":"States","title":"States","text":"and can be called like","category":"page"},{"location":"states/","page":"States","title":"States","text":"RomeoDFT.Euclidean()(state1, state2)","category":"page"},{"location":"states/#Refeference","page":"States","title":"Refeference","text":"","category":"section"},{"location":"states/","page":"States","title":"States","text":"State\ngenerate_Hubbard_occupations\ngenerate_starting_ns_eigenvalue","category":"page"},{"location":"states/#RomeoDFT.State","page":"States","title":"RomeoDFT.State","text":"State\n\nRepresents the local state of a system. This is given by the occupation matrices of the local orbitals. These are usually the valence shells for example to which the +U correction is applied in DFT + U calculations.\n\nExamples\n\nA State can be constructed from the :Hubbard entry from a QE pw output like:\n\nusing DFControl\nusing RomeoDFT\n\no = DFC.FileIO.qe_parse_pw_output(\"<pw_lda_U_output_file_with_Hubbard_blocks>\")\ns = RomeoDFT.State(o[:Hubbard][end])\n\nor using a Vector with MagneticArrays from DFWannier.jl.\n\nusing DFWannier\nusing RomeoDFT\nusing LinearAlgebra\n\nup   = diagm(0 => ones(5))\ndown = diagm(0 => ones(5))\noccs = [DFWannier.ColinMatrix(up, down)]\n\ns = RomeoDFT.State(occs)\n\nUsing RomeoDFT.generate_Hubbard_occupations this can generate the :Hubbard_occupations scf input parameter which will be used as the target during a constrained scf calculation. RomeoDFT.generate_starting_ns_eigenvalue generates the :starting_ns_eigenvalue parameter instead.\n\n\n\n\n\n","category":"type"},{"location":"states/#RomeoDFT.generate_Hubbard_occupations","page":"States","title":"RomeoDFT.generate_Hubbard_occupations","text":"generate_Hubbard_occupations(state)\n\nGenerates the :Hubbard_occupations entry for a constrained scf calculation.\n\n\n\n\n\n","category":"function"},{"location":"states/#RomeoDFT.generate_starting_ns_eigenvalue","page":"States","title":"RomeoDFT.generate_starting_ns_eigenvalue","text":"starting_ns_eigenvalue(state)\n\nGenerates the :starting_ns_eigenvalue entry for a constrained scf calculation.\n\n\n\n\n\n","category":"function"},{"location":"#RomeoDFT.jl","page":"Home","title":"RomeoDFT.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This software package allows one to perform a Robust Occupation Matrix Energy Optimization of a crystalline structure using Density Functional Theory.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It is most useful in materials with complex magnetic configurations or where the local orbital physics of one or more atomic species plays a crucial role in describing the ground state. In many of these cases the energy surface of the DFT functional sports many local minima, which in theory should be extensively explored in order to find the global minimum associated with the ground state.","category":"page"},{"location":"","page":"Home","title":"Home","text":"RomeoDFT.jl facilitates a fully automated global search, exploring many different occupations of local atomic orbitals (target states) in order to reliably determine the ground state without much human interaction.","category":"page"},{"location":"","page":"Home","title":"Home","text":"It uses a patched version of QuantumESPRESSO to run the DFT calculations for each of the target states.","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"\"Press of the button\" determination of the ground state starting from a single template input file with the structure and input parameters\nSeamless interfacing with HPC clusters using RemoteHPC.jl and DFControl.jl\nRun geometry optimization for the initial structure as well as any found self-consistent state associated with a local minimum\nAutomatically and fully self-consistently determine the Hubbard U parameter for the initial structure as well as any self-consistent state associated with a local minimum\nPostProcessing steps including non-selfconsistent, bandstructure and dos/projected dos calculations\nRobust data storage using JLD2.jl with versioning\nCommand line interface facilitating most of the tasks\nBuilt with ECS (see Overseer.jl making it very extensible","category":"page"},{"location":"#Acknowledgements","page":"Home","title":"Acknowledgements","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This work has been funded by H2020-OpenModel and EPFL.","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In order to begin a ground state search, we need a template scf file that will be used to fill out most of the input parameters. An example for MnO can be found in the test/assets/ subdirectory.","category":"page"},{"location":"tutorial/#Creation","page":"Tutorial","title":"Creation","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Once you have such a file for your case (or if you just want to try with MnO), execute:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"cd ~/.julia/dev/RomeoDFT\nromeo searcher create MnO_tutorial test/assets/MnO_scf.in --sleep-time 5","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"the sleep-time flag is used to specify how often the orchestrator polls for new results etc.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In this case, the structure definition will also be taken from the template file. If you want to use a different pw.x input or cif file from which to extract the structure, you can use","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher create MnO_tutorial test/assets/MnO_scf.in -s <structure_file> --sleep-time 5","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"There will be an interactive menu to further set up the Searcher, after which it will be started by the Orchestrator.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For further information on the searcher create command you can run","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher create -h","category":"page"},{"location":"tutorial/#Monitoring","page":"Tutorial","title":"Monitoring","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To get an overview of all the loaded Searchers you can execute","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo orchestrator status","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To view information of a specific Searcher, execute:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher status <name>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where <name> in our case would be MnO_tutorial.","category":"page"},{"location":"tutorial/#Stopping/Starting","page":"Tutorial","title":"Stopping/Starting","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To stop a running Searcher you can execute:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher stop <name>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"and to restart it","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher start <name>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"If no Searcher with name was loaded yet, it will be first.","category":"page"},{"location":"tutorial/#Changing-mode","page":"Tutorial","title":"Changing mode","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"By default a Searcher will start in the search mode. This means that it is actively exploring the energy landscape in order to find unique metastable states. It will continue to do so until the amount of unique states versus trials goes below the value set with the unique-ratio flag during romeo searcher create. Afterwards, it will switch to the postprocess mode, meaning that it will not generate new trials, but simply process all the remaining trials and finish all the pending calculations. If the ratio of unique/trial happens to exceed the unique-ratio value after finishing the postprocessing, the mode will switch back to search.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Sometimes, however, it is desirable to manually switch the mode in order to limit the search or simply stop a searcher from exploring further but cleaning up everything. This can be done using:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher mode <name> <mode>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"See Modes for more info.","category":"page"},{"location":"tutorial/#Plot-States","page":"Tutorial","title":"Plot States","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"An overview of the different states that are found so far can be plotted using","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"romeo searcher plot <name> <outfile>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"which will produce an image and save it to outfile displaying the states.","category":"page"},{"location":"tutorial/#Other-commands","page":"Tutorial","title":"Other commands","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"See CommandLine for more commands.","category":"page"}]
}