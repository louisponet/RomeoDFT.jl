var documenterSearchIndex = {"docs":
[{"location":"components/#Components","page":"Components","title":"Components","text":"","category":"section"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Here we go through the steps required to install all the parts and software that are needed to perform global searches with Quantum Espresso and GOMEO. Afterwards continue with the Tutorial.","category":"page"},{"location":"installation/#Quantum-Espresso","page":"Installation","title":"Quantum Espresso","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Install the version of Quantum Espresso with constrained occupations on each cluster or computer you want to perform the calculations on. This can be done through:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"git clone https://gitlab.com/louisponet/q-e/ qe-occupations\ncd qe-occupations\ngit checkout broyden_constraints","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Followed by the standard commands that are required to install Quantum Espresso.","category":"page"},{"location":"installation/#Julia","page":"Installation","title":"Julia","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Download and install the latest version of Julia on the computer you are going to use to orchestrate the search, e.g. your personal workstation.","category":"page"},{"location":"installation/#RomeoDFT.jl","page":"Installation","title":"RomeoDFT.jl","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Next we install the main driver software which for now is not publicly available, so some steps are required.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Begin by starting a julia REPL session.","category":"page"},{"location":"installation/#First-Time","page":"Installation","title":"First Time","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Execute the following line:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"mkpath(joinpath(homedir(), \".julia/dev\"))","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"then press ; which activates the shell> REPL mode, this allows you to execute standard shell commands like in bash. Execute:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"cd ~/.julia/dev\ngit clone https://github.com/louisponet/RomeoDFT.jl RomeoDFT","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Press backspace and ] which activates the pkg> REPL mode, which allows you to install packages. Execute","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"dev RomeoDFT\nbuild RomeoDFT","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"The package is now being intialized, sit back grab a coffee, this takes a long time.","category":"page"},{"location":"installation/#Update","page":"Installation","title":"Update","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Press ; and execute","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"cd ~/.julia/dev/RomeoDFT\ngit pull","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Then enter the pkg> mode and execute","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"up RomeoDFT\nbuild RomeoDFT","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"again, this will take a while...","category":"page"},{"location":"installation/#GOMEO-executable","page":"Installation","title":"GOMEO executable","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"After having performed the above steps, there should be the gomeo executable in ~/.julia/bin. You can use it as, e.g. ~/.julia/bin/gomeo -h, or add ~/.julia/bin to your PATH, after which you can use it from anywhere simply as gomeo -h.","category":"page"},{"location":"installation/#RemoteHPC-Configuration","page":"Installation","title":"RemoteHPC Configuration","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Next we configure the Servers, Execs and Environments that represent where, what and how things can be ran, in our case the qe-occupations/bin/pw.x executable.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"From a normal shell, execute","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo configure","category":"page"},{"location":"installation/#Server","page":"Installation","title":"Server","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Use the arrow keys and enter to select Server. A Server represents a computer or cluster where executables can be ran. Your local computer is configured automatically so if you want to run everything locally, you can skip this step. Fill in the name, username and domain that correspond to the ssh commands you use to connect to the cluster where you installed qe-occupations (e.g. ssh lponet@cscs.ch would mean username = lponet, domain = cscs.ch). If you installed it on multiple clusters, you can add them one after the other.","category":"page"},{"location":"installation/#Exec","page":"Installation","title":"Exec","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Next select Exec from the menu. Read carefully the documentation that is printed, enter a name (for example qe-occupations), select the Server that corresponds to where you installed Quantum Espresso. After pressing enter, an editor will open with a representation of the Exec with some fields already filled in. Complete the remaining fields, especially path which should be the directory of the qe-occupations/bin/pw.x executable on the remote Server. Again, read the Documentation that is printed below the fields to help filling out the remaining information.","category":"page"},{"location":"installation/#Environment","page":"Installation","title":"Environment","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"An environment essentially represents the #SBATCH commands of a slurm job script, or similar. Carefully read the documentation and fill out the fields similar to above.","category":"page"},{"location":"installation/#Pseudo-Potentials","page":"Installation","title":"Pseudo Potentials","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Finally we set up a set of pseudopotentials. For our intention, we install them locally. First download and unpack your favorite set of pseudopotentials, then configure a PseudoSet. Behind the scenes it will read through the directory you specified and resolve element names for the files in that directory that have the .UPF extension.","category":"page"},{"location":"installation/#Tutorial","page":"Installation","title":"Tutorial","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"In order to begin a ground state search, we need a template scf file that will be used to fill out most of the input parameters. An example for MnO can be found in the test/assets/ subdirectory.","category":"page"},{"location":"installation/#Creation","page":"Installation","title":"Creation","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Once you have such a file for your case (or if you just want to try with MnO), execute:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"cd ~/.julia/dev/RomeoDFT\ngomeo searcher create MnO_tutorial test/assets/MnO_scf.in --sleep-time 5","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"the sleep-time flag is used to specify how often the orchestrator polls for new results etc.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"In this case, the structure definition will also be taken from the template file. If you want to use a different pw.x input or cif file from which to extract the structure, you can use","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher create MnO_tutorial test/assets/MnO_scf.in -s <structure_file> --sleep-time 5","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"There will be an interactive menu to further set up the searcher, after which it will be started by the Orchestrator.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"For further information on the searcher create command you can run","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher create -h","category":"page"},{"location":"installation/#Monitoring","page":"Installation","title":"Monitoring","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"To get an overview of all the loaded Searchers you can execute","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo orchestrator status","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"To view information of a specific Searcher, execute:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher status <name>","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"where <name> in our case would be MnO_tutorial.","category":"page"},{"location":"installation/#Stopping/Starting","page":"Installation","title":"Stopping/Starting","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"To stop a running Searcher you can execute:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher stop <name>","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"and to restart it","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher start <name>","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"If no Searcher with name was loaded yet, it will be first.","category":"page"},{"location":"installation/#Changing-mode","page":"Installation","title":"Changing mode","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"By default a Searcher will start in the search mode. This means that it is actively exploring the energy landscape in order to find unique metastable states. It will continue to do so until the amount of unique states versus trials goes below the value set with the unique-ratio flag during gomeo searcher create. Afterwards, it will switch to the postprocess mode, meaning that it will not generate new trials, but simply process all the remaining trials and finish all the pending calculations. If the ratio of unique/trial happens to exceed the unique-ratio value after finishing the postprocessing, the mode will switch back to search.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Sometimes, however, it is desirable to manually switch the mode in order to limit the search or simply stop a searcher from exploring further but cleaning up everything. This can be done using:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher mode <name> <mode>","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"where <mode> is one of the following:","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"search: Active searching through generating new trials.\npostprocess: No new trials are generated, but all remaining ones are submitted and cleaned up.\ncleanup: No new trials are generated, no new jobs are submitted (i.e. previously generated but not yet submitted trials will not be submitted), and searcher finishes after all currently running calculations finish.","category":"page"},{"location":"installation/#Plot-States","page":"Installation","title":"Plot States","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"An overview of the different states that are found so far can be plotted using","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher plot <name> <outfile>","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"which will produce an image and save it to outfile displaying the states.","category":"page"},{"location":"installation/#Other-commands","page":"Installation","title":"Other commands","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"gomeo searcher unload <name>: Stops, saves and removes a Searcher from the Orchestrator\ngomeo searcher load <name>: Loads a Searcher.","category":"page"},{"location":"#RomeoDFT.jl","page":"Home","title":"RomeoDFT.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation on how to run global searchers with RomeoDFT.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Start with Installation and then follow the Tutorial.","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In order to begin a ground state search, we need a template scf file that will be used to fill out most of the input parameters. An example for MnO can be found in the test/assets/ subdirectory.","category":"page"},{"location":"tutorial/#Creation","page":"Tutorial","title":"Creation","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Once you have such a file for your case (or if you just want to try with MnO), execute:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"cd ~/.julia/dev/RomeoDFT\ngomeo searcher create MnO_tutorial test/assets/MnO_scf.in --sleep-time 5","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"the sleep-time flag is used to specify how often the orchestrator polls for new results etc.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"In this case, the structure definition will also be taken from the template file. If you want to use a different pw.x input or cif file from which to extract the structure, you can use","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher create MnO_tutorial test/assets/MnO_scf.in -s <structure_file> --sleep-time 5","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"There will be an interactive menu to further set up the searcher, after which it will be started by the Orchestrator.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For further information on the searcher create command you can run","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher create -h","category":"page"},{"location":"tutorial/#Monitoring","page":"Tutorial","title":"Monitoring","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To get an overview of all the loaded Searchers you can execute","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo orchestrator status","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To view information of a specific Searcher, execute:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher status <name>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where <name> in our case would be MnO_tutorial.","category":"page"},{"location":"tutorial/#Stopping/Starting","page":"Tutorial","title":"Stopping/Starting","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To stop a running Searcher you can execute:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher stop <name>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"and to restart it","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher start <name>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"If no Searcher with name was loaded yet, it will be first.","category":"page"},{"location":"tutorial/#Changing-mode","page":"Tutorial","title":"Changing mode","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"By default a Searcher will start in the search mode. This means that it is actively exploring the energy landscape in order to find unique metastable states. It will continue to do so until the amount of unique states versus trials goes below the value set with the unique-ratio flag during gomeo searcher create. Afterwards, it will switch to the postprocess mode, meaning that it will not generate new trials, but simply process all the remaining trials and finish all the pending calculations. If the ratio of unique/trial happens to exceed the unique-ratio value after finishing the postprocessing, the mode will switch back to search.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Sometimes, however, it is desirable to manually switch the mode in order to limit the search or simply stop a searcher from exploring further but cleaning up everything. This can be done using:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher mode <name> <mode>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where <mode> is one of the following:","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"search: Active searching through generating new trials.\npostprocess: No new trials are generated, but all remaining ones are submitted and cleaned up.\ncleanup: No new trials are generated, no new jobs are submitted (i.e. previously generated but not yet submitted trials will not be submitted), and searcher finishes after all currently running calculations finish.","category":"page"},{"location":"tutorial/#Plot-States","page":"Tutorial","title":"Plot States","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"An overview of the different states that are found so far can be plotted using","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher plot <name> <outfile>","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"which will produce an image and save it to outfile displaying the states.","category":"page"},{"location":"tutorial/#Other-commands","page":"Tutorial","title":"Other commands","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gomeo searcher unload <name>: Stops, saves and removes a Searcher from the Orchestrator\ngomeo searcher load <name>: Loads a Searcher.","category":"page"}]
}
