# Installation

Here we go through the steps required to install all the parts and software that are needed to perform global searches with Quantum Espresso and GOMEO.
Afterwards continue with the [`Tutorial`](@ref Tutorial).

## Quantum Espresso

Install the [version of Quantum Espresso with constrained occupations](https://gitlab.com/louisponet/q-e/-/tree/broyden_constraints) on each cluster or computer you want to perform the calculations on.
This can be done through:

```bash
git clone https://gitlab.com/louisponet/q-e/ qe-occupations
cd qe-occupations
git checkout broyden_constraints
```
Followed by the [standard commands](https://www.quantum-espresso.org/Doc/user_guide/node7.html) that are required to install Quantum Espresso.

## Julia
Download and install the latest version of [Julia](https://julialang.org/downloads/) on the computer you are going to use to orchestrate the search, e.g. your personal workstation.

## RomeoDFT.jl

Next we install the main driver software which for now is not publicly available, so some steps are required.

Begin by starting a julia REPL session.

### First Time
Execute the following line:
```julia
mkpath(joinpath(homedir(), ".julia/dev"))
```
then press `;` which activates the `shell>` REPL mode, this allows you to execute standard shell commands like in `bash`.
Execute:
```bash
cd ~/.julia/dev
git clone https://github.com/louisponet/RomeoDFT.jl RomeoDFT
```
Press backspace and `]` which activates the `pkg>` REPL mode, which allows you to install packages.
Execute
```
dev RomeoDFT
build RomeoDFT
```
The package is now being intialized, sit back grab a coffee, this takes a long time.

### Update
Press `;` and execute
```
cd ~/.julia/dev/RomeoDFT
git pull
```
Then enter the `pkg>` mode and execute
```
up RomeoDFT
build RomeoDFT
```
again, this will take a while...

## GOMEO executable
After having performed the above steps, there should be the `gomeo` executable in `~/.julia/bin`.
You can use it as, e.g. `~/.julia/bin/gomeo -h`, or add `~/.julia/bin` to your `PATH`, after which you can use it from anywhere simply as `gomeo -h`.

## RemoteHPC Configuration
Next we configure the `Servers`, `Execs` and `Environments` that represent where, what and how things can be ran, in our case the `qe-occupations/bin/pw.x` executable.

From a normal `shell`, execute
```
gomeo configure
```
### Server
Use the arrow keys and enter to select `Server`. A `Server` represents a computer or cluster where executables can be ran.
Your local computer is configured automatically so if you want to run everything locally, you can skip this step.
Fill in the `name`, `username` and `domain` that correspond to the `ssh` commands you use to connect to the cluster where you installed `qe-occupations` (e.g. `ssh lponet@cscs.ch` would mean `username = lponet`, `domain = cscs.ch`).
If you installed it on multiple clusters, you can add them one after the other.

### Exec
Next select `Exec` from the menu. Read carefully the documentation that is printed, enter a `name` (for example `qe-occupations`), select the `Server` that corresponds to where you installed Quantum Espresso.
After pressing enter, an editor will open with a representation of the `Exec` with some fields already filled in.
Complete the remaining fields, especially `path` which should be the directory of the `qe-occupations/bin/pw.x` executable on the remote `Server`.
Again, read the Documentation that is printed below the fields to help filling out the remaining information.

### Environment
An environment essentially represents the #SBATCH commands of a slurm job script, or similar.
Carefully read the documentation and fill out the fields similar to above.

### Pseudo Potentials
Finally we set up a set of pseudopotentials. For our intention, we install them locally.
First download and unpack your favorite set of pseudopotentials, then configure a `PseudoSet`.
Behind the scenes it will read through the directory you specified and resolve element names for the files in that directory that have the `.UPF` extension.

# Tutorial

In order to begin a ground state search, we need a template `scf` file that will be used to fill out most of the input parameters.
An example for `MnO` can be found in the `test/assets/` subdirectory.

## Creation
Once you have such a file for your case (or if you just want to try with `MnO`), execute:
```
cd ~/.julia/dev/RomeoDFT
gomeo searcher create MnO_tutorial test/assets/MnO_scf.in --sleep-time 5
```
the `sleep-time` flag is used to specify how often the orchestrator polls for new results etc.

In this case, the structure definition will also be taken from the template file.
If you want to use a different `pw.x` input or cif file from which to extract the structure, you can use
```
gomeo searcher create MnO_tutorial test/assets/MnO_scf.in -s <structure_file> --sleep-time 5
```
There will be an interactive menu to further set up the `searcher`, after which it will be started by the Orchestrator.

For further information on the `searcher create` command you can run
```
gomeo searcher create -h
```

## Monitoring
To get an overview of all the loaded `Searchers` you can execute
```
gomeo orchestrator status
```
To view information of a specific `Searcher`, execute:
```
gomeo searcher status <name>
```
where `<name>` in our case would be `MnO_tutorial`.

## Stopping/Starting

To stop a running `Searcher` you can execute:
```
gomeo searcher stop <name>
```
and to restart it
```
gomeo searcher start <name>
```
If no `Searcher` with `name` was loaded yet, it will be first.

## Changing mode
By default a `Searcher` will start in the `search` mode.
This means that it is actively exploring the energy landscape in order to find unique metastable states.
It will continue to do so until the amount of unique states versus trials goes below the value set with the `unique-ratio` flag during `gomeo searcher create`.
Afterwards, it will switch to the `postprocess` mode, meaning that it will not generate new trials, but simply process all the remaining trials and finish all the pending calculations.
If the ratio of unique/trial happens to exceed the `unique-ratio` value after finishing the postprocessing, the mode will switch back to `search`.

Sometimes, however, it is desirable to manually switch the mode in order to limit the search or simply stop a searcher from exploring further but cleaning up everything.
This can be done using:
```
gomeo searcher mode <name> <mode>
```
where `<mode>` is one of the following:
- `search`: Active searching through generating new trials.
- `postprocess`: No new trials are generated, but all remaining ones are submitted and cleaned up.
- `cleanup`: No new trials are generated, no new jobs are submitted (i.e. previously generated but not yet submitted trials will not be submitted), and searcher finishes after all currently running calculations finish.

## Plot States
An overview of the different states that are found so far can be plotted using
```
gomeo searcher plot <name> <outfile>
```
which will produce an image and save it to `outfile` displaying the states.

## Other commands
- `gomeo searcher unload <name>`: Stops, saves and removes a `Searcher` from the `Orchestrator`
- `gomeo searcher load <name>`: Loads a `Searcher`.
