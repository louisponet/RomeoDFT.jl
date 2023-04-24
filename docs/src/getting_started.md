# Installation

Here we go through the steps required to install all the parts and software that are needed to perform global searches with Quantum Espresso and ROMEO.
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

Begin by starting a julia REPL session, and execute:
```julia
using Pkg
Pkg.add("RomeoDFT")
```
or press the `]` key to enter the julia `pkg` REPL mode and type `add RomeoDFT`.
The package is now being intialized, sit back grab a coffee, this takes a long time (~5-10 min).

To update to a newer version you can do
```julia
using Pkg
Pkg.update("RomeoDFT")
```
or `up RomeoDFT` in the `pkg` REPL mode.

## ROMEO executable
After having performed the above steps, there should be the `romeo` executable in `~/.julia/bin`.
You can use it as, e.g. `~/.julia/bin/romeo -h`.

!!! note
    It is strongly recommended to add ~/.julia/bin to your PATH environment variable by changing `~/.bashrc` or similar.

## RemoteHPC Configuration
Next we configure the `Servers`, `Execs` and `Environments` that represent where, what and how executables can be ran on the remote clusters.
In our case the `qe-occupations/bin/pw.x` executable.

!!! note
    We will be using the default editor to write the configuration files in the following steps. If the `EDITOR` environment variable is not set, this
    will be `nano`. Again it is advised to set it with `export EDITOR=<vim vscode or any other editor>` or in your `~/.bashrc`.

From inside a shell, execute
```
romeo configure
```

### Server
Use the arrow keys and enter to select `Server`.
A `Server` represents a computer or cluster where executables can be ran.
Your local computer is configured automatically so if you want to run everything locally, you can skip this step.
Give the server a `name` that is used to load it later, and insert your `username` and `domain` that correspond to the `ssh` commands
you use to connect to the cluster where you installed `qe-occupations` (e.g. `ssh lponet@cscs.ch` would mean `username = lponet`, `domain = cscs.ch`).

If you plan to run on multiple clusters, you can add them one after the other. For now, however, each [`Searcher`](@ref) runs on a single `Server`.

If you configured a new `Server`, go out of the configuration at this point and run `romeo server start <name of new server>` then run `romeo configure` again.

### Exec
Next select `Exec` from the menu. Read carefully the documentation that is printed, enter a `name` (for example `qe-occupations`), select the `Server` that corresponds to where you installed the patched `QuantumESPRESSO`.
After pressing enter, an editor will open with a representation of the `Exec` with some fields already filled in.
Complete the remaining fields, especially `path` which should be the **absolute** directory of the `qe-occupations/bin/pw.x` executable on the remote `Server`.
Again, read the Documentation that is printed below the fields to help filling out the remaining information.

### Environment
An environment essentially represents the #SBATCH commands of a slurm job script, or similar.
Carefully read the documentation and fill out the fields similar to above.

### Pseudo Potentials
Finally we set up a set of pseudopotentials. For our intention, we install them locally.
First download and unpack your favorite set of pseudopotentials into a local directoyr, then supply the directory during the configuration of a `PseudoSet`.
Behind the scenes it will read through the directory you specified and resolve element names for the files in that directory that have the `.UPF` extension.

You're now ready to start the [Orchestrator](@ref) and run your first [`Searcher`](@ref).

