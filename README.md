# RomeoDFT

# Installation

Here we go through the steps required to install all the parts and software that are needed to perform global searches with Quantum Espresso and `romeo`.

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
git clone git@github.com:louisponet/RomeoDFT.jl RomeoDFT
```
Press backspace and `]` which activates the `pkg>` REPL mode, which allows you to install packages.
Execute
```
up
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

## romeo executable
After having performed the above steps, there should be the `romeo` executable in `~/.julia/bin`.
You can use it as, e.g. `~/.julia/bin/romeo -h`, or add `~/.julia/bin` to your `PATH`, after which you can use it from anywhere simply as `romeo -h`.

## RemoteHPC Configuration
Next we configure the `Servers`, `Execs` and `Environments` that represent where, what and how things can be ran, in our case the `qe-occupations/bin/pw.x` executable.

From a normal `shell`, execute
```
romeo configure
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
