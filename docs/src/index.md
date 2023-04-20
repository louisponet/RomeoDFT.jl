# RomeoDFT.jl

This software package allows one to perform a **R**obust **O**ccupation **M**atrix **E**nergy **O**ptimization of a crystalline structure using **D**ensity **F**unctional **T**heory.

It is most useful in materials with complex magnetic configurations or where the local orbital physics of one or more atomic species plays a crucial role in describing the ground state.
In many of these cases the energy surface of the DFT functional sports many _local minima_, which in theory should be extensively explored in order to find the _global minimum_ associated with the _ground state_.

`RomeoDFT.jl` facilitates a fully automated global search, exploring many different occupations of local atomic orbitals (`target states`) in order to reliably determine the ground state without much human interaction.

It uses a patched version of [QuantumESPRESSO](https://gitlab.com/louisponet/q-e) to run the DFT calculations for each of the target states.

## Features
- "Press of the button" determination of the ground state starting from a single template input file with the structure and input parameters
- Seamless interfacing with HPC clusters using [RemoteHPC.jl](https://github.com/louisponet/RemoteHPC.jl) and [DFControl.jl](https://github.com/louisponet/DFControl.jl)
- Run geometry optimization for the initial structure as well as any found self-consistent state associated with a local minimum
- Automatically and fully self-consistently determine the Hubbard U parameter for the initial structure as well as any self-consistent state associated with a local minimum
- PostProcessing steps including non-selfconsistent, bandstructure and dos/projected dos calculations
- Robust data storage using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) with versioning
- Command line interface facilitating most of the tasks
- Built with ECS (see [Overseer.jl](https://github.com/louisponet/Overseer.jl) making it very extensible


## Acknowledgements
This work has been funded by H2020-OpenModel and EPFL.
