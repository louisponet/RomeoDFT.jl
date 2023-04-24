# Tutorial

In order to begin a ground state search, we need a template `scf` file that will be used to fill out most of the input parameters.
An example for `MnO` can be found in the `test/assets/` subdirectory.

## Creation
Once you have such a file for your case (or if you just want to try with `MnO`), execute:
```
cd ~/.julia/dev/RomeoDFT
romeo searcher create MnO_tutorial test/assets/MnO_scf.in --sleep-time 5
```
the `sleep-time` flag is used to specify how often the orchestrator polls for new results etc.

In this case, the structure definition will also be taken from the template file.
If you want to use a different `pw.x` input or cif file from which to extract the structure, you can use
```
romeo searcher create MnO_tutorial test/assets/MnO_scf.in -s <structure_file> --sleep-time 5
```
There will be an interactive menu to further set up the [`Searcher`](@ref), after which it will be started by the [Orchestrator](@ref).

For further information on the `searcher create` command you can run
```
romeo searcher create -h
```

## Monitoring
To get an overview of all the loaded [Searchers](@ref) you can execute
```
romeo orchestrator status
```
To view information of a specific [`Searcher`](@ref), execute:
```
romeo searcher status <name>
```
where `<name>` in our case would be `MnO_tutorial`.

## Stopping/Starting

To stop a running [`Searcher`](@ref) you can execute:
```
romeo searcher stop <name>
```
and to restart it
```
romeo searcher start <name>
```
If no [`Searcher`](@ref) with `name` was loaded yet, it will be first.

## Changing mode
By default a [`Searcher`](@ref) will start in the `search` mode.
This means that it is actively exploring the energy landscape in order to find unique metastable states.
It will continue to do so until the amount of unique states versus trials goes below the value set with the `unique-ratio` flag during `romeo searcher create`.
Afterwards, it will switch to the `postprocess` mode, meaning that it will not generate new trials, but simply process all the remaining trials and finish all the pending calculations.
If the ratio of unique/trial happens to exceed the `unique-ratio` value after finishing the postprocessing, the mode will switch back to `search`.

Sometimes, however, it is desirable to manually switch the mode in order to limit the search or simply stop a searcher from exploring further but cleaning up everything.
This can be done using:
```
romeo searcher mode <name> <mode>
```

See [Modes](@ref) for more info.

## Plot States
An overview of the different states that are found so far can be plotted using
```
romeo searcher plot <name> <outfile>
```
which will produce an image and save it to `outfile` displaying the states.

## Other commands
See [CommandLine](@ref) for more commands.
