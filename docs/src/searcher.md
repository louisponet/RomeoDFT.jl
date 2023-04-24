# Searchers
```@meta
CurrentModule = RomeoDFT
```

A [`Searcher`](@ref) is the main object representing all the data related with a global search. Most interactions can be performed from the command line using the `romeo searcher` subset of commands.
Each [`Searcher`](@ref) is assigned a label or `name` with which it can be identified.

## Data
The `name` also references the sub directory where all the data will be stored, which can be retrieved using `dirname`.
This data consists of:
- `$(DATABASE_VERSION)/ledger.jld2`: stores all the search data that is kept in memory during the search. See [`load`](@ref) and [`save`](@ref)
- `job_backups`: a directory storing all the input and output data of the various calculations ran during the search
- `log.log`: file that keeps a log of everything that happens during the search. The amount of information logged can be tuned using different levels passed to `--verbosity <n>` flag of `romeo searcher create` or `romeo searcher start`
- `HTTP.log`: file holding all HTTP related logs
- `report.out`: a report generated when a [`Searcher`](@ref) is finished

## CommandLine
Most interaction with a [`Searcher`](@ref) can be done through the command line using the `romeo searcher` subcommands (use `-h` to show more info on commands and flags).
The most important ones are:
- `romeo searcher create`: interactive wizard to initialize and start a new [`Searcher`](@ref), see [`setup_search`](@ref)
- `romeo searcher start`: starts a [`Searcher`](@ref)
- `romeo searcher stop`: stops a [`Searcher`](@ref)
- `romeo searcher status`: shows a status overview
- `romeo searcher unload`: unloads a [`Searcher`](@ref) from the orchestrator
- `romeo searcher log`: opens the `log.log` file in the editor specified by the `\$EDITOR` environment variable
- `romeo searcher plot`: creates a standard [`State`](@ref) plot
- `romeo searcher mode`: sets the mode, see [Modes](@ref) for more information

## Life Cycle
The lifecycle of a [`Searcher`](@ref) can be summarized by the following steps:
1. Potential geometry and U parameter optimization when `--hp-base` or `--relax-base` are passed to `romeo searcher create`
2. Initial unconstrained `scf` calculation to determine the number of electrons in the target manifold.
3. Random search with `--nrand` trials
4. Based on the unique [`States`](@ref State) found in `3.` (determined by `--unique-thr`) perform midpoint intersections (`n = (n1 + n2)/2`) between each pair of unique states
5. repeat `4.` for each newly found unique state
6. if no new intersections can be created and the stopcondition is not yet reached go to step `3.`
7. if the stop condition determined by `--stopping-unique-ratio` and `--stopping-n-generations` is reached switch to `postprocess` mode
8. if all jobs are finished and the stop condition is still satisfied, stop the [`Searcher`](@ref) and print results to `report.out`, otherwise switch back to `3.`

## Modes
There several execution modes that a [`Searcher`](@ref) can be in:
- `search`: search mode characterized by performing intersection and random search, and postprocessing on unique [`States`](@ref State).
- `postprocess`: when the stopcondition is satisfied, no new intersections or random trials will be generated, and only `postprocessing` actions will be taken
- `cleanup`: no new trials or postprocessing jobs will be created and submitted, only cleanup of the currently running jobs will be performed, after which the [`Searcher`](@ref) will finish
- `manual`: only used when performing manual tasks and you don't want any searching or postprocessing to happen automatically

These can be set using `romeo searcher mode`. 

## Reference
```@docs
Searcher
setup_search
load(::Searcher)
save(::Searcher)
sssp_distance
```
