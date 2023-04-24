# Orchestrator
```@meta
CurrentModule = RomeoDFT
```

The [Orchestrator](@ref) can be thought of as a locally running daemon that handles running the various [Searchers](@ref).
There are three main commands that are useful to interact with it (from a standard `shell`):
- `romeo orchestrator start`: starts the [Orchestrator](@ref) and loads and starts the searchers that were previously running
- `romeo orchestrator stop`: stops and saves the currently loaded [Searchers](@ref) and gracefully stops the [Orchestrator](@ref)
- `romeo orchestrator status`: shows an overview of the currently loaded [Searchers](@ref)

For more information use `-h` after these commands in the shell.

## Advanced
Behind the scenes the [Orchestrator](@ref) runs in a [RemoteREPL](https://github.com/c42f/RemoteREPL.jl), which can be accessed to interact directly with
the loaded [`Searchers`](@ref Searcher).
This can be achieved by starting a julia REPL and executing the following commands:

```julia
using RomeoDFT
connect_orchestrator()
```

You will see that now by pressing `>` you get a pink julia prompt, which means you're executing commands inside the [Orchestrator](@ref) REPL.
The useful commands here are:

- [`start_searcher("<searcher_name>"; verbosity = <number>)`](@ref start_searcher): starts a [`Searcher`](@ref) with the provided name with a given logging `verbosity` (higher = more logged).
- [`stop_searcher("<searcher_name>")`](@ref stop_searcher): stops a [`Searcher`](@ref) with the name
- [`load_searcher`](@ref): simply loads a [`Searcher`](@ref) from disk if no function is provided as the first argument. See [`load_searcher documentation`](@ref load_searcher) for more info.

## Reference
```@docs
start_searcher
stop_searcher
load_searcher
```
