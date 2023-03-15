using Pkg
using TOML
using UUIDs
const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "RomeoDFT") :
                   abspath(Base.DEPOT_PATH[1], "config", "RomeoDFT")

config_path(path...) = joinpath(CONFIG_DIR, path...)
mkpath(config_path())

paths = ["config.toml"]
if !ispath(config_path(paths[1]))
    conf = Dict("searchers_directory" => config_path("searchers"))
    open(config_path(paths[1]), "w") do f
        TOML.print(f, conf)
    end
end

if !haskey(ENV, "CI")
    Pkg.activate(CONFIG_DIR)
    Pkg.add(["Plots","LaTeXStrings","UnicodePlots","Revise"])
    if !haskey(Pkg.dependencies(), UUIDs.UUID("87c4fabc-abb4-4467-86a6-1748b5c259fe"))
        Pkg.add(url="https://github.com/louisponet/RomeoDFT.jl")
    end
    Pkg.update()

    if !occursin("config/RomeoDFT", pwd())
        Pkg.activate(".")
        using RomeoDFT
        RomeoDFT.comonicon_install()
    end
end
using RemoteHPC
if !RemoteHPC.exists(Environment("default"))
    save(Environment(name="default"))
end
