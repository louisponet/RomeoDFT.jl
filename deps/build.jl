using Pkg
using TOML
const CONFIG_DIR = occursin("cache", first(Base.DEPOT_PATH)) ?
                   abspath(Base.DEPOT_PATH[2], "config", "RomeoDFT") :
                   abspath(Base.DEPOT_PATH[1], "config", "RomeoDFT")

config_path(path...) = joinpath(CONFIG_DIR, path...)
mkpath(config_path())

paths = ["config.toml"]

Pkg.activate(CONFIG_DIR)
Pkg.add(["Plots","LaTeXStrings","UnicodePlots","Revise"])
Pkg.add(url="https://github.com/louisponet/RomeoDFT.jl")
Pkg.update()

conf = Dict("searchers_directory" => config_path("searchers"))
open(config_path(paths[1]), "w") do f
    TOML.print(f, conf)
end

using RemoteHPC
if !RemoteHPC.exists(local_server(), Environment("default"))
    save(local_server(), Environment(name="default"))
end
Pkg.activate(".")
using RomeoDFT
RomeoDFT.comonicon_install()
