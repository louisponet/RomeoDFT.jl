using Pkg
using TOML
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
    Pkg.add(url="https://github.com/louisponet/RomeoDFT.jl")
    Pkg.update()


    Pkg.activate(".")
    using RomeoDFT
    RomeoDFT.comonicon_install()
end
using RemoteHPC
if !RemoteHPC.exists(Environment("default"))
    save(Environment(name="default"))
end
