using Documenter
using RomeoDFT
using DocThemeIndigo

indigo = DocThemeIndigo.install(Configurations)

makedocs(;
    modules = [RomeoDFT],
    format = Documenter.HTML(
        prettyurls = !("local" in ARGS),
        canonical="https://Louis Ponet.github.io/RomeoDFT.jl",
        assets=String[indigo],
    ),
    pages = [
        "Home" => "index.md",
    ],
    repo = "https://github.com/Louis Ponet/RomeoDFT.jl",
    sitename = "RomeoDFT.jl",
)

deploydocs(; repo = "https://github.com/Louis Ponet/RomeoDFT.jl")
