push!(LOAD_PATH,"../src/")

#using SpatialOmicsGPs
using Documenter, SpatialOmicsGPs

#DocMeta.setdocmeta!(SpatialOmicsGPs, :DocTestSetup, :(using SpatialOmicsGPs); recursive=true)

makedocs(;
    modules=[SpatialOmicsGPs],
    authors="tmichoel <11647967+tmichoel@users.noreply.github.com> and contributors",
    repo="https://github.com/tmichoel/SpatialOmics.jl/blob/{commit}{path}#{line}",
    sitename="SpatialOmics.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tmichoel.github.io/SpatialOmics.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "FaST-LMM" => "FaSTLMM.md",
        "List of functions" => "listfunctions.md"
    ],
)

deploydocs(;
    repo="github.com/tmichoel/SpatialOmics.jl",
    devbranch="main",
)
