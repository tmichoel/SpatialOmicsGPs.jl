push!(LOAD_PATH,"../src/")

#using SpatialOmicsGPs
using Documenter, SpatialOmicsGPs

#DocMeta.setdocmeta!(SpatialOmicsGPs, :DocTestSetup, :(using SpatialOmicsGPs); recursive=true)

makedocs(;
    modules=[SpatialOmicsGPs],
    sitename="SpatialOmics.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://tmichoel.github.io/SpatialOmicsGPs.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Introduction" => "index.md",
        "FaST-LMM" => "FaST-LMM.md",
        "List of functions" => "listfunctions.md"
    ],
)

deploydocs(;
    repo="github.com/tmichoel/SpatialOmicsGPs.jl.git",
    devbranch="main",
)
