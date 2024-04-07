push!(LOAD_PATH,"../src/")

#using SpatialOmicsGPs
using Documenter, SpatialOmicsGPs

#DocMeta.setdocmeta!(SpatialOmicsGPs, :DocTestSetup, :(using SpatialOmicsGPs); recursive=true)

makedocs(;
    modules=[SpatialOmicsGPs],
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://tmichoel.github.io/SpatialOmics.jl",
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
    repo="github.com/tmichoel/SpatialOmics.jl.git",
    devbranch="main",
)
