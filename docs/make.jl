using FETI
using Documenter

DocMeta.setdocmeta!(FETI, :DocTestSetup, :(using FETI); recursive=true)

makedocs(;
    modules=[FETI],
    authors="Jared Frazier <cscidev001@gmail.com> and contributors",
    sitename="FETI.jl",
    format=Documenter.HTML(;
        canonical="https://jfdev001.github.io/FETI.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="https://github.com/jfdev001/FETI.jl",
)
