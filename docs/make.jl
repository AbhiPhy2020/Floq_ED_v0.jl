using Floq_ED_v0
using Documenter

DocMeta.setdocmeta!(Floq_ED_v0, :DocTestSetup, :(using Floq_ED_v0); recursive=true)

makedocs(;
    modules=[Floq_ED_v0],
    authors="Abhishek Kumar",
    sitename="Floq_ED_v0.jl",
    format=Documenter.HTML(;
        canonical="https://abhiphys2020.github.io/Floq_ED_v0.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/abhiphys2020/Floq_ED_v0.jl",
    devbranch="main",
)
