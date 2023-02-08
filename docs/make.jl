using AdditiveCellCom
using Documenter

makedocs(;
    modules=[AdditiveCellCom],
    authors="Matt Karikomi <mattkarikomi@gmail.com> and contributors",
    repo="https://github.com/mkarikom/AdditiveCellCom.jl/blob/{commit}{path}#L{line}",
    sitename="AdditiveCellCom.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mkarikom.github.io/AdditiveCellCom.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mkarikom/AdditiveCellCom.jl",
)
