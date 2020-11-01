using GCom
using Documenter

makedocs(;
    modules=[GCom],
    authors="Matt Karikomi <mattkarikomi@gmail.com> and contributors",
    repo="https://github.com/mkarikom/GCom.jl/blob/{commit}{path}#L{line}",
    sitename="GCom.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://mkarikom.github.io/GCom.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/mkarikom/GCom.jl",
)
