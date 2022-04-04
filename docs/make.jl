using Interesso
using Documenter

DocMeta.setdocmeta!(Interesso, :DocTestSetup, :(using Interesso); recursive=true)

makedocs(;
    modules=[Interesso],
    authors="astroEduardo <72969764+astroEduardo@users.noreply.github.com> and contributors",
    repo="https://github.com/JuDO-dev/Interesso.jl/blob/{commit}{path}#{line}",
    sitename="Interesso.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://JuDO-dev.github.io/Interesso.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuDO-dev/Interesso.jl",
    devbranch="dev",
)
