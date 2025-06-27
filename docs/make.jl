using Interesso
using Documenter
using DocumenterInterLinks

DocMeta.setdocmeta!(Interesso, :DocTestSetup, :(using Interesso); recursive=true)

const _PAGES = [
    "Home" => "index.md",
    "Examples" => "examples.md",
    "API Reference" => [
        "reference/interpolants.md",
        "reference/points.md",
        "reference/methods.md",
        "reference/bounds.md",
        "reference/intervals.md",
    ],
    "Changelog" => "changelog.md",
]

const _LINKS = InterLinks(
    "DOI" => "https://judo.dev/DynOptInterface.jl/dev/objects.inv"
)

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
    pages=_PAGES,
    plugins=[_LINKS],
)

deploydocs(;
    repo="github.com/JuDO-dev/Interesso.jl",
    devbranch="dev",
)