using StateSpaceIFO
using Documenter

DocMeta.setdocmeta!(StateSpaceIFO, :DocTestSetup, :(using StateSpaceIFO); recursive=true)

makedocs(;
    modules=[StateSpaceIFO],
    authors="Jeffrey Wack <52295204+jeffwack@users.noreply.github.com> and contributors",
    sitename="StateSpaceIFO.jl",
    format=Documenter.HTML(;
        canonical="https://jeff-wack-bot.github.io/StateSpaceIFO.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/jeff-wack-bot/StateSpaceIFO.jl",
    devbranch="main",
)
