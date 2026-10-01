using StateSpaceIFO
using Documenter

DocMeta.setdocmeta!(StateSpaceIFO, :DocTestSetup, :(using StateSpaceIFO); recursive = true)

makedocs(;
    modules = [StateSpaceIFO],
    authors = "Jeffrey Wack <52295204+jeffwack@users.noreply.github.com> and contributors",
    sitename = "StateSpaceIFO.jl",
    format = Documenter.HTML(;
        canonical = "https://jeff-wack-bot.github.io/StateSpaceIFO",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Conventions" => "conventions.md",
        "KLMTV in state space" => "klmtv_state_space.md",
        "Filter cavities as Blaschke factors" => "pole_placement.md",
        "Readout design" => "readout_design.md",
        "Chain scattering" => "chain_scattering.md",
        "Factoring the KLMTV plant" => "dual_factorization.md",
        "Quadratic invariance" => "quadratic_invariance.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo = "github.com/jeff-wack-bot/StateSpaceIFO",
    devbranch = "main",
)
