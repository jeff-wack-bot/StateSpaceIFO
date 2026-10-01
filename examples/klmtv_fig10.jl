# Reproduce KLMTV Fig. 10: parameters of the two variational-readout filter
# cavities vs (Λ/γ)⁴ = 2·I_o/I_SQL.
#
# The poles of the output all-pass 𝔄(s) = (D + iΛ⁴)/(D − iΛ⁴), D = s² + γ²s,
# in units γ = 1 and s = Ω², are
#
#     p_{I,II} = (1/2)·(−1 ± √(1 + 8i·I_o/I_SQL))
#
# and each pole is one cavity: ω = √p with Im ω < 0, half-bandwidth δ = −Im ω,
# offset ξ = −Re ω / δ.
#
# Run:  julia --project=examples examples/klmtv_fig10.jl

using StateSpaceIFO
using CairoMakie

x = 10 .^ range(-2, 2; length = 400)     # (Λ/γ)⁴ = 2·I_o/I_SQL, KLMTV's x-axis
I0 = x ./ 2                              # I_o/I_SQL

cavs = klmtv_cavities.(I0)
ξI = [c[1].ξ for c in cavs]
δI = [c[1].δ for c in cavs]
ξII = [c[2].ξ for c in cavs]
δII = [c[2].δ for c in cavs]

fig = Figure(size = (520, 420))
ax = Axis(
    fig[1, 1];
    xscale = log10,
    xlabel = L"(\Lambda/\gamma)^4 = 2 I_o / I_\mathrm{SQL}",
    ylabel = "filter parameters"
)
lines!(ax, x, ξI; label = L"\xi_\mathrm{I}")
lines!(ax, x, δI; label = L"\delta_\mathrm{I}/\gamma")
lines!(ax, x, ξII; label = L"\xi_\mathrm{II}")
lines!(ax, x, δII; label = L"\delta_\mathrm{II}/\gamma")
axislegend(ax; position = :lt)

save(joinpath(@__DIR__, "klmtv_fig10.png"), fig)
println("wrote klmtv_fig10.png")
