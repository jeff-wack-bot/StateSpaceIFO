# H∞ vs H² filter-cavity design on the KLMTV plant.
#
# Run (about a minute):  julia --project=examples examples/klmtv_hinf.jl
#
# Setup (γ = G = m = 1, I₀/I_SQL = 1, so 𝒦 = 2/(Ω²(1+Ω²)) ):
#   readout coordinate      q(Ω) = cot(θ − α_rot(Ω)),   α_rot = Σ arctan(ξ ± Ω/δ)
#   excess over the ideal   R(Ω) = S_h(Ω)/S_h^min(Ω) = (q(Ω) − 𝒦(Ω))² + 1
# so, measured against the ideal variational curve,
#   H²  design  =  weighted L²  approximation of 𝒦 by q
#   H∞  design  =  unweighted L∞ (Chebyshev) approximation of 𝒦 by q
# with NO weight left over in the H∞ case — R is already dimensionless.

using StateSpaceIFO
using Printf

# ---------------- objectives ----------------
band(lo, hi, n) = exp.(range(log(lo), log(hi); length = n))

"Inspiral-weighted integrated excess: ∫W·S_h / ∫W·S_h^min  (a pure ratio)."
function J2(p, Ωs, wts)
    num = 0.0; den = 0.0
    for (i, Ω) in enumerate(Ωs)
        r = excess_ratio(p, Ω)
        !isfinite(r) && return Inf
        num += wts[i] * r; den += wts[i]
    end
    return num / den
end

"Worst-case excess: max_Ω R(Ω)."
function Jinf(p, Ωs)
    m = 0.0
    for Ω in Ωs
        r = excess_ratio(p, Ω)
        !isfinite(r) && return Inf
        r > m && (m = r)
    end
    return m
end

# ---------------- global search: coarse sweep + multistart local refine ----------------
# The objective is multi-modal in (ξ, δ, θ) — see klmtv_state-space.jl section D — so a
# single shrinking-grid descent lands in whichever basin it starts in.  Sweep first,
# then refine the best several basins independently and keep the winner.

grid1(lo, hi, n, lg) = lg ? exp.(range(log(lo), log(hi); length = n)) :
    collect(range(lo, hi; length = n))

"Local shrinking-box descent from a starting box.  Slow shrink (×2/3 per round) so
that the box does not collapse before the basin is located."
function descend(f, lo, hi, logscale; rounds = 22, npts = 9)
    lo = copy(lo); hi = copy(hi); best = (Inf, copy(lo))
    for _ in 1:rounds
        ax = [grid1(lo[i], hi[i], npts, logscale[i]) for i in eachindex(lo)]
        for idx in Iterators.product(ax...)
            x = collect(idx); v = f(x)
            v < best[1] && (best = (v, x))
        end
        x0 = best[2]
        for i in eachindex(lo)
            if logscale[i]
                w = (log(hi[i]) - log(lo[i])) / 3
                lo[i], hi[i] = exp(log(x0[i]) - w), exp(log(x0[i]) + w)
            else
                w = (hi[i] - lo[i]) / 3
                lo[i], hi[i] = x0[i] - w, x0[i] + w
            end
        end
    end
    return best
end

function optimize(f, lo, hi, logscale; coarse = 25, nkeep = 20, halo = 4, seeds = [])
    # (1) a coarse-to-fine descent over the WHOLE box,
    # (2) descents seeded from the best coarse-grid points, each with a generous halo,
    # (3) descents seeded from supplied points (e.g. the H² optimum warm-starting H∞).
    # Keep the winner.  (1) alone gets trapped; (2) alone can cut off the true basin.
    best = descend(f, lo, hi, logscale; rounds = 22, npts = 15)
    ax = [grid1(lo[i], hi[i], coarse, logscale[i]) for i in eachindex(lo)]
    pts = Tuple{Float64, Vector{Float64}}[]
    for idx in Iterators.product(ax...)
        x = collect(idx); v = f(x)
        isfinite(v) && push!(pts, (v, x))
    end
    sort!(pts; by = t -> t[1])
    for x0 in vcat([p[2] for p in pts[1:min(end, nkeep)]], collect(seeds))
        l = similar(x0); h = similar(x0)
        for i in eachindex(x0)
            if logscale[i]
                w = halo * (log(hi[i]) - log(lo[i])) / (coarse - 1)
                l[i] = max(lo[i], exp(log(x0[i]) - w)); h[i] = min(hi[i], exp(log(x0[i]) + w))
            else
                w = halo * (hi[i] - lo[i]) / (coarse - 1)
                l[i] = max(lo[i], x0[i] - w); h[i] = min(hi[i], x0[i] + w)
            end
        end
        v, x = descend(f, l, h, logscale)
        v < best[1] && (best = (v, x))
    end
    return best
end

println("="^74)
println("0.  Sanity: TWO cavities reach the pointwise optimum exactly (R ≡ 1),")
println("    so H² and H∞ agree there and both reduce to pole placement.")
println("="^74)
cav = klmtv_cavities(1.0)
p2 = [cav[1].ξ, cav[1].δ, cav[2].ξ, cav[2].δ, π / 2]
@printf(
    "  KLMTV parameters:  ξ_I = %+.6f  δ_I = %.6f   ξ_II = %+.6f  δ_II = %.6f  θ = π/2\n",
    p2[1], p2[2], p2[3], p2[4]
)
Ωchk = band(0.02, 50.0, 400)
@printf("  max |R(Ω) − 1| over Ω/γ ∈ [0.02, 50] : %.3e\n", maximum(abs(excess_ratio(p2, Ω) - 1) for Ω in Ωchk))
println("  ⟹ at n = deg 𝒦 the two design norms cannot disagree. The H²/H∞ question")
println("     only has content BELOW the exact order.")

# ---------------- the interesting case: ONE cavity ----------------
for (lo, hi) in [(0.3, 3.0), (0.1, 10.0)]
    Ωs = band(lo, hi, 801)
    wts = Ωs .^ (-7 / 3) .* [
        i == 1 ? Ωs[2] - Ωs[1] : i == length(Ωs) ? Ωs[end] - Ωs[end - 1] :
            (Ωs[i + 1] - Ωs[i - 1]) / 2 for i in eachindex(Ωs)
    ]
    println()
    println("="^74)
    @printf("1.  ONE filter cavity, band Ω/γ ∈ [%.1f, %.1f]\n", lo, hi)
    println("="^74)

    lo3 = [-8.0, 0.02, 0.02π]; hi3 = [8.0, 50.0, 0.98π]; ls3 = [false, true, false]
    v2, x2 = optimize(p -> J2(p, Ωs, wts), lo3, hi3, ls3)
    vinf, xi = optimize(p -> Jinf(p, Ωs), lo3, hi3, ls3; seeds = [x2])
    # self-check: each design must win on its own objective, or the search failed
    (Jinf(xi, Ωs) <= Jinf(x2, Ωs) + 1.0e-9) ||
        println("  !! search failure: H∞ design is worse in H∞ than the H² design")
    (J2(x2, Ωs, wts) <= J2(xi, Ωs, wts) + 1.0e-9) ||
        println("  !! search failure: H² design is worse in H² than the H∞ design")

    @printf(
        "  %-22s %9s %9s %9s   %10s %10s\n",
        "design", "ξ", "δ/γ", "θ/π", "∫W (dB)", "worst (dB)"
    )
    for (nm, x) in [("H²  (weighted L²)", x2), ("H∞  (Chebyshev)", xi)]
        @printf(
            "  %-22s %+9.4f %9.4f %9.4f   %10.4f %10.4f\n",
            nm, x[1], x[2], x[3] / π, 10log10(J2(x, Ωs, wts)), 10log10(Jinf(x, Ωs))
        )
    end

    # equioscillation check for the H∞ design
    Rs = [excess_ratio(xi, Ω) for Ω in Ωs]
    peaks = [(Ωs[i], Rs[i]) for i in 2:(length(Ωs) - 1) if Rs[i] >= Rs[i - 1] && Rs[i] >= Rs[i + 1]]
    push!(peaks, (Ωs[1], Rs[1])); push!(peaks, (Ωs[end], Rs[end]))
    sort!(peaks; by = t -> -t[2])
    @printf("\n  H∞ design: interior/edge maxima of R(Ω), largest first (equioscillation)\n")
    for (Ω, r) in peaks[1:min(end, 5)]
        @printf(
            "      Ω/γ = %7.3f   R = %8.5f   (%+7.4f dB above the ideal curve)\n",
            Ω, r, 10log10(r)
        )
    end
    nflat = count(t -> t[2] > 0.99 * peaks[1][2], peaks)
    @printf("  maxima within 1%% of the peak: %d  (equiripple ⟹ ≥ 2)\n", nflat)

    # where each design is worst, and how they trade
    @printf("\n  cross-comparison (dB relative to the ideal variational curve):\n")
    @printf(
        "      H² design : ∫W = %+7.4f   worst = %+7.4f\n",
        10log10(J2(x2, Ωs, wts)), 10log10(Jinf(x2, Ωs))
    )
    @printf(
        "      H∞ design : ∫W = %+7.4f   worst = %+7.4f\n",
        10log10(J2(xi, Ωs, wts)), 10log10(Jinf(xi, Ωs))
    )
    @printf(
        "      H∞ buys %.4f dB of worst case for %.4f dB of band-average.\n",
        10log10(Jinf(x2, Ωs)) - 10log10(Jinf(xi, Ωs)),
        10log10(J2(xi, Ωs, wts)) - 10log10(J2(x2, Ωs, wts))
    )
end

println()
println("="^74)
println("2.  Feasibility reading: γ* = min max_Ω S_h/S_target answers 'can this")
println("    design goal be met?' directly.  Target = c × (ideal curve):")
println("="^74)
Ωs = band(0.3, 3.0, 801)
lo3 = [-8.0, 0.02, 0.02π]; hi3 = [8.0, 50.0, 0.98π]; ls3 = [false, true, false]
_, xw = optimize(p -> J2(p, Ωs, Ωs .^ (-7 / 3)), lo3, hi3, ls3)
vinf, xi = optimize(p -> Jinf(p, Ωs), lo3, hi3, ls3; seeds = [xw])
@printf(
    "  best worst-case excess with 1 cavity on Ω/γ ∈ [0.3,3] : %.4f  (%+.4f dB)\n",
    vinf, 10log10(vinf)
)
for c in [1.0, 1.2, 1.5, 2.0]
    @printf(
        "    goal S_h ≤ %.1f × ideal  ⟹  γ* = %.4f  %s\n",
        c, vinf / c, vinf / c <= 1 ? "FEASIBLE with 1 cavity" : "INFEASIBLE — needs cavity #2"
    )
end
println("\ndone.")
