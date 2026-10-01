# KLMTV from a 4-state model, and the convexity of the readout-design problem.
#
# Run:  julia --project=examples examples/klmtv_state-space.jl
#
# The plant (γ = G = m = L = 1), states (c1, c2, q, p):
#   ċ1 = -γ c1 + √(2γ) a1                     (amplitude quadrature, undriven)
#   ċ2 = -γ c2 + √(2γ) a2 + G (q + L h)       (phase quadrature reads displacement)
#   q̇  = p/m
#   ṗ  = G c1                                 (radiation pressure ∝ amplitude)
#   b  = -a + √(2γ) c
#
# The two G's are equal *because* the coupling is Hamiltonian (H_int = -G c1 q);
# that equality is the physical-realizability condition, not a modelling choice.

using StateSpaceIFO
using LinearAlgebra, Printf

# ------------------------------------------------------------------ #
# A.  The state space reproduces KLMTV Eq. (16).                      #
# ------------------------------------------------------------------ #
println("="^72)
println("A.  State space  ⟹  KLMTV Eq. (16):  Gn = e^{-2iβ} M,  M = [[1,0],[-𝒦,1]]")
println("="^72)
P0 = klmtv_plant()
@printf("%10s %14s %14s %12s %12s\n", "Ω/γ", "𝒦 (model)", "𝒦 (analytic)", "|Gn[1,2]|", "arg Gn11/2")
maxerrK = 0.0; maxerrβ = 0.0; maxoff = 0.0; maxdiag = 0.0
for Ω in [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0]
    Ga, _, Gh = transfer_matrices(P0, Ω)
    Kmod = real(backaction_gain(Ga))
    Kan = klmtv_K(Ω)
    βmod = arm_phase(Ga)                 # s = +iΩ ⟹ our phase is e^{-2iβ}
    global maxerrK = max(maxerrK, abs(Kmod - Kan) / Kan)
    global maxerrβ = max(maxerrβ, abs(βmod - atan(Ω)))
    global maxoff = max(maxoff, abs(Ga[1, 2]))
    global maxdiag = max(maxdiag, abs(Ga[2, 2] - Ga[1, 1]))
    @printf("%10.3f %14.6e %14.6e %12.2e %12.6f\n", Ω, Kmod, Kan, abs(Ga[1, 2]), βmod)
end
@printf("\n  max rel. error in 𝒦      : %.3e\n", maxerrK)
@printf("  max abs. error in β      : %.3e   (β = arctan Ω/γ)\n", maxerrβ)
@printf("  max |Gn[1,2]|  (must be 0): %.3e\n", maxoff)
@printf("  max |Gn[2,2]-Gn[1,1]|    : %.3e\n", maxdiag)

# ------------------------------------------------------------------ #
# B.  Variational readout = inverse Fisher information.               #
#         S_h^min(Ω) = 1 / (Gh† Σ⁻¹ Gh),   Σ = Σ_ports G G†           #
#     and the optimal readout direction is  v ∝ Gh† Σ⁻¹  ∝ (𝒦, 1).    #
# ------------------------------------------------------------------ #
println()
println("="^72)
println("B.  S_h^min = 1/(Gh† Σ⁻¹ Gh)   vs   KLMTV  h²_SQL/2𝒦  (= 1/|Gh[2]|²)")
println("="^72)
@printf("%10s %16s %16s %12s %14s\n", "Ω/γ", "1/(Gh†Σ⁻¹Gh)", "1/|Gh[2]|²", "rel.err", "v1/v2  vs  𝒦")
maxerrS = 0.0; maxerrdir = 0.0
for Ω in [0.05, 0.2, 0.5, 1.0, 2.0, 5.0, 20.0]
    Ga, _, Gh = transfer_matrices(P0, Ω)
    Σ = Ga * Ga'
    Smin = sh_min(Gh, Σ)
    Spred = 1 / abs2(Gh[2])
    v = optimal_readout(Gh, Σ)
    ratio = real(v[1] / v[2])
    global maxerrS = max(maxerrS, abs(Smin - Spred) / Spred)
    global maxerrdir = max(maxerrdir, abs(ratio - klmtv_K(Ω)) / klmtv_K(Ω))
    @printf("%10.3f %16.8e %16.8e %12.2e %14.8f\n", Ω, Smin, Spred, abs(Smin - Spred) / Spred, ratio)
end
@printf("\n  max rel. error, S_h^min  : %.3e\n", maxerrS)
@printf("  max rel. error, v1/v2 = 𝒦: %.3e\n", maxerrdir)

# ------------------------------------------------------------------ #
# C.  Which combination of sideband phases is the rotation?           #
#     α± = 2 atan(ξ ± Ω/δ);  claim: rotation = (α₊+α₋)/2 (even in Ω)  #
#     and e^{i(α₊+α₋)} = (s - p̄)/(s - p),  s = Ω², p = (Δ - iδ)².     #
# ------------------------------------------------------------------ #
println()
println("="^72)
println("C.  Detuned cavity: which sideband combination is the quadrature rotation?")
println("="^72)
Δ, δ = 0.7, 0.4
ξ = -Δ / δ
p = (Δ - im * δ)^2
errsum = 0.0; errdif = 0.0
for Ω in [0.1, 0.35, 0.9, 1.7, 3.3]
    αp, αm = sideband_phases(ξ, δ, Ω)
    blaschke = (Ω^2 - conj(p)) / (Ω^2 - p)
    global errsum = max(errsum, abs(cis(αp + αm) - blaschke))
    global errdif = max(errdif, abs(cis(αp - αm) - blaschke))
end
@printf("  max |e^{i(α₊+α₋)} - (s-p̄)/(s-p)| = %.3e   <- SUM matches the Blaschke factor\n", errsum)
@printf("  max |e^{i(α₊-α₋)} - (s-p̄)/(s-p)| = %.3e   <- DIFFERENCE does not\n", errdif)
# tuned cavity sanity check: Δ = 0 must give zero quadrature rotation
rot0 = maximum(abs(2rotation_angle(0.0, 0.4, Ω)) for Ω in [0.1, 1.0, 5.0])
@printf("  tuned cavity (Δ=0): max |α₊+α₋| = %.3e  (rotation vanishes, as it must)\n", rot0)

# ------------------------------------------------------------------ #
# D.  Convexity in the readout coordinate vs. in the hardware.        #
#     S_h(Ω) = h²_SQL/(2𝒦) · [ (q - 𝒦)² + 1 ],   q = cot(readout ang) #
#     J(q) = ∫ W S_h dΩ  is a convex quadratic functional of q.       #
# ------------------------------------------------------------------ #
println()
println("="^72)
println("D.  J = ∫ W(Ω) S_h dΩ :  convex quadratic in the readout coordinate q")
println("="^72)

const Ωlo, Ωhi, NΩ = 0.1, 10.0, 2001
const Ωg = exp.(range(log(Ωlo), log(Ωhi); length = NΩ))
const dΩ = [i == 1 ? Ωg[2] - Ωg[1] : i == NΩ ? Ωg[end] - Ωg[end - 1] : (Ωg[i + 1] - Ωg[i - 1]) / 2 for i in 1:NΩ]
Wg = Ωg .^ (-7 / 3)                    # inspiral-weighted band
Kg = klmtv_K.(Ωg)
hS2 = [1 / abs2(transfer_matrices(P0, Ω)[3][2]) * 2Kg[i] for (i, Ω) in enumerate(Ωg)]   # h²_SQL from the model
μ = Wg .* hS2 ./ (2 .* Kg) .* dΩ    # positive measure

"J for an arbitrary readout coordinate profile q(Ω)."
Jofq(q) = sum(μ .* ((q .- Kg) .^ 2 .+ 1))

# A NESTED basis of real-rational readout coordinates with FIXED poles, chosen so
# that the exact answer 𝒦 = Λ⁴/(Ω²(γ²+Ω²)) is NOT in the span of any finite piece.
const POLES = [0.19, 4.7, 0.83, 1.9, 0.42, 3.1, 1.3, 0.28]     # γ = 1 deliberately absent
function basis(nb)
    Φ = zeros(NΩ, nb)
    for k in 1:nb, i in 1:NΩ
        Φ[i, k] = 1 / (Ωg[i]^2 * (POLES[k]^2 + Ωg[i]^2))
    end
    return Φ
end

println("  n (basis)   min eig(H)     J*(n)         J*/J_opt      excess [dB]")
for nb in 1:6
    Φ = basis(nb)
    H = Φ' * (μ .* Φ)                # Gram matrix = ½·Hessian of J
    g = Φ' * (μ .* Kg)
    c = H \ g
    q = Φ * c
    r = Jofq(q) / Jofq(Kg)
    @printf(
        "  %6d     %+.3e   %.8e   %.8f   %+8.4f\n",
        nb, minimum(eigvals(Symmetric(H))), Jofq(q), r, 10log10(r)
    )
end
@printf("\n  exact optimum J(q=𝒦) = %.8e  (attained by 2 lossless filter cavities)\n", Jofq(Kg))
println("  min eig(H) > 0 at every order ⟹ J is a strictly convex quadratic in the")
println("  basis coefficients; the QP optimum is one linear solve, no search.")

# --- the same objective in the hardware coordinates (one filter cavity) ---
println()
println("  Hardware coordinates (ONE detuned cavity + homodyne θ): J(ξ, δ, θ)")
const sub = 1:10:NΩ                                   # coarser grid for the 2-D scan
function Jhw(ξ, δ, θ)
    tot = 0.0
    @inbounds for i in sub
        Ω = Ωg[i]
        ψ = θ - rotation_angle(ξ, δ, Ω)
        sψ = sin(ψ)
        abs(sψ) < 1.0e-12 && return Inf
        q = cos(ψ) / sψ
        tot += μ[i] * ((q - Kg[i])^2 + 1)
    end
    return tot
end

Jopt_sub = sum(μ[i] for i in sub)          # J at q = 𝒦 on the same subgrid

# (i) best achievable with ONE cavity: coarse scan then shrinking-grid refinement
function best_one_cavity()
    ξlo, ξhi = -6.0, 6.0
    δlo, δhi = 0.03, 30.0
    θlo, θhi = 0.02π, 0.98π
    best = (Inf, 0.0, 0.0, 0.0)
    for _ in 1:8
        for θ in range(θlo, θhi; length = 15),
                ξ in range(ξlo, ξhi; length = 25),
                δ in exp.(range(log(δlo), log(δhi); length = 25))
            v = Jhw(ξ, δ, θ)
            v < best[1] && (best = (v, ξ, δ, θ))
        end
        _, ξ0, δ0, θ0 = best
        wξ, wδ, wθ = (ξhi - ξlo) / 6, (log(δhi) - log(δlo)) / 6, (θhi - θlo) / 6
        ξlo, ξhi = ξ0 - wξ, ξ0 + wξ
        δlo, δhi = exp(log(δ0) - wδ), exp(log(δ0) + wδ)
        θlo, θhi = max(1.0e-3, θ0 - wθ), min(π - 1.0e-3, θ0 + wθ)
    end
    return best
end
best = best_one_cavity()
@printf("  best 1-cavity readout : ξ = %+.3f  δ/γ = %.3f  θ/π = %.3f\n", best[2], best[3], best[4] / π)
@printf(
    "                          J = %.6e   (%.4f × J_opt,  %+.3f dB excess)\n",
    best[1], best[1] / Jopt_sub, 10log10(best[1] / Jopt_sub)
)
@printf("  (2 cavities reach J_opt = %.6e exactly — pole placement)\n", Jopt_sub)

# (ii) a CERTIFICATE of non-convexity: two points whose midpoint beats the chord
println("\n  Certificate that J is non-convex in the hardware coordinates (ξ, δ, θ):")
found = false
for (p1, p2) in [
        ((-2.0, 0.5, π / 2), (2.0, 0.5, π / 2)),
        ((-1.5, 0.3, π / 2), (1.5, 3.0, π / 2)),
        ((0.5, 0.2, 0.3π), (0.5, 5.0, 0.7π)),
        ((-3.0, 1.0, 0.4π), (3.0, 1.0, 0.6π)),
    ]
    J1, J2 = Jhw(p1...), Jhw(p2...)
    mid = ((p1[1] + p2[1]) / 2, (p1[2] + p2[2]) / 2, (p1[3] + p2[3]) / 2)
    Jm = Jhw(mid...)
    if isfinite(J1) && isfinite(J2) && Jm > (J1 + J2) / 2
        @printf("    x₁ = (ξ %+.2f, δ %.2f, θ/π %.2f)  J₁ = %.4e\n", p1[1], p1[2], p1[3] / π, J1)
        @printf("    x₂ = (ξ %+.2f, δ %.2f, θ/π %.2f)  J₂ = %.4e\n", p2[1], p2[2], p2[3] / π, J2)
        @printf(
            "    midpoint                            J  = %.4e  >  (J₁+J₂)/2 = %.4e  ✗convex\n",
            Jm, (J1 + J2) / 2
        )
        global found = true
        break
    end
end
found || println("    (no violation among the probe pairs tried)")
println("\n  ⟹ the SAME objective is a positive-definite quadratic in q and")
println("     non-convex in (ξ, δ, θ). Convexity is a property of the coordinates.")

# ------------------------------------------------------------------ #
# E.  With loss, the optimal readout direction leaves the real axis.  #
# ------------------------------------------------------------------ #
println()
println("="^72)
println("E.  When can a LOSSLESS PASSIVE filter chain reach the optimum?")
println("    Such a chain gives w = (scalar phase)·(REAL 2-vector).  It attains the")
println("    optimum iff the optimal direction v ∝ Gh†Σ⁻¹ is real UP TO A PHASE.")
println("    Diagnostic:  ρ ≡ |sin(arg v₁ − arg v₂)| ∈ [0,1];  ρ = 0 ⟺ realizable.")
println("="^72)

println("\n  (a) TUNED arm cavity, loss port co-located with the input mirror:")
@printf("  %8s %8s %12s %12s %14s %14s %9s\n", "γl/γ", "Ω/γ", "ρ", "||χ|-1|", "S_min(any)", "S_min(real)", "gap [dB]")
for γl in [1.0e-3, 1.0e-2, 1.0e-1], Ω in [0.2, 1.0, 5.0]
    g = readout_gap(klmtv_plant(γl = γl), Ω)
    @printf("  %8.0e %8.2f %12.3e %12.3e %14.6e %14.6e %9.5f\n", γl, Ω, g.ρ, g.χdef, g.Sany, g.Sreal, 10log10(g.Sreal / g.Sany))
end
println("  ⟹ ρ = 0: KLMTV's lossless passive 2-cavity chain stays EXACTLY optimal")
println("     under this loss.  Two structural reasons, both visible in (A,B,C,D):")
println("     Σ stays real (the loss column is a real multiple of Ga + I), and the")
println("     signal enters one quadrature only, Gh ∝ (0,1).")

println("\n  (b) DETUNED arm cavity (Δc ≠ 0):")
@printf("  %8s %8s %12s %12s %14s %14s %9s\n", "Δc/γ", "Ω/γ", "ρ", "||χ|-1|", "S_min(any)", "S_min(real)", "gap [dB]")
for Δc in [0.3, 1.0], Ω in [0.2, 0.5, 1.0, 3.0]
    g = readout_gap(klmtv_plant(γl = 1.0e-2, Δc = Δc), Ω)
    @printf("  %8.2f %8.2f %12.3e %12.3e %14.6e %14.6e %9.5f\n", Δc, Ω, g.ρ, g.χdef, g.Sany, g.Sreal, 10log10(g.Sreal / g.Sany))
end
println("  ⟹ ρ ≠ 0.  Note Σ is STILL real here — what goes complex is the direction")
println("     of the SIGNAL vector Gh, because detuning puts strain into both")
println("     quadratures with unequal phase.  No lossless PASSIVE filter chain of")
println("     any order attains the optimum.  Closing the gap requires leaving that")
println("     class: either deliberate asymmetric attenuation (pays vacuum) or an")
println("     ACTIVE, squeezing element in the chain (does not).")
println()
println("done.")
