# klmtv_jspectral.jl — factorizing the regularized KLMTV plant.
#
# Run:  julia --project=examples examples/klmtv_jspectral.jl
#
# ---------------------------------------------------------------------------
# THE GENERALIZED PLANT.
#
# The readout is a row K = (K1, K2) acting on the dark-port output b.  A common
# scalar factor on K is free (offline calibration), so fix it by K2 = 1 and let
# the single design variable be the readout coordinate  q = K1/K2.
# With M = [[1,0],[-Kp,1]],
#
#     S_h / S_h^min  =  || (q,1) M ||^2  =  |q - Kp|^2 + 1  =:  R(Omega),
#
# because S_h^min = 1/|g|^2 exactly (section 1 below), so the target weight
# cancels the signal transfer and leaves no dynamic weight at all.
#
# Fixing K2 = 1 means the compensator sees only b1.  That is the whole content
# of the port count:
#
#     z = weighted error   dim m = 1        u = q*b1          dim p = 1
#     w = (a1, a2)         dim r = 2        y = b1            dim q = 1
#
# so m = p and P21 is NOT square: this is Kimura's *dual* two-block case (7.30),
# H = DCHAIN(P), Phi = DHM(H; K), and Theorem 7.5 applies.  With the level
# normalization P_gam = diag(1/lev, 1) P of (7.36),
#
#     H(lev) = [ lev   Kp   -1 ]        (m+q) x (m+r) = 2 x 3
#              [  0     1    0 ]
#
# and  H J_{1,2} H~ = Omega J_{1,1} Omega~  with  det(H J H~) = 1 - lev^2.
#
# Theorem 7.5: solvable iff H = Omega * Psi (dual (J_{1,2}, J_{1,1})-lossless
# factorization, Theorem 6.12), and then  q = DHM(Omega^{-1}; S),  S in BH-inf.
#
# Predictions tested:
#   (a) det(H J H~) = 1 - lev^2, so the signature flips exactly at lev = 1, and
#       infeasibility is reported by the inertia condition (6.50) — a
#       certificate, reached before any Riccati equation is solved;
#   (b) the S = 0 design achieves sup R <= lev^2;
#   (c) q -> Kp as lev -> 1: the variational readout, recovered as a limit of
#       factorizations rather than by minimizing a Rayleigh quotient;
#   (d) the lev -> 1 filter is KLMTV's two cavities, Eq. (89);
#   (e) with a free test mass the factorization does not exist (Lemma 5.5), and
#       the suspension that repairs it is itself what makes q complex — the
#       reality defect.
# ---------------------------------------------------------------------------

using StateSpaceIFO
using LinearAlgebra, Printf

hdr(s) = (println(); println("="^78); println(s); println("="^78))
ok(b) = b ? "OK  " : "FAIL"
band(lo, hi, n) = exp10.(range(log10(lo), log10(hi); length = n))

# ======================================================================
# 1.  The regularized plant, and the collapse of the target weight
# ======================================================================

hdr("1.  The plant, and why the performance weight collapses to a constant")

let μ = 0.02, k = 0.02
    P = klmtv_plant(μ = μ, k = k)
    Kp(s) = ponderomotive_gain(s; μ, k)
    Kss = ponderomotive_gain_ss(; μ, k)
    ekp = ess = ehs = 0.0
    for Ω in band(0.02, 50, 400)
        T = evalfr(P.Gn, im * Ω)
        ekp = max(ekp, abs(backaction_gain(T) - Kp(im * Ω)) / abs(Kp(im * Ω)))
        Σ = T * transpose(evalfr(P.Gn, -im * Ω))
        gh = evalfr(P.Gh, im * Ω)
        ess = max(ess, abs(real(gh' * (Σ \ gh))[1] - 2 / (1 + Ω^2)) * (1 + Ω^2) / 2)
        ehs = max(ehs, abs(evalfr(Kss, im * Ω)[1] - Kp(im * Ω)) / abs(Kp(im * Ω)))
    end
    @printf("  Kp from -Gn[2,1]/Gn[1,1] vs closed form (rel)   : %.2e  %s\n", ekp, ok(ekp < 1.0e-9))
    @printf("  canonical realization of Kp        (rel)        : %.2e  %s\n", ehs, ok(ehs < 1.0e-9))
    @printf("  1/S_h^min = Gh' Sigma^-1 Gh  vs  |g|^2  (rel)   : %.2e  %s\n", ess, ok(ess < 1.0e-9))
end
println(
    """
    The third line is load-bearing.  M is unit triangular, so M^{-1}(0,1)' = (0,1)'
    whatever Kp is: S_h^min = 1/|g|^2 EXACTLY, and the suspension does not move it.
    Choosing S_target = S_h^min therefore makes the error weight W = 1/g and
    W*Gh = (0,1)' a constant.  The performance specification enters the plant as
    the constant `lev` in the D-matrix, not as a dynamic weight — which is why
    everything below is a signature computation rather than a loop shaping."""
)

# ======================================================================
# 2.  det(H J H~) = 1 - lev^2
# ======================================================================

hdr("2.  The signature flips exactly at lev = 1")

let Kss = ponderomotive_gain_ss(μ = 0.02, k = 0.02), J = Jsig(1, 2)
    @printf("  %-8s %-24s %-20s %s\n", "lev", "max|det(HJH~)-(1-lev^2)|", "inertia of H J H~", "(6.50) inertia of DJD'")
    for lev in [0.5, 0.9, 0.999, 1.001, 1.05, 1.5, 5.0]
        H = dual_chain_plant(Kss, lev); ed = 0.0; sig = Set{Tuple{Int, Int}}()
        for Ω in band(0.02, 50, 300)
            T = evalfr(H, im * Ω) * J * transpose(evalfr(H, -im * Ω))
            ed = max(ed, abs(det(T) - (1 - lev^2)))
            λ = eigvals(Hermitian((T + T') / 2))
            push!(sig, (count(>(0), λ), count(<(0), λ)))
        end
        E = signature_sqrt(H.D * J * transpose(H.D), 1, 1)
        @printf(
            "  %-8.3f %-24.2e %-20s %s\n", lev, ed,
            join(["($a,$b)" for (a, b) in sig], " "),
            E === nothing ? "wrong -> INFEASIBLE" : "(1,1) -> feasible"
        )
    end
end

# ======================================================================
# 3.  The factorization H = Omega * Psi
# ======================================================================

hdr("3.  Dual (J_{1,2}, J_{1,1})-lossless factorization, Theorem 6.12")

let lev = 1.05
    d = dual_factor_design(lev; μ = 0.02, k = 0.02)
    @printf("  lev = %.3f      (all errors relative)\n", lev)
    @printf("  Omega * Psi = H                       : %.2e  %s\n", d.eH, ok(d.eH < 1.0e-8))
    @printf("  H J_{1,2} H~ = Omega J_{1,1} Omega~   : %.2e  %s\n", d.eJ, ok(d.eJ < 1.0e-8))
    r = dual_lossless_residual(d.Ψ, 1, 1; Ωs = band(0.02, 50, 200))
    @printf("  Psi D J D' = J'  (4.74)               : %.1e  %s\n", r.dD, ok(r.dD < 1.0e-9))
    @printf("  Psi J_{1,2} Psi~ = J_{1,1}  (4.71)    : %.1e  %s\n", r.dFreq, ok(r.dFreq < 1.0e-7))
    Ωm, iΩm = minreal(d.Ω), minreal(d.iΩ)
    @printf("  Omega unimodular: deg %d\n", nx(Ωm))
    @printf(
        "      Omega   poles %s   all stable: %s\n",
        join([@sprintf("%+.4f%+.4fi", real(λ), imag(λ)) for λ in poles(Ωm)], " "),
        ok(isstable(Ωm.A))
    )
    @printf(
        "      Omega^-1 poles %s   all stable: %s\n",
        join([@sprintf("%+.4f%+.4fi", real(λ), imag(λ)) for λ in poles(iΩm)], " "),
        ok(isstable(iΩm.A))
    )
    @printf(
        "  Psi: deg %d, poles %s\n", nx(d.Ψ),
        join([@sprintf("%+.4f%+.4fi", real(λ), imag(λ)) for λ in poles(d.Ψ)], " ")
    )
    @printf(
        "  Theorem 6.12 as literally stated: min eig Y = %+.4e, min eig Ybar = %+.4e, sigma(Y Ybar) = %.4f\n",
        d.minY, d.minȲ, d.ρ
    )
    println(
        """
        CAVEAT, stated plainly.  Omega is verified above against the definition:
        unimodular, and Omega*Psi = H with Psi dual (J,J')-UNITARY.  But conditions
        (ii)/(iii) of Theorem 6.12 as transcribed (Y >= 0, Ybar >= 0) do NOT
        hold: both come out rank-1 NEGATIVE semidefinite.  The structural reason is
        visible in the plant: the only nonzero column of B sits in the w1 slot, which
        carries -1 in J_{1,2}, so B J B' = -Bk Bk' <= 0 and the Riccati (6.62) is an
        'anti-LQR' whose stabilizing solution is <= 0.  Two readings are open --- a
        sign/convention error in the transcription, or a genuine failure of (ii)/(iii)
        for this plant --- and it is not settled which.  Consequence: Psi is verified
        UNITARY but not verified LOSSLESS, so Theorem 4.15 cannot be invoked and the
        COMPLETENESS of the S-parameterization is not established here.  Everything
        claimed below is instead verified by direct evaluation of the S = 0 design."""
    )
end

# ======================================================================
# 4.  Feasibility
# ======================================================================

hdr("4.  Feasibility: lev* = 1 as an inertia certificate")

@printf("  %-9s  %-46s  %-11s %s\n", "lev", "factorization", "sup R", "sqrt(sup R)")
for lev in [0.9, 0.99, 1.0, 1.00001, 1.0001, 1.001, 1.01, 1.1, 1.5, 3.0, 10.0]
    d = dual_factor_design(lev; μ = 0.02, k = 0.02)
    if d.ok
        @printf("  %-9.5f  %-46s  %-11.7f %.7f\n", lev, "exists", d.supR, sqrt(d.supR))
    else
        @printf("  %-9.5f  %-46s\n", lev, d.why)
    end
end
println(
    """
    Below lev = 1 the matrix D J D' has inertia (0,2) instead of (1,1), so (6.50)
    has no solution E and Theorem 6.12 returns INFEASIBLE before touching a
    Riccati equation.  That is a certificate.  Compare klmtv_hinf.jl section 2,
    which located the same boundary by multistart search and needed a self-check
    to stop reporting infeasible optima as optimal.
    sqrt(sup R) <= lev everywhere, with equality approached only near lev ~ 1.1:
    S = 0 is one point of the solution set, not the equalizing one, so it
    generally does BETTER than the level it was asked for.  For large lev the
    S = 0 design saturates near sup R ~ 1.92 rather than degrading to lev^2 --
    loosening the specification stops buying anything once S = 0 is no longer
    binding.  The level is a guarantee, not a prediction of what S = 0 achieves."""
)

# ======================================================================
# 5.  q -> Kp
# ======================================================================

hdr("5.  As lev -> 1, the readout coordinate converges to Kp")

let Ωs = band(0.05, 20, 400), Kp(s) = ponderomotive_gain(s; μ = 0.02, k = 0.02)
    @printf("  %-10s  %-16s  %-16s %s\n", "lev - 1", "max |q - Kp|", "max |q-Kp|/|Kp|", "sup R - 1")
    for δ in [1.0, 0.3, 0.1, 3.0e-2, 1.0e-2, 3.0e-3, 1.0e-3, 1.0e-4, 1.0e-5]
        d = dual_factor_design(1 + δ; μ = 0.02, k = 0.02, Ωs = Ωs)
        d.ok || (println("  infeasible at lev-1 = $δ"); continue)
        ea = maximum(abs(d.q[i] - Kp(im * Ωs[i])) for i in eachindex(Ωs))
        er = maximum(abs(d.q[i] - Kp(im * Ωs[i])) / abs(Kp(im * Ωs[i])) for i in eachindex(Ωs))
        @printf("  %-10.1e  %-16.3e  %-16.3e %.3e\n", δ, ea, er, d.supR - 1)
    end
end
println(
    """
    q -> Kp is the variational condition  u ∝ (Kp, 1)  of KLMTV Eqs. (58)-(59).
    At lev = 1 exactly, H J H~ is singular and the solution set collapses to that
    one direction."""
)

# ======================================================================
# 6.  The filter cavities
# ======================================================================

hdr("6.  Cavity parameters from the lev -> 1 design")

println("  Reference — analytic free-mass Kp = 2/(x^2 + x), x = Omega^2, so N = 2,")
println("  D = x^2 + x, and the two roots of D - iN are the KLMTV cavities:")
for (j, c) in enumerate(cavities_from_rational([2.0], [0.0, 1.0, 1.0]))
    @printf("     cavity %d   xi = %+.10f   delta/gam = %.10f\n", j, c.ξ, c.δ)
end
println("     KLMTV Eq. (89):   xi_I  = +1.7671310136   d_I  = 0.5156013216")
println("                       xi_II = -0.2772297408   d_II = 1.3017526994")

println(
    """
    Now from the numerics.  q(iOm) at lev = 1 + 1e-5 is fitted to N0/(D0 + D1 x + x^2)
    in x = Omega^2 by linear least squares on Re q, and the cavities re-extracted.
    The suspension (mu, k) is relaxed toward the free mass along the way."""
)
let Ωs = band(0.05, 20, 800), δlev = 1.0e-5
    lo, hi = Ωs[1], Ωs[end]
    @printf("  band Omega/gam in [%.2f, %.1f];  a row is marked (*) when the suspension\n", lo, hi)
    @printf("  resonance Om_0 = sqrt(k/m) falls INSIDE that band, which contaminates it.\n\n")
    @printf(
        "  %-9s %-9s %-9s %-13s  %-26s %-26s\n",
        "mu", "k", "Om_0", "Im-q cost/dB", "cavity I  (xi, delta)", "cavity II (xi, delta)"
    )
    for (mu, k) in [
            (0.3, 0.3), (0.1, 0.1), (0.03, 0.03), (0.01, 0.01), (3.0e-3, 3.0e-3),
            (1.0e-3, 1.0e-3), (3.0e-4, 3.0e-4), (1.0e-4, 1.0e-4),
        ]
        Kp(s) = ponderomotive_gain(s; μ = mu, k = k)
        d = dual_factor_design(1 + δlev; μ = mu, k = k, Ωs = Ωs)
        d.ok || (println("  infeasible at mu=$mu"); continue)
        # cost of being restricted to a REAL q, i.e. to lossless passive optics:
        # the achievable excess if the imaginary part must be discarded.
        Rre = [abs2(real(d.q[i]) - Kp(im * Ωs[i])) + 1 for i in eachindex(Ωs)]
        imq = 10log10(maximum(Rre) / d.supR)
        # relative-weighted least squares, so the high-Omega tail (where |q| -> 0)
        # does not dominate the fit
        x = Ωs .^ 2; rq = real.(d.q)
        M = hcat(ones(length(x)), -rq, -rq .* x)        # (N0, D0, D1)
        wt = 1.0 ./ (abs.(rq) .* x .^ 2 .+ 1.0e-12)
        sol = (wt .* M) \ (wt .* (rq .* x .^ 2))
        cav = cavities_from_rational([sol[1]], [sol[2], sol[3], 1.0])
        if length(cav) == 2
            Ω0 = sqrt(k)
            @printf(
                "  %-9.0e %-9.0e %-9.3f %-13.1f  (%+.7f, %.7f)  (%+.7f, %.7f) %s\n",
                mu, k, Ω0, imq, cav[1].ξ, cav[1].δ, cav[2].ξ, cav[2].δ,
                (lo <= Ω0 <= hi) ? "(*)" : ""
            )
        else
            @printf("  %-9.0e %-9.0e %-13.2e  degenerate fit (%d roots)\n", mu, k, imq, length(cav))
        end
    end
end
println(
    """
    The `Im-q cost` column is the reality defect rho, priced.  A lossless PASSIVE
    chain can realize only a REAL q, so the imaginary part of the optimal q is
    exactly the piece passive optics cannot build; the column is 10log10 of the
    excess incurred by discarding it.  The suspension that Lemma 5.5 forces us to
    add in order to have a factorization at all is itself what makes q complex, so
    the two requirements pull against each other and this is the exchange rate
    between them.  It goes to zero with the suspension, which is the quantitative
    version of the finding that KLMTV's tuned, free-mass plant has rho = 0 exactly.
    Read only the unmarked rows: where Om_0 sits inside the band, q has a genuine
    suspension resonance in it and both columns measure that instead."""
)

# ======================================================================
# 7.  The free mass
# ======================================================================

hdr("7.  Free test mass: no factorization (Lemma 5.5)")

@printf("  %-9s %-9s  %-18s  %s\n", "mu", "k", "poles of H on jw", "dual_jj_factor(H, 1, 1)")
for (mu, k) in [(0.0, 0.0), (0.1, 0.0), (0.0, 0.1), (0.02, 0.02)]
    H = dual_chain_plant(ponderomotive_gain_ss(μ = mu, k = k), 1.05)
    non = count(abs(real(λ)) < 1.0e-8 for λ in poles(H))
    f = dual_jj_factor(H, 1, 1; strict = false)
    @printf("  %-9g %-9g  %-18d  %s\n", mu, k, non, f.ok ? "exists" : f.why)
end
println(
    """
    Damping alone leaves the rigid-translation pole at s = 0; a spring alone
    leaves an undamped oscillator pair on the axis; only a damped suspension
    admits a factorization.  This is Kimura's standing hypothesis for Theorems
    7.4/7.5 — no poles or zeros on the jw-axis — biting on a physical plant."""
)

println("\ndone.")
