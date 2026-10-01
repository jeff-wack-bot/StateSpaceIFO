# jlossless.jl — validation of the chain-scattering toolbox, and which optical
# elements are J-lossless.
#
# Run:  julia --project=examples examples/jlossless.jl
#
# Sections:
#   C. validation against Kimura's published worked examples
#   D. optics: which optical elements are J-lossless, and with which J
#   E. KLMTV: the filter-cavity chain as the J-lossless factor

using StateSpaceIFO
using LinearAlgebra, Printf

hdr(s) = (println(); println("="^76); println(s); println("="^76))
ok(b) = b ? "OK  " : "FAIL"

hdr("C.  Reproduction of Kimura's worked examples")

# --- Example 4.4 / 4.6: J-lossless matrices, Theorem 4.5 test ---
Θ44 = SS([-2.0 0; 0 1], [-4.0 0; 0 2], [1.0 0; 0 1], [1.0 0; 0 1])
res = jjlossless_residual(Θ44, 1, 1)
@printf("Ex 4.4  Theta = diag((s-2)/(s+2), (s+1)/(s-1))  J_{1,1}-lossless\n")
@printf(
    "        D'JD-J' = %.1e   D'JC+B'P = %.1e   sup||T*JT-J'|| = %.1e   min eig P = %+.4f\n",
    res.dD, res.dLyap, res.dFreq, res.minP
)
@printf(
    "        book says P = diag(1/4, 1/2) > 0 : %s\n",
    ok(
        isapprox(
            lyap(transpose(Θ44.A), transpose(Θ44.C) * Jsig(1, 1) * Θ44.C),
            [0.25 0; 0 0.5]; atol = 1.0e-10
        )
    )
)

Θ46 = SS([0.0 -1; -2 -1], [-1.0 -1; -5 / 2 3 / 2], [1.0 5 / 4; -1 3 / 4], [1.0 0; 0 1])
res = jjlossless_residual(Θ46, 1, 1)
@printf("Ex 4.6  degree-2 J-lossless example\n")
@printf(
    "        D'JD-J' = %.1e   D'JC+B'P = %.1e   sup||T*JT-J'|| = %.1e   min eig P = %+.4f  (book: P = diag(1,1/2))\n",
    res.dD, res.dLyap, res.dFreq, res.minP
)

# --- Example 6.3 / 6.4: (J,J')-lossless factorization of an unstable G ---
G64 = SS([-2.0 0; 0 1], [-2.0 2; 0 2], [2.0 0; 0 2], [1.0 -1; 0 1])
f = jj_factor(G64, 1, 1)
@printf("\nEx 6.3/6.4  G = [ (s-2)/(s+2)  -(s-2)/(s+2) ; 0  (s+3)/(s-1) ]\n")
@printf("        factorization exists: %s   %s\n", ok(f.ok), f.why)
if f.ok
    @printf("        X    = %s   (book: [1 0; 0 0])\n", replace(string(round.(f.X; digits = 6)), "\n" => ""))
    @printf("        Xbar = %s   (book: [0 0; 0 1/2])\n", replace(string(round.(f.X̄; digits = 6)), "\n" => ""))
    # book (6.7):  Theta = diag((s-2)/(s+2), (s+1)/(s-1)),  Pi = [1 -1; 0 (s+3)/(s+1)]
    Πref(s) = [1 -1; 0 (s + 3) / (s + 1)]
    local eΘ = 0.0
    local eΠ = 0.0
    for Ω in [0.1, 0.7, 1.3, 5.0, 20.0]
        s = im * Ω
        eΘ = max(eΘ, norm(evalfr(f.Θ, s) * evalfr(f.Π, s) - evalfr(G64, s)))
        # Theta, Pi are unique only up to a constant J'-unitary N (Kimura p.130),
        # so compare the N-invariant product Pi' J' Pi against the book's (6.7).
        eΠ = max(
            eΠ, norm(
                evalfr(f.Π, s)' * Jsig(1, 1) * evalfr(f.Π, s)
                    - Πref(s)' * Jsig(1, 1) * Πref(s)
            )
        )
    end
    r = jjlossless_residual(f.Θ, 1, 1)
    @printf("        Theta*Pi = G          : max err %.2e  %s\n", eΘ, ok(eΘ < 1.0e-9))
    @printf(
        "        Theta is (J,J')-lossless: D-cond %.1e  Lyap %.1e  freq %.1e  min eig P %+.3f\n",
        r.dD, r.dLyap, r.dFreq, r.minP
    )
    @printf("        Pi'J'Pi matches book (6.7): max err %.2e  %s\n", eΠ, ok(eΠ < 1.0e-8))
    @printf(
        "        Pi unimodular (poles %s, zeros of Pi = poles of Pi^-1)\n",
        replace(string(round.(real.(poles(f.Π)); digits = 4)), "\n" => "")
    )
end

# --- Example 6.5: a tall (J_{21}, J_{11})-lossless factorization ---
s2 = sqrt(2)
G65 = SS(
    [3.0 0; 0 -3], [-s2 4 - s2; s2 / 2 s2 / 2], [0.0 -4; 1 -2; 2 -1],
    [s2 / 2 s2 / 2; s2 / 2 s2 / 2; 0 1]
)
f5 = jj_factor(G65, 2, 1)
@printf("\nEx 6.5  tall system, (J_21, J_11)-lossless factorization\n")
@printf("        exists: %s  %s\n", ok(f5.ok), f5.why)
if f5.ok
    @printf("        X    = %s   (book: [1/14 0; 0 2])\n", replace(string(round.(f5.X; digits = 6)), "\n" => ""))
    @printf("        Xbar = %s   (book: [2 0; 0 0])\n", replace(string(round.(f5.X̄; digits = 6)), "\n" => ""))
    local e = 0.0
    for Ω in [0.05, 0.4, 1.0, 3.0, 11.0]
        e = max(e, norm(evalfr(f5.Θ, im * Ω) * evalfr(f5.Π, im * Ω) - evalfr(G65, im * Ω)))
    end
    r = jjlossless_residual(f5.Θ, 2, 1)
    @printf(
        "        Theta*Pi = G: %.2e  %s ;  Theta (J,J')-lossless: %.1e / %.1e / %.1e, min eig P %+.3f\n",
        e, ok(e < 1.0e-8), r.dD, r.dLyap, r.dFreq, r.minP
    )
end

# ======================================================================
# D.  Which optical elements are J-lossless, and with which J
# ======================================================================

hdr("D.  Optical elements as J-unitary chain matrices")

println("D1.  A lossless *passive* element is inner: G~G = I  (J = I, the r = 0 case).")
for (δ, Δ) in [(1.0, 0.0), (0.5156013216, -0.5156013216 * 1.7671310136), (1.3, 2.0)]
    G = filter_cavity(δ, Δ)
    e = maximum(norm(evalfr(G, im * Ω)' * evalfr(G, im * Ω) - I) for Ω in exp10.(range(-2, 2; length = 400)))
    @printf("     delta=%.4f  Delta=%+.4f   sup||G*G - I|| = %.2e   %s\n", δ, Δ, e, ok(e < 1.0e-10))
end

println("\nD2.  Quadrature rotation and common phase.  For a detuned cavity the 2x2")
println("     quadrature matrix is exp(-i*alpha_com) * R(alpha_rot) with")
println("     alpha_pm = 2 atan(xi +- Omega/delta),  xi = -Delta/delta;")
println("     alpha_rot = (alpha_+ + alpha_-)/2 (even in Omega),")
println("     alpha_com = (alpha_+ - alpha_-)/2 (odd).")
let δ = 0.7, Δ = -0.9, ξ = -Δ / δ
    e = 0.0
    for Ω in [0.05, 0.3, 1.0, 2.5, 8.0]
        αrot, αcom = rotation_angle(ξ, δ, Ω), common_phase(ξ, δ, Ω)
        R = [cos(αrot) -sin(αrot); sin(αrot) cos(αrot)]
        e = max(e, norm(evalfr(filter_cavity(δ, Δ), im * Ω) - exp(-im * αcom) * R))
    end
    @printf("     max || G(iOmega) - exp(-i a_com) R(a_rot) || = %.2e   %s\n", e, ok(e < 1.0e-10))
    println("     (the minus sign on a_com is the s = +i*Omega convention, the same one that")
    println("      turns KLMTV's exp(+2i beta) into exp(-2i beta))")
end

println("\nD3.  A degenerate OPA (squeezer) is NOT J-unitary in the quadrature basis")
println("     (it is symplectic), but IS J-unitary with J = diag(1,-1) in the")
println("     (a, a^dagger) basis — the Bogoliubov metric.  This is the indefinite")
println("     signature that parametric gain contributes, distinct from the")
println("     *directional* signature that carries the H-infinity level gamma.")
let rsq = 0.8
    Sq = [exp(rsq) 0.0; 0.0 exp(-rsq)]                    # quadrature picture
    W = [1 im; 1 -im] / sqrt(2)                          # (x1,x2) -> (a, a^dagger)
    Bg = W * Sq * inv(W)
    Jb = Diagonal([1.0, -1.0])
    @printf("     quadrature:  S'S - I           = %.3f   (not unitary)\n", norm(Sq'Sq - I))
    @printf(
        "     quadrature:  S' Sigma S - Sigma = %.1e   (symplectic, Sigma = [0 1; -1 0])\n",
        norm(transpose(Sq) * [0.0 1; -1 0] * Sq - [0.0 1; -1 0])
    )
    @printf(
        "     (a,adag)  :  S' J S - J        = %.1e   J = diag(1,-1)   %s\n",
        norm(Bg' * Jb * Bg - Jb), ok(norm(Bg' * Jb * Bg - Jb) < 1.0e-12)
    )
end

# ======================================================================
# E.  KLMTV: the filter chain as the J-lossless (here: inner) factor
# ======================================================================

hdr("E.  KLMTV filter cavities by pole splitting (Lemma 4.9 / Potapov)")

cav = klmtv_cavities(1.0)
@printf("KLMTV Eqs. (89) at I0/ISQL = 1 (published: xi_I=1.7671310136 d_I=0.5156013216,\n")
@printf("                               xi_II=-0.2772297408 d_II=1.3017526994)\n")
for (k, c) in enumerate(cav)
    @printf("     cavity %-2d  xi = %+.10f   delta/gamma = %.10f\n", k, c.ξ, c.δ)
end

# Build the cascade of the two cavities and check it is inner of degree 4 in Omega
Fchain = filter_cavity(cav[1].δ, cav[1].Δ) * filter_cavity(cav[2].δ, cav[2].Δ)
r = jjlossless_residual(Fchain, 2, 2; Ωs = exp10.(range(-2, 2; length = 300)))
@printf(
    "\ntwo-cavity chain: %d states;  inner (J = J' = I): sup||F*F - I|| = %.2e  %s\n",
    nx(Fchain), r.dFreq, ok(r.dFreq < 1.0e-9)
)
@printf("   poles: %s\n", join([@sprintf("%+.4f%+.4fi", real(λ), imag(λ)) for λ in poles(Fchain)], "  "))

# Potapov split: recover the individual cavities from the cascade.
println("\nLemma 4.9 pole splitting — recover cavity I from the cascade by selecting its poles:")
p1 = poles(filter_cavity(cav[1].δ, cav[1].Δ))
sel = λ -> minimum(abs.(λ .- p1)) < 1.0e-6
Θ1, Θ2 = potapov_split(Fchain, 2, 2, sel)
e1 = maximum(
    norm(evalfr(Θ1, im * Ω) * evalfr(Θ2, im * Ω) - evalfr(Fchain, im * Ω))
        for Ω in exp10.(range(-2, 2; length = 200))
)
@printf("   Theta1*Theta2 = F : %.2e  %s\n", e1, ok(e1 < 1.0e-8))
@printf(
    "   Theta1 poles: %s   (cavity I: %s)\n",
    join([@sprintf("%+.4f%+.4fi", real(λ), imag(λ)) for λ in poles(Θ1)], " "),
    join([@sprintf("%+.4f%+.4fi", real(λ), imag(λ)) for λ in p1], " ")
)
e2 = maximum(
    norm(evalfr(Θ1, im * Ω)' * evalfr(Θ1, im * Ω) - I)
        for Ω in exp10.(range(-2, 2; length = 200))
)
@printf("   Theta1 inner: sup||T*T - I|| = %.2e  %s\n", e2, ok(e2 < 1.0e-9))

# --- the marginal-pole obstruction ------------------------------------
hdr("E2.  The marginal pole: where the Riccati route stops (Lemma 5.5)")

println("Free test mass => A has a double eigenvalue at 0.  Theorem 5.2's Riccati (5.3)")
println("is applied to the plant pair (A,B) with J = diag(1,-1,-1) on (a1,a2,h).\n")
@printf("   %-8s %-8s  %-14s  %s\n", "mu", "k", "poles on jw", "J-lossless conjugator")
for (μ, k) in [(0.0, 0.0), (1.0e-1, 0.0), (0.0, 1.0e-2), (1.0e-2, 1.0e-2), (0.2, 1.0)]
    P = klmtv_plant(μ = μ, k = k)
    λ = eigvals(P.A)
    onaxis = count(abs(real(l)) < 1.0e-8 for l in λ)
    Θc, _ = jlossless_conjugator(P.A, [P.Ba P.Bh], Jsig(1, 2))
    @printf(
        "   %-8g %-8g  %-14d  %s\n", μ, k, onaxis,
        Θc === nothing ? "no   (Lemma 5.5)" : "yes"
    )
end
println(
    """
    Reading.  With a free test mass the plant has a double pole ON the jw axis,
    and Lemma 5.5 says the conjugation Riccati (5.3) then has no stabilizing
    solution: the whole Chapter 5/6 route is unavailable at the plant level, not
    merely ill-conditioned.  The table also shows that neither repair alone
    works — damping leaves the rigid-translation pole at s = 0, a spring alone
    leaves an undamped oscillator pair on the axis.  Both are needed, i.e. a
    *damped* suspension (or, equivalently, a band weight that rolls off at DC and
    is folded into the plant by weight augmentation before factorizing)."""
)

println("\ndone.")
