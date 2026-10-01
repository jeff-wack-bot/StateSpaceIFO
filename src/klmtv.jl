# The KLMTV single-arm plant as a 4-state quadrature model, and the readout
# quantities computed from it.

"""
    klmtv_plant(; γ = 1, G = 1, m = 1, L = 1, μ = 0, k = 0, γl = 0, Δc = 0)

Dark-port KLMTV plant in the two-photon quadrature picture, states
`x = (c₁, c₂, q, p)` (two intracavity quadratures, test-mass displacement and
momentum):

    ċ₁ = -γt c₁ - Δc c₂          + √(2γ) a₁ + √(2γl) l₁
    ċ₂ = -γt c₂ + Δc c₁ + G q    + √(2γ) a₂ + √(2γl) l₂ + G L h
    q̇  = p/m
    ṗ  = G c₁ - k q - μ p
    b  = -a + √(2γ) c,           γt = γ + γl

`γ` is the input-coupler half-bandwidth, `γl` an intracavity loss rate with
its own vacuum input `l`, `Δc` the arm detuning (a symplectic rotation), and
`k`, `μ` an optional suspension spring and damping.  The radiation-pressure
and phase-modulation couplings are the same `G` because both come from
`H_int = -G c₁ q`.

Returns a NamedTuple with the matrices `A, Ba, Bl, Bh, C, Da` and the
systems `Gn` (vacuum `a` → `b`), `Gl` (loss vacuum → `b`) and `Gh`
(strain → `b`).
"""
function klmtv_plant(; γ = 1.0, G = 1.0, m = 1.0, L = 1.0, μ = 0.0, k = 0.0, γl = 0.0, Δc = 0.0)
    γt = γ + γl
    A = [
        -γt  -Δc   0.0   0.0;
        Δc  -γt    G     0.0;
        0.0  0.0   0.0   1 / m;
        G    0.0  -k    -μ
    ]
    Ba = [sqrt(2γ) 0.0; 0.0 sqrt(2γ); 0.0 0.0; 0.0 0.0]
    Bl = [sqrt(2γl) 0.0; 0.0 sqrt(2γl); 0.0 0.0; 0.0 0.0]
    Bh = reshape([0.0, G * L, 0.0, 0.0], 4, 1)
    C = [sqrt(2γ) 0.0 0.0 0.0; 0.0 sqrt(2γ) 0.0 0.0]
    Da = -Matrix{Float64}(I, 2, 2)
    return (
        A = A, Ba = Ba, Bl = Bl, Bh = Bh, C = C, Da = Da,
        Gn = SS(A, Ba, C, Da), Gl = SS(A, Bl, C, zeros(2, 2)), Gh = SS(A, Bh, C, zeros(2, 1)),
    )
end

"""
    transfer_matrices(P, Ω) -> (Gn, Gl, Gh)

Transfer matrices of the plant `P` (from [`klmtv_plant`](@ref)) at sideband
frequency `Ω`, `s = +iΩ`.
"""
transfer_matrices(P, Ω) = (evalfr(P.Gn, im * Ω), evalfr(P.Gl, im * Ω), evalfr(P.Gh, im * Ω))

"""
    klmtv_K(Ω; γ = 1, G = 1, m = 1)

KLMTV's back-action gain for a free test mass,
`𝒦(Ω) = 2γG²/(mΩ²(γ² + Ω²))` (so `Λ⁴ = 2γG²/m`).
"""
klmtv_K(Ω; γ = 1.0, G = 1.0, m = 1.0) = 2γ * G^2 / (m * Ω^2 * (γ^2 + Ω^2))

"""
    ponderomotive_gain(s; γ = 1, G = 1, m = 1, μ = 0, k = 0)

The back-action gain as a rational function of the Laplace variable,
`𝒦(s) = 2γG²/((s² - γ²)(ms² + μs + k))`, including a suspension.  On the
axis, `ponderomotive_gain(iΩ)` equals [`klmtv_K`](@ref)`(Ω)` when
`μ = k = 0`.  It equals `-Gn[2,1]/Gn[1,1]` of the tuned, lossless plant.
"""
ponderomotive_gain(s; γ = 1.0, G = 1.0, m = 1.0, μ = 0.0, k = 0.0) =
    2γ * G^2 / ((s^2 - γ^2) * (m * s^2 + μ * s + k))

"""
    ponderomotive_gain_ss(; γ = 1, G = 1, m = 1, μ = 0, k = 0) -> SS

Controllable-canonical realization of [`ponderomotive_gain`](@ref), with
denominator `m s⁴ + μ s³ + (k - γ²m) s² - γ²μ s - γ²k`.
"""
function ponderomotive_gain_ss(; γ = 1.0, G = 1.0, m = 1.0, μ = 0.0, k = 0.0)
    a = [-γ^2 * k, -γ^2 * μ, k - γ^2 * m, μ] ./ m
    Ak = [0.0 1 0 0; 0 0 1 0; 0 0 0 1; -a[1] -a[2] -a[3] -a[4]]
    Bk = reshape([0.0, 0, 0, 1], 4, 1)
    Ck = [2γ * G^2 / m 0.0 0 0]
    return SS(Ak, Bk, Ck, zeros(1, 1))
end

"""
    backaction_gain(Gn)

`𝒦 = -Gn[2,1]/Gn[1,1]` read off a 2×2 noise transfer matrix — the
definition that applies when `𝒦` has no closed form.
"""
backaction_gain(Gn::AbstractMatrix) = -Gn[2, 1] / Gn[1, 1]

"""
    arm_phase(Gn)

`β = -½ arg Gn[1,1]`.  With `s = +iΩ` the arm phase appears as `e^{-2iβ}`.
"""
arm_phase(Gn::AbstractMatrix) = -angle(Gn[1, 1]) / 2

# --- readout -----------------------------------------------------------------

"""
    sh_min(Gh, Σ)

Minimum strain-referred noise over all readout directions,
`S_h^min = 1/(Gh† Σ⁻¹ Gh)`, with `Σ = Σ_ports G G†` the output noise
covariance at one frequency (unit vacuum).
"""
sh_min(Gh, Σ) = real(1 / (Gh' * (Σ \ Gh))[1])

"""
    optimal_readout(Gh, Σ)

The readout row `v ∝ Gh† Σ⁻¹` that attains [`sh_min`](@ref), as a vector.
"""
optimal_readout(Gh, Σ) = vec(Gh' / Σ)

"""
    reality_defect(v)

`ρ = |sin(arg v₁ - arg v₂)| ∈ [0, 1]`.  A lossless passive filter chain
produces readout directions that are real up to a common phase, so it can
attain the optimum iff `ρ = 0`.
"""
function reality_defect(v)
    (abs(v[1]) < 1.0e-300 || abs(v[2]) < 1.0e-300) && return 0.0
    return abs(imag(v[1] * conj(v[2]))) / (abs(v[1]) * abs(v[2]))
end

"""
    double_angle_phasor(z)

`χ(z) = (z₁ + iz₂)/(z₁ - iz₂)`.  For a real direction at angle `ψ`,
`χ = e^{2iψ}`; `χ(R(α)z) = e^{2iα}χ(z)`; and `|χ(v)| = 1 ⟺ ρ(v) = 0`.
"""
double_angle_phasor(z) = (z[1] + im * z[2]) / (z[1] - im * z[2])

"""
    sh_min_real(Gh, Σ)

Minimum strain-referred noise over **real** readout directions `u ∝ (q, 1)`,
`q ∈ ℝ ∪ {∞}`, i.e. what a lossless passive filter chain can reach.  Solved
exactly from the stationarity condition, a quadratic in `q`.
"""
function sh_min_real(Gh, Σ)
    Σr = real(Σ); Hr = real(Gh * Gh')
    a, b, c = Σr[1, 1], Σr[1, 2], Σr[2, 2]
    p, r, s = Hr[1, 1], Hr[1, 2], Hr[2, 2]
    nd = Tuple{Float64, Float64}[(a, p)]                # q → ∞  (u = (1,0))
    qs = Float64[]
    A2, B2, C2 = (a * r - b * p), (a * s - c * p), (b * s - c * r)
    if abs(A2) > 1.0e-280
        disc = B2^2 - 4A2 * C2
        disc >= 0 && append!(qs, [(-B2 + sqrt(disc)) / (2A2), (-B2 - sqrt(disc)) / (2A2)])
    elseif abs(B2) > 1.0e-280
        push!(qs, -C2 / B2)
    end
    for q in qs
        push!(nd, (a * q^2 + 2b * q + c, p * q^2 + 2r * q + s))
    end
    return minimum(d <= 0 ? Inf : n / d for (n, d) in nd)
end

"""
    readout_gap(P, Ω) -> (ρ, χdef, Sany, Sreal)

At sideband frequency `Ω`, for the plant `P`: the reality defect `ρ` of the
optimal readout, `χdef = ||χ(v)| - 1|`, the best noise over all readout
directions `Sany` ([`sh_min`](@ref)) and over real directions only `Sreal`
([`sh_min_real`](@ref)).  The noise covariance includes the loss port.
"""
function readout_gap(P, Ω)
    Ga, Gl, Gh = transfer_matrices(P, Ω)
    Σ = Ga * Ga' + Gl * Gl'
    v = optimal_readout(Gh, Σ)
    χdef = abs(abs(double_angle_phasor(v)) - 1)
    return (ρ = reality_defect(v), χdef = χdef, Sany = sh_min(Gh, Σ), Sreal = sh_min_real(Gh, Σ))
end

"""
    excess_ratio(p, Ω; K = klmtv_K(Ω))

`R = S_h/S_h^min = (q - 𝒦)² + 1` for a chain of lossless filter cavities
followed by a homodyne at angle `θ`, with `p = [ξ₁, δ₁, ξ₂, δ₂, …, θ]` and
readout coordinate `q = cot(θ - Σⱼ α_rot,j(Ω))`.  Returns `Inf` where the
readout is blind to the signal.
"""
function excess_ratio(p, Ω; K = klmtv_K(Ω))
    θ = p[end]
    a = 0.0
    for j in 1:2:(length(p) - 1)
        a += rotation_angle(p[j], p[j + 1], Ω)
    end
    ψ = θ - a
    sψ = sin(ψ)
    abs(sψ) < 1.0e-14 && return Inf
    return (cos(ψ) / sψ - K)^2 + 1
end

# --- the dual chain-scattering representation, and its factorization --------

"""
    dual_chain_plant(Kss, lev) -> SS

Dual chain-scattering representation of the KLMTV readout problem with the
calibration `K₂ ≡ 1` and level normalization `diag(1/lev, 1)`:

    H(lev) = [ lev   𝒦   -1 ]
             [  0    1    0 ]

`Kss` is a realization of `𝒦(s)`, e.g. [`ponderomotive_gain_ss`](@ref).
Then `det(H J_{1,2} H~) = 1 - lev²`.
"""
function dual_chain_plant(Kss::SS, lev)
    n = nx(Kss)
    return SS(
        Kss.A, [zeros(n, 1)  Kss.B  zeros(n, 1)], [Kss.C; zeros(1, n)],
        [lev 0.0 -1.0; 0.0 1.0 0.0]
    )
end

"""
    dual_factor_design(lev; γ = 1, G = 1, m = 1, μ, k, Ωs, S = 0, strict = false)

Factor `H(lev)` ([`dual_chain_plant`](@ref)) by [`dual_jj_factor`](@ref) and
read off the readout coordinate `q = DHM(Ω⁻¹; S)` and the achieved excess
ratio `R(Ω) = |q - 𝒦|² + 1` on the grid `Ωs`.  A suspension (`μ > 0`,
`k > 0`) is required for the factorization to exist.

Returns `(ok, why, Ω, Ψ, iΩ, Y, Ȳ, minY, minȲ, ρ, Ωs, q, R, eH, eJ, supR)`,
where `eH` and `eJ` are the relative errors of `ΩΨ = H` and
`H J H~ = Ω J′ Ω~` over the grid.
"""
function dual_factor_design(
        lev; γ = 1.0, G = 1.0, m = 1.0, μ, k,
        Ωs = exp10.(range(log10(0.02), log10(50); length = 600)), S = 0.0, strict = false
    )
    Kss = ponderomotive_gain_ss(; γ, G, m, μ, k)
    Kp(s) = ponderomotive_gain(s; γ, G, m, μ, k)
    H = dual_chain_plant(Kss, lev)
    f = dual_jj_factor(H, 1, 1; strict = strict)
    f.ok || return (ok = false, why = f.why)
    eH = eJ = 0.0
    q = ComplexF64[]; R = Float64[]
    for Ω in Ωs
        s = im * Ω
        Om, Ps = evalfr(f.Ω, s), evalfr(f.Ψ, s)
        Hs = evalfr(H, s)
        sc = max(norm(Hs), 1.0)
        eH = max(eH, norm(Om * Ps - Hs) / sc)
        eJ = max(
            eJ, norm(
                Hs * Jsig(1, 2) * transpose(evalfr(H, -s))
                    - Om * Jsig(1, 1) * transpose(evalfr(f.Ω, -s))
            ) / sc^2
        )
        qs = dhm(evalfr(f.iΩ, s), fill(S, 1, 1), 1)[1]
        push!(q, qs); push!(R, abs2(qs - Kp(s)) + 1)
    end
    return (
        ok = true, why = "", Ω = f.Ω, Ψ = f.Ψ, iΩ = f.iΩ, Y = f.Y, Ȳ = f.Ȳ,
        minY = f.minY, minȲ = f.minȲ, ρ = f.ρ,
        Ωs = Ωs, q = q, R = R, eH = eH, eJ = eJ, supR = maximum(R),
    )
end
