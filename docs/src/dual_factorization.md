# Factoring the KLMTV plant

Script: `examples/klmtv_jspectral.jl`, on top of the toolbox of
[Chain scattering and J-lossless factorization](@ref). Equation and theorem
numbers are Kimura's.

## Result

KLMTV's two filter cavities come out of a dual (J,J′)-lossless factorization of
the chain-scattering representation of the interferometer, with no ansatz, no
trigonometry and no optimization. Against KLMTV Eq. (89) at `I₀/I_SQL = 1`:

| | `ξ_I` | `δ_I/γ` | `ξ_II` | `δ_II/γ` |
|---|---|---|---|---|
| KLMTV Eq. (89) | +1.7671310 | 0.5156013 | −0.2772297 | 1.3017527 |
| factorization, `μ = k = 10⁻³` | +1.7698072 | 0.5145005 | −0.2776593 | 1.2989507 |
| factorization, `μ = k = 10⁻⁴` | +1.7686190 | 0.5147564 | −0.2775913 | 1.2993195 |

Agreement is 0.1–0.2 %, limited by the finite suspension and the finite
`lev − 1` used, not by the factorization, which is exact. In addition:

1. **Feasibility is an inertia test.** The optimal level `lev* = 1` is returned
   as a signature flip of a constant matrix, before any Riccati equation is
   solved.
2. **`det(HJH~) = 1 − lev²` identically.**
3. **The reality defect is priced**: the suspension that Lemma 5.5 forces is
   itself what makes the optimal readout complex.

## The generalized plant

**The calibration constraint fixes the port count.** A common scalar factor on
the readout row `K = (K₁, K₂)` is free. Spend it by setting `K₂ ≡ 1`, leaving one
design variable, the readout coordinate `q = K₁/K₂`. With
`M = [[1,0],[−𝒦,1]]`,

```
R(Ω) := S_h/S_h^min = ‖(q, 1) M‖² = |q − 𝒦|² + 1.
```

Fixing `K₂ = 1` means the compensator sees only `b₁`:

```
z = weighted error   dim m = 1        u = q·b₁          dim p = 1
w = (a₁, a₂)         dim r = 2        y = b₁            dim q = 1
```

`P₂₁` is 1×2, so `CHAIN(P)` does not exist — but `m = p`, so this is Kimura's
dual two-block case (7.30): `H = DCHAIN(P)`, `Φ = DHM(H; K)`, Theorem 7.5.
Which of Kimura's cases applies is decided by how the free scalar is spent.

**The performance weight collapses to a constant.** Take
`S_target = S_h^min`. Because `M` is unit triangular, `M⁻¹(0,1)ᵀ = (0,1)ᵀ`
whatever `𝒦` is, so `S_h^min = 1/|g|²` exactly, with `g(s) = √(2γ)GL/(s + γ)`
the signal transfer — untouched by the suspension (verified to `9e−11`). The
error weight is `W = 1/g` and `W·G_h = (0,1)ᵀ` is constant: the specification
enters as one number in the `D` matrix.

**The plant.** With a suspension (required, see below),

```
𝒦(s) = 2γG²/((s² − γ²)(m s² + μ s + k))
```

([`ponderomotive_gain`](@ref), realized by [`ponderomotive_gain_ss`](@ref)) —
note the anti-stable pole at `s = +γ`, from dividing out the all-pass
`e^{−2iβ}`. Verified against `−G_n[2,1]/G_n[1,1]` to `1.3e−15`. With the level
normalization `diag(1/lev, 1)` of (7.36), the dual chain-scattering
representation is ([`dual_chain_plant`](@ref))

```
H(lev) = [ lev   𝒦   −1 ]
         [  0    1    0 ]
```

## Feasibility is a signature flip

```
H J_{1,2} H~ = − [ 1 + 𝒦𝒦~ − lev²    𝒦  ]      det = 1 − lev²,
                 [       𝒦~           1  ]
```

verified to `~1e−10` across four decades. The `(2,2)` entry is `−1`, so the
signature is `(1,1)` iff `lev > 1`. Kimura's condition (6.50) evaluated on the
constant `D` matrix alone, `D J_{1,2} Dᵀ = diag(lev² − 1, −1)`, gives the same
answer:

| `lev` | inertia of `HJH~` | (6.50) |
|---|---|---|
| 0.500, 0.900, 0.999 | (0,2) | infeasible |
| 1.001, 1.050, 5.000 | (1,1) | feasible |

So `lev* = 1` — nothing beats `S_h^min` — is an inertia count on a 2×2 constant
matrix.

## The factorization

[`dual_factor_design`](@ref) calls [`dual_jj_factor`](@ref) and evaluates the
readout. At `μ = k = 0.02`, `lev = 1.05` (relative errors):

```
Ω · Ψ = H                             1.8e−13
H J_{1,2} H~ = Ω J_{1,1} Ω~           1.9e−13
Ψ  D J Dᵀ = J′            (4.74)      2.0e−15
Ψ J_{1,2} Ψ~ = J_{1,1}    (4.71)      2.0e−13

Ω  unimodular, degree 4; Ω and Ω⁻¹ poles −1, −1, −0.0100 ± 0.1411i, all stable
Ψ  degree 2, poles −1, +1
```

`Ω` carries the four stable poles of `HJH~` (mirror images of the arm pole and
the two suspension poles); `Ψ` carries the `±γ` pair, absorbing the anti-stable
pole of `𝒦`.

**The readout.** Theorem 7.5 gives `q = DHM(Ω⁻¹; S)` for a free scalar
contractive stable `S`. With `S = 0`:

| `lev` | factorization | `sup R` | `√(sup R)` |
|---|---|---|---|
| 0.9, 0.99, 1.0 | infeasible (inertia) | | |
| 1.00001 | exists | 1.0000000 | 1.0000000 |
| 1.001 | exists | 1.0000681 | 1.0000340 |
| 1.01 | exists | 1.0059160 | 1.0029536 |
| 1.1 | exists | 1.2095176 | 1.0997807 |
| 1.5 | exists | 1.6583843 | 1.2877827 |
| 3.0 | exists | 1.8730491 | 1.3685939 |
| 10.0 | exists | 1.9198896 | 1.3856008 |

`√(sup R) ≤ lev` throughout. `S = 0` is one point of the solution set, not the
equalizing one, so it generally does better than the level asked for, and it
saturates near `sup R ≈ 1.92` for large `lev`.

**`q → 𝒦` as `lev → 1`:**

| `lev − 1` | `max |q − 𝒦|` | `sup R − 1` |
|---|---|---|
| 1e−1 | 4.577e−01 | 2.095e−01 |
| 1e−2 | 7.686e−02 | 5.907e−03 |
| 1e−3 | 8.243e−03 | 6.794e−05 |
| 1e−4 | 8.303e−04 | 6.894e−07 |
| 1e−5 | 8.309e−05 | 6.904e−09 |

`|q − 𝒦| ∝ (lev − 1)^{1/2}` and `sup R − 1 = |q − 𝒦|²`. The variational
condition `u ∝ (𝒦, 1)` reappears as the `lev → 1` limit of a family of
factorizations.

## The cavities

`q` at `lev = 1 + 10⁻⁵` is fitted to `N₀/(D₀ + D₁x + x²)` in `x = Ω²` by
relative-weighted least squares on `Re q`, and the cavities extracted as the
roots of `D − iN` ([`cavities_from_rational`](@ref)). `(*)` marks rows where the
suspension resonance `Ω₀ = √(k/m)` falls inside the band `Ω/γ ∈ [0.05, 20]` and
contaminates the fit.

```
mu      k       Om_0    Im-q cost/dB   cavity I  (xi, delta)      cavity II (xi, delta)
3e-01   3e-01   0.548   19.8           (+3.1934552, 0.0651613)  (-0.4441789, 0.1747194) (*)
1e-01   1e-01   0.316   35.3           (+3.1353010, 0.0582897)  (-0.5088111, 0.1446949) (*)
3e-02   3e-02   0.173   51.5           (+2.4731207, 0.1456704)  (-0.4162940, 0.3550536) (*)
1e-02   1e-02   0.100   65.9           (+2.4202537, 0.1595315)  (-0.4130254, 0.3861788) (*)
3e-03   3e-03   0.055   81.7           (+2.6454924, 0.0947191)  (-0.4553108, 0.2283163) (*)
1e-03   1e-03   0.032   32.9           (+1.7698072, 0.5145005)  (-0.2776593, 1.2989507)
3e-04   3e-04   0.017   15.9           (+1.7688484, 0.5147156)  (-0.2775974, 1.2992866)
1e-04   1e-04   0.010    6.0           (+1.7686190, 0.5147564)  (-0.2775913, 1.2993195)
```

The unmarked rows converge monotonically to KLMTV Eq. (89).

**The reality defect, priced.** The `Im-q cost` column is `10 log₁₀` of the
excess incurred by discarding `Im q`, i.e. by insisting on a real `q`, which is
all a lossless passive chain can realize. It falls `32.9 → 15.9 → 6.0 dB` as
`μ = k` goes `10⁻³ → 3·10⁻⁴ → 10⁻⁴`, consistent with `ρ = 0` exactly for the
tuned free-mass plant. The suspension that Lemma 5.5 requires for the
factorization to exist is itself what makes the optimal readout complex; for
KLMTV the tension is benign, since both errors vanish together as the
suspension is relaxed.

## The free mass

```
mu     k      poles of H on jw   dual_jj_factor(H, 1, 1)
0      0      2                  Riccati (6.62) has no stabilizing solution
0.1    0      1                  Riccati (6.62) has no stabilizing solution
0      0.1    2                  Riccati (6.62) has no stabilizing solution
0.02   0.02   0                  exists
```

Only a damped suspension admits a factorization (Lemma 5.5, on the factored
object itself).

## Known limitation

`Ω` is verified directly against the definition: unimodular, `ΩΨ = H`, and `Ψ`
dual (J,J′)-**unitary**. But conditions (ii)/(iii) of Theorem 6.12 as
transcribed do **not** hold for this plant: `Y` and `Ȳ` come out rank-1
negative semidefinite (`min eig Y = −0.205`, `min eig Ȳ = −4.085`) where the
theorem requires `≥ 0`. Structurally, the only nonzero column of `B` sits in
the `w₁` slot, which carries `−1` in `J_{1,2}`, so `BJBᵀ = −B_kB_kᵀ ≤ 0` and the
Riccati (6.62) is an "anti-LQR" whose stabilizing solution is `≤ 0`. The
transcription of (6.63) for `Ψ` likewise disagrees with `Ω⁻¹H`, so the code
computes `Ψ = Ω⁻¹H` and checks it.

Whether this is a sign/convention error in the transcription or a genuine
failure of (ii)/(iii) for this plant is not settled. Consequently `Ψ` is
verified unitary but not lossless, Theorem 4.15 cannot be invoked, and the
completeness of the `S`-parameterization is not established. Everything above is
verified by direct evaluation of the `S = 0` design and does not depend on it.
This is why [`dual_factor_design`](@ref) uses `strict = false`.
