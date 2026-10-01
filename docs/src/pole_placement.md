# Filter cavities as Blaschke factors

Scripts: `examples/klmtv_state-space.jl` (section C), `examples/klmtv_fig10.jl`.
Source: KLMTV Secs. IV B and V (variational output); cross-checked against
Purdue & Chen 2002, App. A.

Starting from KLMTV's input–output relation, the optimal frequency-dependent
readout filter is found by **pole placement**: the four cavity parameters of
KLMTV Eq. (89) are the roots of one quadratic.

Throughout, write

```
s ≡ Ω²
```

because every quantity here is a rational function of `Ω²`, and one lossless
detuned filter cavity is exactly one first-order pole in `s`.

## Input–output relation and readout

KLMTV Eq. (16) in two-photon matrix notation is

```
b = e^{2iβ} M a + e^{iβ} d h,    M = [[1, 0], [−𝒦, 1]],    d = (√(2𝒦)/h_SQL)·(0, 1)ᵀ.
```

A homodyne readout in direction `u = (u₁, u₂)ᵀ` measures `y = uᵀb`. Referring
the noise to `h` and dropping scalar phases,

```
S_h = (h²_SQL/2𝒦) · [(u₁ − 𝒦u₂)² + u₂²]/u₂²  ≥  h²_SQL/2𝒦,
```

with equality iff `u ∝ (𝒦, 1)ᵀ`: the readout must be orthogonal to the
back-action column `(1, −𝒦)ᵀ` of the plant. This is KLMTV Eqs. (58)–(59)
(`ζ = Φ = arccot 𝒦`) without the inverse trigonometric function.

## The design condition as a transfer-function equation

We want a filter rotation `α_rot(Ω)` and homodyne angle `θ` with
`R(−α_rot) u_θ ∥ (𝒦, 1)ᵀ` at every `Ω`. The condition holds only mod π, so
encode a direction by its double-angle phasor
`χ(z) = (z₁ + iz₂)/(z₁ − iz₂)`, for which `|χ| = 1`, `χ(−z) = χ(z)` and
`χ(R(α) z) = e^{2iα} χ(z)`. Then `χ((𝒦, 1)) = −(1 − i𝒦)/(1 + i𝒦)` and, with
`𝔄(s) ≡ e^{2iα_rot(Ω)}` the chain's rotation phasor,

```
𝔄(s) = −e^{2iθ} (1 + i𝒦)/(1 − i𝒦).
```

Since `𝒦` is a real rational function of `s`, the right side is a rational
all-pass function of `s`.

## One lossless detuned cavity = one Blaschke factor in `s`

A lossless cavity with half-bandwidth `δ > 0` and detuning `Δ` phase-shifts the
sidebands by `α± = 2 arctan(ξ ± Ω/δ)`, `ξ = −Δ/δ`. Each `e^{iα±}` is a
first-order all-pass in `Ω` with its pole at `Ω = Δ − iδ`; multiplying them, the
odd powers cancel and

```
e^{i(α₊ + α₋)} = e^{2iα_rot} = (s − p̄)/(s − p),     p = ω², ω = Δ − iδ.
```

Numerically the *sum* `α₊ + α₋` matches the Blaschke factor to `4e−16`; the
difference misses by `1.6`. A tuned cavity (`Δ = 0`) gives zero rotation.

**Dictionary and realizability.** `ω ↦ ω²` maps the open lower half `Ω`-plane
bijectively onto `ℂ \ [0, ∞)`. So `p` is the pole of a physical lossless filter
cavity iff `p ∉ [0, ∞)`, and then the cavity is unique:

```
ω = √p with Im ω < 0,   δ = −Im ω,   Δ = Re ω,   ξ = −Re ω/δ
```

([`cavity_from_pole`](@ref)). An `n`-cavity chain gives a degree-`n` all-pass
`∏ (s − p̄ⱼ)/(s − pⱼ)` normalized to `𝔄(∞) = 1`, and conversely every such
all-pass with no poles on `[0, ∞)` factors uniquely this way. Filter synthesis
is pole placement. ([`filter_cavity`](@ref) builds the cavity as a state-space
system; cascading them with `*` builds the chain.)

## The synthesis theorem

Let `𝒦(s) = N(s)/D(s)` in lowest terms, `n ≡ max(deg N, deg D)`. The design
condition reads `𝔄(s) = −e^{2iθ} (D + iN)/(D − iN)`.

1. **The constant is fixed at high frequency.** Each Blaschke factor → 1 and
   `𝒦 → 0` as `s → ∞`, so `−e^{2iθ} = 1`: `θ = π/2 (mod π)`, KLMTV Eq. (90).
2. **The cavity count is a degree count.** `𝔄(s) = (D + iN)/(D − iN)`, whose
   poles are the `n` roots of `D − iN`. The number of filter cavities is the
   degree of `𝒦` as a rational function of `Ω²`.
3. **Realizability is automatic.** `D` and `N` are real with no common root, so
   `D − iN` has no real roots; every pole is off `[0, ∞)` and is realized by
   exactly one cavity.
4. **The parameters are the roots**, mapped through [`cavity_from_pole`](@ref).

[`cavities_from_rational`](@ref) implements steps 2–4 for arbitrary `N`, `D`.
This is Purdue & Chen's App. A statement ("match the roots of the polynomials")
with a count, existence and uniqueness.

## The KLMTV case

`𝒦 = Λ⁴/(s² + γ²s)` with `Λ⁴ = 2(I₀/I_SQL)γ⁴`: `N` has degree 0, `D` degree 2,
so **two filter cavities**, and the arm-cavity poles `s(s + γ²)` cancel out of
`𝔄` entirely.

```
𝔄(s) = (s² + γ²s + iΛ⁴)/(s² + γ²s − iΛ⁴),      θ = π/2,
p_{I,II} = (γ²/2)(−1 ± √(1 + 8i·I₀/I_SQL))
```

([`klmtv_filter_poles`](@ref), [`klmtv_cavities`](@ref)). KLMTV's `P ≡ 4Λ⁴/γ⁴`
is `8·I₀/I_SQL`, their `Q = (1 + √(1+P²))/2` is the squared real part of
`√(1 + iP)`, and their Eqs. (89) are `δ = −Im √p`, `ξ = −Re √p/δ` written out;
the trigonometric identities of their App. C reduce to one quadratic formula and
two complex square roots. In units `γ = 1`:

| `I₀/I_SQL` | `ξ_I` | `δ_I/γ` | `ξ_II` | `δ_II/γ` |
|---|---|---|---|---|
| 1/4 | 1.4041853665 | 0.3741199545 | −0.1681173890 | 1.0812267357 |
| 1   | 1.7671310136 | 0.5156013216 | −0.2772297408 | 1.3017526994 |
| 10  | 2.1681653525 | 0.8512981061 | −0.3687198181 | 2.0643323900 |

These agree with KLMTV's closed forms to double precision, which also settles
the factor-of-2 question in Eq. (88).

```@example
using StateSpaceIFO
klmtv_cavities(0.25)
```

`examples/klmtv_fig10.jl` sweeps `(Λ/γ)⁴ = 2I₀/I_SQL` over four decades and
reproduces KLMTV Fig. 10.

## Scope

This is lossless synthesis theory: `𝔄` unimodular *is* the lossless
assumption. Filter-cavity loss makes the Blaschke factors contractive and is not
modelled here.
