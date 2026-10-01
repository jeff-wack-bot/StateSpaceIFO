# Readout design: convexity, H² and H∞

Scripts: `examples/klmtv_state-space.jl` (section D), `examples/klmtv_hinf.jl`.

## Convex in the readout coordinate, not in the hardware

With the readout coordinate `q(Ω)` of [The KLMTV plant in state space](@ref),
a band-weighted noise objective is

```
J(q) = ∫ μ(Ω) [(q − 𝒦)² + 1] dΩ,     μ ≡ W h²_SQL/2𝒦 > 0,
```

a positive-definite quadratic functional minimized at `q = 𝒦`. The *same
objective* in hardware coordinates `(ξ, δ, θ)` of one filter cavity is badly
non-convex, because `q = cot(θ − α_rot(ξ, δ; Ω))` wraps the variables through
two arctangents and a cotangent.

Measured with the inspiral weight `W ∝ Ω^{−7/3}` on `Ω/γ ∈ [0.1, 10]`,
`I₀/I_SQL = 1`. First, `q` expanded in a nested fixed-pole basis
`1/(Ω²(p_k² + Ω²))` with generic poles (γ deliberately absent, so `𝒦` is in no
finite span): the optimum is one linear solve.

| basis size | min eig of Hessian | `J*/J_opt` | excess [dB] |
|---|---|---|---|
| 1 | `+8.0e+06` | 83.357972 | +19.2095 |
| 2 | `+5.8e−01` | 1.157422 | +0.6349 |
| 3 | `+2.0e−03` | 1.000789 | +0.0034 |
| 4 | `+1.4e−05` | 1.000060 | +0.0003 |
| 6 | `+5.5e−08` | 1.000002 | +0.0000 |

Then one physical cavity, by refined grid search: best
`ξ = +0.992, δ/γ = 0.989, θ/π = 0.497`, `J = 1.4481 J_opt` (**+1.608 dB**). The
script also prints a non-convexity certificate — two points whose midpoint is
worse than the chord:

```
x₁ = (ξ 0.50, δ 0.20, θ/π 0.30)   J₁ = 7.8787e+03
x₂ = (ξ 0.50, δ 5.00, θ/π 0.70)   J₂ = 8.4151e+03
midpoint                          J  = 8.2796e+03 > (J₁ + J₂)/2 = 8.1469e+03
```

Two numbers to keep: one cavity instead of two costs ≈1.6 dB on this
objective, and a convex fit with three generic fixed poles is within 0.003 dB of
optimal. Convexity is a property of the coordinates.

## H² vs H∞

With unit vacuum, `S_h(Ω) = ‖T(iΩ)‖²` for the closed-loop row `T`, so the H∞
norm is the peak of the strain-noise PSD and the weighted H² norm is the
band-integrated weighted noise. Raw `sup S_h` is dominated by band edges; the
usable object is the ratio to a target curve. For KLMTV, measure against the
ideal variational curve:

```
R(Ω) ≡ S_h/S_h^min = (q − 𝒦)² + 1
```

([`excess_ratio`](@ref), for a chain `p = [ξ₁, δ₁, …, θ]`). `R` is
dimensionless, so **H² design = weighted L² approximation of `𝒦` by `q`; H∞
design = unweighted Chebyshev approximation of `𝒦` by `q`**.

**At the exact order the question is empty.** With KLMTV's two cavities,
`max |R − 1| = 0` over `Ω/γ ∈ [0.02, 50]`: every objective is minimized by the
same design. The distinction has content only below the exact order.

**One cavity**, band `Ω/γ ∈ [0.3, 3]`, by multistart search (excess over ideal):

| design | ξ | δ/γ | θ/π | band-average | worst case |
|---|---|---|---|---|---|
| H² (weighted L²) | +1.8367 | 0.5733 | 0.6807 | **+0.072 dB** | +0.871 dB |
| H∞ (Chebyshev) | +1.5047 | 0.6767 | 0.6240 | +0.129 dB | **+0.271 dB** |

H∞ buys 0.60 dB of worst case for 0.058 dB of band-average, with two genuinely
different cavities. The H∞ solution equioscillates as Chebyshev theory predicts:
the three largest maxima of `R` are 1.06428, 1.06417, 1.06159 at
`Ω/γ = 3.00, 0.822, 0.370`.

Widen the band to `Ω/γ ∈ [0.1, 10]` and H² gives `+0.244 / +3.863 dB`, H∞
`+0.252 / +3.798 dB`: the trade nearly vanishes and the H∞ solution does not
equioscillate — a single maximum at the top band edge (`R = 2.40` at
`Ω/γ = 10`, next peak 1.11) carries the objective. Over two decades one cavity
is edge-limited rather than ripple-limited; the band width, not the norm, is the
binding decision.

**Feasibility reading.** `γ* = min_designs sup_band S_h/S_target`, so `γ* ≤ 1`
means the goal is achievable. One cavity on `Ω/γ ∈ [0.3, 3]`, target
`c × (ideal curve)`:

| goal | `γ*` | verdict |
|---|---|---|
| `S_h ≤ 1.0 × ideal` | 1.0643 | infeasible — needs a second cavity |
| `S_h ≤ 1.2 × ideal` | 0.8869 | feasible |
| `S_h ≤ 1.5 × ideal` | 0.7095 | feasible |
| `S_h ≤ 2.0 × ideal` | 0.5321 | feasible |

The search is a coarse sweep plus multistart shrinking-box descent and needs a
warm start (the H² optimum seeds the H∞ search) and a self-check that each
design wins on its own objective. [Factoring the KLMTV plant](@ref) obtains the
`lev* = 1` feasibility boundary as an inertia certificate instead.

## Marginal pole

`𝒦 ∝ 1/Ω²` has a double pole on the imaginary axis from the free test mass, so
`𝒦 ∉ L∞` and the unweighted H² norm diverges. Restricting to a band (as above)
repairs it. The same obstruction reappears as a hard theorem in the
factorization route; see [Chain scattering and J-lossless factorization](@ref).
