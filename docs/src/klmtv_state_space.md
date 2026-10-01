# The KLMTV plant in state space

Script: `examples/klmtv_state-space.jl`.

## The plant

[`klmtv_plant`](@ref) builds the dark-port plant with states
`x = (c₁, c₂, q, p)`: two intracavity quadratures, test-mass displacement and
momentum, in the rotating frame linearized about the carrier:

```
ċ₁ = −γc₁              + √(2γ)a₁          A  = [ −γ  0  0   0 ]   B_a = [ √2γ  0  ]
ċ₂ = −γc₂ + Gq         + √(2γ)a₂ + GLh         [  0 −γ  G   0 ]         [  0  √2γ ]
q̇  = p/m                                       [  0  0  0  1/m]         [  0   0  ]
ṗ  = Gc₁                                       [  G  0  0   0 ]         [  0   0  ]
b  = −a + √(2γ)c        C = √(2γ)[I₂ 0₂],  D_a = −I,  B_h = (0, GL, 0, 0)ᵀ
```

Keyword arguments add an intracavity loss rate `γl` with its own vacuum port,
an arm detuning `Δc` (a symplectic rotation of `(c₁, c₂)`), and a suspension
spring `k` and damping `μ` (`ṗ = Gc₁ − kq − μp`).

Three structural facts are visible in the matrices:

1. **The two `G`s are the same number.** `c₁ → ṗ` (radiation pressure) and
   `q → ċ₂` (phase modulation) share a coefficient because both come from
   `H_int = −G c₁ q`. That equality is physical realizability, not a modelling
   choice.
2. **The coupling graph `c₁ → p → q → c₂` is acyclic** — there is no `c₂ → force`
   edge because the cavity is on resonance. This is why the input–output matrix
   is a unit-triangular shear `M = [[1,0],[−𝒦,1]]`. Detuning adds a cycle and an
   optical spring.
3. **Signal and back-action enter at the same node (`q`) by different paths**:
   `h` directly, `a₁` through the cavity pole first. Hence `signal ∝ √(2𝒦)`,
   `back-action ∝ 𝒦`.

## KLMTV Eq. (16) by one matrix inversion

`G_n(s) = D_a + C(sI − A)⁻¹B_a` ([`transfer_matrices`](@ref)) gives

```
G_n = e^{−2iβ}·[[1, 0], [−𝒦(Ω), 1]],   G_h = √(2γ)GL/(s+γ)·(0,1)ᵀ,
e^{−2iβ} = (γ−s)/(γ+s),  β = arctan(Ω/γ),  𝒦 = 2γG²/(mΩ²(γ²+Ω²)).
```

| check | result |
|---|---|
| `𝒦` from `−G_n[2,1]/G_n[1,1]` vs `2γG²/(mΩ²(γ²+Ω²))` | rel. err `5.4e−16` |
| `β` from `−½ arg G_n[1,1]` vs `arctan(Ω/γ)` | abs. err `1.1e−16` |
| `G_n[1,2]`, `G_n[2,2] − G_n[1,1]` | exactly `0` |

## The optimal readout

Measure `y = w(iΩ) b` with `w` a 1×2 row (filter then homodyne). Minimizing the
Rayleigh quotient `S_h = wΣw†/|wG_h|²` over the direction of `w` gives

```
S_h^min(Ω) = 1/(G_h† Σ⁻¹ G_h),     attained at  v ∝ G_h† Σ⁻¹
```

([`sh_min`](@ref), [`optimal_readout`](@ref)). `G_h†Σ⁻¹G_h` is the Fisher
information per unit bandwidth and `v` is inverse-covariance weighting: the
variational readout is a matched filter in quadrature space. Equivalently
`S_h^min` is a Schur complement of `Σ`. With `Σ = MM†` and `G_h ∝ (0,1)` it
reproduces KLMTV Eqs. (58)–(59): `v ∝ (𝒦, 1)`, `S_h^min = h²_SQL/2𝒦`.

| check | result |
|---|---|
| `1/(G_h†Σ⁻¹G_h)` vs `1/|G_h[2]|²` | rel. err `4.1e−11` |
| `v₁/v₂` vs `𝒦` | rel. err `8.2e−16` |

With `w ∝ (𝒦, 1)`, `wM = (0, 1)`: the optimal readout sees only `a₂`.
Back-action is cancelled outright — `a₁` appears in `b₁` undisturbed by the
signal, so it is a measured disturbance subtracted by feedforward. This is exact
disturbance decoupling. A readout with `n_y` independent output quadratures and
`n_s` signals can null at most `n_y − n_s` vacuum channels; a single homodyne has
`n_y = 2`, `n_s = 1`, hence exactly one. That is why variational readout removes
back-action completely and does nothing about loss.

## What lossless passive filters can build: the reality defect

A lossless passive element acts on the sidebands as `a± → r± a±` with
`|r±| = 1`, so any lossless passive chain followed by a homodyne at `θ` produces

```
w = e^{iα_com} · (real unit vector at angle θ − α_rot).
```

The design freedom is the real profile `α_rot(Ω)`, and the natural scalar
variable is the **readout coordinate** `q(Ω) ≡ cot(θ − α_rot) = w₁/w₂`, in which

```
S_h = (h²_SQL/2𝒦)·[(q − 𝒦)² + 1].
```

So a lossless passive chain attains the optimum iff `G_h†Σ⁻¹` is real up to an
overall phase. The diagnostic is the **reality defect**

```
ρ(Ω) ≡ |sin(arg v₁ − arg v₂)| ∈ [0, 1],     ρ = 0 ⟺ realizable
```

([`reality_defect`](@ref)), one 2×2 solve per frequency, computable before any
filter exists. Encoding a direction by its double-angle phasor
`χ(z) = (z₁ + iz₂)/(z₁ − iz₂)` ([`double_angle_phasor`](@ref)),
`|χ(v)| = 1 ⟺ Im(v₁v̄₂) = 0 ⟺ ρ = 0`: the reality defect and "the target must be
a rational all-pass" (see [Filter cavities as Blaschke factors](@ref)) are the
same condition.

[`readout_gap`](@ref) returns `ρ`, `||χ| − 1|`, the best noise over all
directions and the best over real directions ([`sh_min_real`](@ref)):

- **Tuned arm, loss at the input coupler** (`γl/γ` up to 0.1): `ρ < 4e−16` and a
  gap of 0 dB at every frequency tested. KLMTV's lossless two-cavity chain stays
  exactly optimal under this loss. Two reasons are readable off the matrices:
  the loss column is a real multiple of `G_n + I`, so `Σ` stays real; and the
  signal enters one quadrature only.
- **Detuned arm** (`γl = 0.01γ`): `ρ = 0.15` at `(Δc, Ω) = (0.3γ, γ)`, `0.46` at
  `(γ, γ)`, `≳ 0.97` at `Ω = 3γ`; the penalty for a real readout reaches
  **0.589 dB** at `Δc = Ω = γ`. `Σ` is still real; what goes complex is the
  direction of the signal vector `G_h`, because detuning writes strain into both
  quadratures with unequal phase.

When `ρ ≠ 0` no lossless passive filter chain of any order reaches the optimum.
Closing the gap requires either deliberate asymmetric attenuation of the
sidebands (which pays vacuum) or an active, squeezing element in the chain.
