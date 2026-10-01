# Chain scattering and J-lossless factorization

Script: `examples/jlossless.jl`. Equation and theorem numbers are from
H. Kimura, *Chain-Scattering Approach to H∞ Control*, Birkhäuser 1997
(DOI 10.1007/978-0-8176-8331-3).

## Why this framework fits optics

Chain scattering is how optics is already written:

| optics | Kimura |
|---|---|
| scattering matrix of an element (fields in → fields out) | `Σ`, the input/output form (4.1) |
| transfer/ABCD matrix (left- vs right-going fields) | `Θ = CHAIN(Σ)` (4.2), [`chain`](@ref) |
| cascading elements along a beam | product `Θ₁Θ₂` (Fig. 4.3), `*` on [`SS`](@ref) |
| terminating a port (homodyne readout, a load) | `HM(Θ; S)` (4.82), [`hm`](@ref) |
| lossless element | `Θ` is `(J,J′)`-unitary (Lemma 4.4) |
| one detuned cavity = one Blaschke factor | one degree-1 lattice section (5.29) |
| splitting a filter chain into its cavities | Lemma 4.9, Potapov pole splitting, [`potapov_split`](@ref) |
| "inner = cavities, outer = calibration" | `G = ΘΠ`, `Θ` lossless, `Π` unimodular (Def. 6.1) |

Kimura §6.1: `(J,J′)`-lossless factorization `G = ΘΠ` generalizes inner–outer
factorization (the case `r = 0`, `J = I`) and spectral factorization.

## What is implemented

| function | Kimura |
|---|---|
| [`chain`](@ref), [`hm`](@ref), [`dhm`](@ref) | (4.19), (4.82), (4.91) |
| [`ric`](@ref) | stabilizing ARE solution via the ordered Schur form of the Hamiltonian |
| [`jlossless_conjugator`](@ref) | Theorem 5.2, (5.3)–(5.5) |
| [`signature_sqrt`](@ref) | `M = EᵀJE`, (6.18) |
| [`jj_factor`](@ref) | Theorem 6.6, (6.19)–(6.33); subsumes 6.5 |
| [`dual_jj_factor`](@ref) | Theorem 6.12, (6.50)–(6.66) |
| [`jjlossless_residual`](@ref) | Theorem 4.5, (4.55)–(4.56), plus the frequency test (4.50) |
| [`dual_lossless_residual`](@ref) | Lemma 4.11, (4.71), (4.74)–(4.75) |
| [`potapov_split`](@ref) | Lemma 4.9, (4.69)–(4.70) |

One implementation note that is not cosmetic: (6.32) is deliberately
non-minimal (2n states; Kimura reduces Example 6.4 by hand). On the unreduced
realization the Theorem 4.5 Lyapunov equation is singular, because `σ(−Āᵀ)` and
`σ(A + BF)` contain sign-reversed pairs. The residual tests therefore run
[`minreal`](@ref) first.

## Validation against Kimura's worked examples

| example | check | result |
|---|---|---|
| Ex 4.4 | `P` from Theorem 4.5 vs book's `diag(1/4, 1/2)` | matches, `1e−16` |
| Ex 4.6 | Theorem 4.5 residuals; `min eig P` vs book's `diag(1, 1/2)` | `6.6e−16`; `+0.5000` |
| Ex 6.3/6.4 | `X` vs `[1 0; 0 0]`, `X̄` vs `[0 0; 0 1/2]` | exact |
| Ex 6.3/6.4 | `ΘΠ = G`; `Θ` (J,J′)-lossless; `Π` unimodular (poles −2, −1) | `9.1e−16`; `4.5e−15`; yes |
| Ex 6.5 (tall) | `X` vs `[1/14 0; 0 2]`, `X̄` vs `[2 0; 0 0]`; `ΘΠ = G` | exact; `1.1e−15` |

`Θ` and `Π` are unique only up to a constant `J′`-unitary factor (Kimura
p. 130), so Ex 6.4 is compared through the invariant `Π*J′Π` (`2.0e−15`).
These checks are in the test suite.

## Optical elements

**A lossless passive element is inner** (`J = I`, the `r = 0` case). For
[`filter_cavity`](@ref):

```
δ=1.0000  Δ=+0.0000   sup‖G*G − I‖ = 1.9e−15
δ=0.5156  Δ=−0.9111   sup‖G*G − I‖ = 1.2e−15      (KLMTV cavity I)
δ=1.3000  Δ=+2.0000   sup‖G*G − I‖ = 1.4e−15
```

and `max ‖G(iΩ) − e^{−iα_com} R(α_rot)‖ = 5.3e−16`.

**A degenerate squeezer** `diag(e^r, e^{−r})`, `r = 0.8`, is not unitary in
quadratures (`‖SᵀS − I‖ = 4.033`), is symplectic (`‖SᵀΣS − Σ‖ = 0`,
`Σ = [0 1; −1 0]`), and is `J`-unitary with `J = diag(1, −1)` in the `(a, a†)`
basis. See [Conventions](@ref) on the two signatures.

## The KLMTV filter chain, split into its cavities

Building the two-cavity KLMTV chain at `I₀/I_SQL = 1` and recovering cavity I
from the cascade by Potapov pole splitting:

```
two-cavity chain: 4 states;  inner (J = J′ = I): sup‖F*F − I‖ = 2.2e−15
   poles: −1.3018±0.3609i,  −0.5156±0.9111i

Lemma 4.9 split on the cavity-I poles:
   Θ₁Θ₂ = F              : 2.4e−15
   Θ₁ poles              : −0.5156±0.9111i   (= cavity I exactly)
   Θ₁ inner              : sup‖T*T − I‖ = 2.9e−15
```

Lemma 4.9 is the matrix Blaschke/Potapov factorization, so it does for a
quadrature-rotating chain what the scalar argument of
[Filter cavities as Blaschke factors](@ref) does for a scalar all-pass.

## The marginal pole (Lemma 5.5)

Kimura states Theorems 7.4 and 7.7 for plants with no poles or zeros on the
jω-axis, and Lemma 5.5 says why: if `A` has a jω eigenvalue, the conjugation
Riccati (5.3) has no stabilizing solution. The free test mass gives a double
pole at `s = 0`. Applying Theorem 5.2 to the plant pair `(A, [B_a B_h])` with
`J = diag(1, −1, −1)` on `(a₁, a₂, h)`, with a suspension
(`ṗ = Gc₁ − kq − μp`):

| `μ` | `k` | poles on jω | J-lossless conjugator |
|---|---|---|---|
| 0 | 0 | 2 | no |
| 0.1 | 0 | 1 | no |
| 0 | 0.01 | 2 | no |
| 0.01 | 0.01 | 0 | yes |
| 0.2 | 1 | 0 | yes |

Neither repair alone works: damping leaves the rigid-translation pole at
`s = 0`, a spring alone leaves an undamped oscillator pair on the axis. Both are
needed — a damped suspension, or equivalently a band weight that rolls off at DC
folded into the plant before factorizing. So the factorization route cannot use
a free test mass.

## Scope

Kimura's theory is classical. The losslessness conditions (4.55)–(4.56) match
quantum physical-realizability conditions in form, but nothing here checks that
a constructed factor is realizable as a quantum network with vacuum-driven
added channels.
