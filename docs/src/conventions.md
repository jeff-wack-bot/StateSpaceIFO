# Conventions

These are the conventions the code in this package uses. Where they differ from
a source, the difference is stated.

## Quadratures and ordering

Two-photon quadrature amplitudes at sideband frequency `Ω`, collected into
column vectors: `a = (a₁, a₂)ᵀ` is the dark-port input (amplitude, phase),
`b = (b₁, b₂)ᵀ` the dark-port output, and

```
R(α) = [cos α  −sin α;
        sin α   cos α]
```

## `s = +iΩ`

Every transfer function is evaluated at `s = +iΩ`. One consequence: phases that
KLMTV write as `e^{+2iβ}` appear here as `e^{−2iβ}`. For the arm cavity,

```
e^{−2iβ} = (γ − s)/(γ + s),     β = arctan(Ω/γ),
```

and in a lossless detuned filter cavity the common phase appears as
`e^{−iα_com}`. Scalar phases cancel from every spectral density, so this never
changes a noise number.

## Rates

`γ` (and `δ` for a filter cavity) is a **half-bandwidth**: the amplitude decay
rate, so the input coupling is `√(2γ)`. KLMTV's offset of a filter cavity is
`ξ = −Δ/δ` with `Δ = ω_res − ω₀` the detuning. Unless stated otherwise the
examples use units `γ = 1`, and set the optomechanical coupling `G`, the mass `m`
and the strain calibration `L` to 1.

## Spectral normalization

Vacuum inputs have unit, uncorrelated spectra (single-sided, KLMTV Eqs. (25)–(26)).
The output noise covariance at one frequency is `Σ = Σ_ports G_port G_port†`, and
the strain-referred noise of a readout row `w` is

```
S_h = (w Σ w†)/|w G_h|².
```

Because `S_h` is invariant under `w → λ(Ω) w` — a common factor on the
photocurrent is offline post-processing — only the readout *direction* is a
design variable.

## Back-action gain

The KLMTV back-action gain and arm phase are read off the noise transfer matrix,

```
𝒦 = −G_n[2,1]/G_n[1,1]          β = −½ arg G_n[1,1]
```

([`backaction_gain`](@ref), [`arm_phase`](@ref)). This is a ratio of
transfer-function entries, so it does not depend on quadrature units or
normalization constants. For the free-mass plant it equals
`𝒦 = 2γG²/(mΩ²(γ² + Ω²))` ([`klmtv_K`](@ref)), i.e. `Λ⁴ = 2γG²/m`. With a free
test mass the mechanical response is `x̃ = −F̃/(mΩ²)`, which is the minus sign in
`M = [[1, 0], [−𝒦, 1]]` with `𝒦 > 0`.

## Filter-cavity phases

A lossless filter cavity shifts the two sidebands by

```
α± = 2 arctan(ξ ± Ω/δ)
```

([`sideband_phases`](@ref)). The factor of 2 corrects KLMTV Eq. (88), which prints
a bare `arctan`; with the factor, the reconstruction reproduces KLMTV's own
Eq. (89) to machine precision (Purdue & Chen, *Phys. Rev. D* **66**, 122004,
2002, App. A note the same erratum). The pair recombines into a common phase and
a quadrature rotation,

```
α_rot = (α₊ + α₋)/2   (even in Ω, zero for a tuned cavity)
α_com = (α₊ − α₋)/2   (odd in Ω)
G(iΩ) = e^{−iα_com} R(α_rot)
```

([`rotation_angle`](@ref), [`common_phase`](@ref)). The rotation is the *sum*;
this is checked to `5e−16` against the state-space model of the cavity.

## Two signatures

Two different indefinite metrics appear in the chain-scattering code, for
unrelated reasons, and they must not be confused:

- **Directional**, `J_{mr} = diag(I_m, −I_r)` ([`Jsig`](@ref)). It comes from
  splitting port variables by direction of travel and appears in every
  chain-scattering representation. After level normalization it is the metric in
  which a performance specification is stated.
- **Bogoliubov**, `J = diag(1, −1)` in the `(a, a†)` basis. A degenerate squeezer
  `diag(e^r, e^{−r})` is symplectic in the quadrature basis but not unitary;
  in the `(a, a†)` basis it is `J`-unitary with this `J`.

A passive lossless element is *inner*: `J = I`, the `r = 0` case. A chain
containing parametric gain carries both signatures, and choosing the wrong one is
a silent error — the Riccati equations still solve and the factor is `J`-unitary
for the wrong `J`. The code in this package only factors chains of passive
elements and the KLMTV plant, for which the directional signature is the
relevant one.
