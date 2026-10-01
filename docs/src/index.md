```@meta
CurrentModule = StateSpaceIFO
```

# StateSpaceIFO

State-space and chain-scattering tools for the quantum noise of
gravitational-wave interferometers.

The package treats an interferometer as a linear system `(A, B, C, D)` in the
two-photon quadrature picture and computes the things one wants from it with
linear algebra alone: the best achievable strain noise, the readout that
attains it, whether lossless filter cavities can build that readout, and the
cavities themselves. The worked example throughout is the conventional
interferometer with variational readout of Kimble, Levin, Matsko, Thorne and
Vyatchanin (KLMTV, *Phys. Rev. D* **65**, 022002, 2001), which every routine
here reproduces.

It depends only on the `LinearAlgebra` and `Printf` standard libraries.

## Contents

| | |
|---|---|
| [Conventions](@ref) | quadratures, `s = +iΩ`, spectral normalization, the two `J` signatures |
| [The KLMTV plant in state space](@ref) | the 4-state plant, KLMTV Eq. (16), `S_h^min`, the reality defect |
| [Filter cavities as Blaschke factors](@ref) | one lossless detuned cavity = one pole in `s = Ω²`; KLMTV Eq. (89) as a quadratic |
| [Readout design: convexity, H² and H∞](@ref) | convex in the readout coordinate, not in the hardware; one cavity vs two |
| [Chain scattering and J-lossless factorization](@ref) | Kimura's toolbox, validated on his worked examples |
| [Factoring the KLMTV plant](@ref) | KLMTV's cavities from a dual (J,J′)-lossless factorization |
| [Quadratic invariance of a double-passed element](@ref) | an exact QI test, applied to a bidirectional internal filter |
| [API](@ref) | docstrings |

## Quick start

```@example quick
using StateSpaceIFO

# KLMTV's two variational-readout filter cavities at I₀/I_SQL = 1 (units γ = 1)
klmtv_cavities(1.0)
```

```@example quick
# The optimal strain noise of the 4-state plant at Ω = γ, and the readout direction
P = klmtv_plant()
Gn, _, Gh = transfer_matrices(P, 1.0)
Σ = Gn * Gn'
(sh_min(Gh, Σ), optimal_readout(Gh, Σ))
```

## Examples

The `examples/` directory holds runnable scripts that print the numbers quoted
in these pages. Run them with the examples environment:

```
julia --project=examples examples/klmtv_state-space.jl
```

| script | contents | runtime |
|---|---|---|
| `klmtv_state-space.jl` | the plant, `S_h^min`, convexity in `q` vs `(ξ, δ, θ)`, the reality defect | seconds |
| `klmtv_hinf.jl` | one-cavity H² vs H∞ designs by multistart search; feasibility | ~1 min |
| `klmtv_fig10.jl` | reproduces KLMTV Fig. 10 (needs CairoMakie) | seconds |
| `jlossless.jl` | Kimura's worked examples; optical elements; Potapov split; the marginal pole | seconds |
| `klmtv_jspectral.jl` | dual factorization of the regularized KLMTV plant | seconds |
| `qi_test.jl` | quadratic-invariance test of a double-passed internal element | seconds |
