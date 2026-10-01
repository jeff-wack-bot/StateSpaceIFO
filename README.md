# StateSpaceIFO

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jeff-wack-bot.github.io/StateSpaceIFO/dev/)
[![Build Status](https://github.com/jeff-wack-bot/StateSpaceIFO/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/jeff-wack-bot/StateSpaceIFO/actions/workflows/CI.yml?query=branch%3Amain)

State-space and chain-scattering tools for the quantum noise of
gravitational-wave interferometers, in pure `LinearAlgebra`:

- a minimal state-space type with series connection, para-Hermitian conjugate
  and Kalman minimal realization;
- Kimura's chain-scattering toolbox — `CHAIN`, `HM`/`DHM`, J-lossless
  conjugation, primal and dual (J,J′)-lossless factorization, losslessness
  tests and Potapov pole splitting — validated on Kimura's worked examples;
- lossless filter cavities and the pole ↔ cavity dictionary in `s = Ω²`;
- the 4-state KLMTV plant, `S_h^min = 1/(G_h†Σ⁻¹G_h)`, the optimal readout and
  its reality defect;
- an exact quadratic-invariance test for pointwise controller patterns.

```julia
using StateSpaceIFO
klmtv_cavities(1.0)   # KLMTV Eq. (89): the two variational-readout filter cavities
```

The `examples/` directory reproduces the numbers in the documentation:

```
julia --project=examples examples/klmtv_state-space.jl
```
