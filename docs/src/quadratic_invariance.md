# Quadratic invariance of a double-passed element

Script: `examples/qi_test.jl`.

## The test

Rotkowitz & Lall (*IEEE TAC* **51**, 274, 2006): a subspace `S` of controllers
is *quadratically invariant* (QI) under the plant `G = P₂₂` iff `KGK ∈ S` for
all `K ∈ S`. QI implies the constraint `K ∈ S` convexifies in the Youla
parameter; Lessard & Lall proved the converse; Zheng, Furieri, Kamgarpour &
Papachristodoulou (*IEEE TAC* 2020, arXiv:1907.06256) showed the same condition
is necessary and sufficient for SLS and IOP.

For a **pointwise pattern** (the same structure at every frequency), write
`K = Σᵢ cᵢEᵢ` over a basis of `S`:

```
K G K = Σᵢ cᵢ² (Eᵢ G Eᵢ) + Σ_{i<j} cᵢcⱼ (Eᵢ G Eⱼ + Eⱼ G Eᵢ).
```

Since the `cᵢ` are free, QI holds iff every coefficient matrix lies in `S` —
an exact, finite, deterministic test ([`qi_defect`](@ref)). The defect is the
largest relative part of those matrices falling outside `S`
([`projector`](@ref)): 0 means QI, 1 means the products lie entirely outside.
[`qi_basis`](@ref) supplies the patterns used below.

## The geometry

Unfold a cavity containing an internal filter element that the field traverses
twice, in opposite directions. Absent backscatter the two traversals are
dynamically uncoupled:

```
y₁ = field arriving at the element from the readout side
u₁ = field leaving it toward the arm                     u₁ = F y₁   (forward)
y₂ = field arriving from the arm side
u₂ = field leaving it toward the readout                 u₂ = F y₂   (backward)
```

so the controller is block diagonal with equal blocks, `K = diag(F, F)`. The
surrounding interferometer connects each controller output to the *other*
input — `y₂ = A u₁` (arm reflection, all-pass × ponderomotive shear) and
`y₁ = R u₂` (signal-extraction mirror) — so the plant is purely block
anti-diagonal,

```
G = [ 0  R ]
    [ A  0 ]
```

and `KGK = [[0, FRF], [FAF, 0]]` is block anti-diagonal for every `K ∈ S`,
meeting the block-diagonal `S` only at the origin.

## What the test returns

```
case                                          defect    verdict
-- positive controls --
unconstrained S, full plant                   0.000     QI
out-of-loop readout (G = P22 = 0)             0.000     QI
block-diagonal S, block-diagonal plant        0.000     QI
single-pass internal element                  0.000     QI

-- bidirectional cases --
double-passed element, TIED  diag(F,F)        1.000     NOT QI (maximal)
double-passed element, UNTIED diag(F1,F2)     1.000     NOT QI (maximal)
...with the SRM removed (R = 0)               1.000     NOT QI (maximal)
...with the arm removed (A = 0, unphysical)   1.000     NOT QI (maximal)
purely back-reflecting element [0 X; Y 0]     0.000     QI
```

The positive controls show the test discriminates. The out-of-loop case
(`P₂₂ = 0`, `K·0·K = 0 ∈ S`) recovers that the readout layer of an
out-of-loop filter is convex. For the double-passed element the projection of
`KGK` onto `S` is identically zero, not merely small.

Readings:

1. **The tie is not the culprit.** Untying the two passes gives defect 1.000,
   identical to the tied case. What breaks QI is bidirectionality itself.
2. Removing the signal-extraction mirror, or the arm, does not help; any one
   coupling path suffices.
3. **Single-passing the element restores QI**, because then there is no block
   structure to violate.
4. In chain-scattering terms the element appears twice in the unfolded cascade,
   `Θ_F · Θ_arm · Θ_F`, so the closed loop is quadratic rather than affine in the
   design variable. Non-QI and "the factor appears twice" are the same
   statement.

## Scope

The test covers a tie that is a genuine subspace (two blocks constrained
equal). A multiplicative tie, such as an OPA pair with `G_f G_b = 1`, is not a
subspace, and QI theory does not apply to it.
