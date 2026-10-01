# Kimura's chain-scattering / (J,J')-lossless toolbox.
#
# Equation and theorem numbers are from
#   H. Kimura, *Chain-Scattering Approach to H∞ Control*, Birkhäuser 1997:
#
#   CHAIN / DCHAIN                (4.19), (4.38)     scattering -> transfer
#   HM  (homographic transf.)     (4.82), (4.85)     termination = feedback
#   Theorem 4.5                   (4.55), (4.56)     state-space (J,J')-losslessness
#   Lemma 4.9                                        Potapov pole-splitting factorization
#   Theorem 5.2                   (5.3)-(5.6)        J-lossless conjugation
#   Theorem 6.5                   (6.18)-(6.22)      (J,J')-lossless factorization, G stable
#   Theorem 6.6                   (6.29)-(6.33)      general case, two Riccati equations
#   Theorem 6.12                  (6.50)-(6.66)      dual factorization

# --- CHAIN and HM ----------------------------------------------------------

"""
    chain(Σ, m, p) -> SS

`CHAIN(Σ)` in state space, Kimura (4.19).  `Σ` is partitioned with outputs
`a₁` (`m` rows) / `a₂` (rest) and inputs `b₁` (rest) / `b₂` (`p` columns);
requires `D₂₁` invertible.
"""
function chain(Σ::SS, m::Int, p::Int)
    A, B, C, D = Σ.A, Σ.B, Σ.C, Σ.D
    B1, B2 = B[:, 1:(end - p)], B[:, (end - p + 1):end]
    C1, C2 = C[1:m, :], C[(m + 1):end, :]
    D11, D12 = D[1:m, 1:(end - p)], D[1:m, (end - p + 1):end]
    D21, D22 = D[(m + 1):end, 1:(end - p)], D[(m + 1):end, (end - p + 1):end]
    iD21 = inv(D21)
    return SS(
        A - B1 * iD21 * C2,
        [B2 - B1 * iD21 * D22   B1 * iD21],
        [C1 - D11 * iD21 * C2; -iD21 * C2],
        [
            D12 - D11 * iD21 * D22   D11 * iD21;
            -iD21 * D22            iD21
        ]
    )
end

"""
    hm(Θ, S, m, p)

Homographic transformation `HM(Θ; S) = (Θ₁₁S + Θ₁₂)(Θ₂₁S + Θ₂₂)⁻¹`, Kimura (4.82),
on constant matrices (evaluate `Θ` at a frequency first).
"""
function hm(Θ::AbstractMatrix, S::AbstractMatrix, m::Int, p::Int)
    T11, T12 = Θ[1:m, 1:p], Θ[1:m, (p + 1):end]
    T21, T22 = Θ[(m + 1):end, 1:p], Θ[(m + 1):end, (p + 1):end]
    return (T11 * S + T12) / (T21 * S + T22)
end

"""
    dhm(Ψ, S, m)

Dual homographic transformation `DHM(Ψ; S) = -(Ψ₁₁ - SΨ₂₁)⁻¹(Ψ₁₂ - SΨ₂₂)`,
Kimura (4.91), on constant matrices.
"""
function dhm(Ψ::AbstractMatrix, S::AbstractMatrix, m::Int)
    P11, P12 = Ψ[1:m, 1:m], Ψ[1:m, (m + 1):end]
    P21, P22 = Ψ[(m + 1):end, 1:m], Ψ[(m + 1):end, (m + 1):end]
    return -(P11 - S * P21) \ (P12 - S * P22)
end

# --- Riccati ---------------------------------------------------------------

"""
    ric(A, R, Q) -> X or nothing

Stabilizing solution of `XA + AᵀX - XRX + Q = 0`, i.e. `A - RX` stable.
Computed from the stable invariant subspace of the Hamiltonian
`[A -R; -Q -Aᵀ]`.  Returns `nothing` if no stabilizing solution exists
(including when the Hamiltonian has eigenvalues on the imaginary axis).
"""
function ric(A, R, Q)
    n = size(A, 1)
    H = [A  -R; -Q  -transpose(A)]
    any(abs(real(λ)) < 1.0e-9 for λ in eigvals(H)) && return nothing
    F = schur(H)
    sel = real.(F.values) .< 0
    count(sel) == n || return nothing
    ordschur!(F, sel)
    U = F.Z[:, 1:n]
    U1, U2 = U[1:n, :], U[(n + 1):end, :]
    rank(U1) == n || return nothing
    X = U2 / U1
    X = real((X + transpose(X)) / 2)
    isstable(A - R * X) || return nothing
    return X
end

# --- J-lossless conjugation and factorization --------------------------------

"""
    jlossless_conjugator(A, B, J; stabilizing = true) -> (Θ, X)

Kimura Theorem 5.2.  Stabilizing J-lossless conjugator of the pair `(A, B)`:
`X` solves `XA + AᵀX - XBJBᵀX = 0` with `A - BJBᵀX` stable, and
`Θ(s) = [-Aᵀ | XB ; -JBᵀ | I]` (5.5), taking `D_c = I`.  For the
anti-stabilizing conjugator pass `stabilizing = false`.  Returns
`(nothing, nothing)` when the Riccati equation has no suitable solution —
by Lemma 5.5 this always happens when `A` has an eigenvalue on the jω-axis.
"""
function jlossless_conjugator(A, B, J; stabilizing = true)
    n = size(A, 1)
    if stabilizing
        X = ric(A, B * J * transpose(B), zeros(n, n))
    else
        # anti-stabilizing: X s.t. A - BJB'X is anti-stable  <=>  stabilizing for -A
        Y = ric(-A, -B * J * transpose(B), zeros(n, n))
        X = Y === nothing ? nothing : -Y
    end
    X === nothing && return (nothing, nothing)
    Θ = SS(-transpose(A), X * B, -J * transpose(B), Matrix{Float64}(I, size(B, 2), size(B, 2)))
    return (Θ, X)
end

"""
    signature_sqrt(M, m, r) -> E or nothing

Nonsingular `E` with `M = Eᵀ J_{mr} E` (Kimura (6.18)), or `nothing` if the
inertia of the symmetric matrix `M` is not (`m` positive, `r` negative).
"""
function signature_sqrt(M, m, r)
    Ms = (M + transpose(M)) / 2
    F = eigen(Symmetric(Ms))
    λ, V = F.values, F.vectors
    ord = sortperm(λ, rev = true)          # positives first, to match diag(I_m, -I_r)
    λ, V = λ[ord], V[:, ord]
    (count(>(0), λ) == m && count(<(0), λ) == r) || return nothing
    return Diagonal(sqrt.(abs.(λ))) * transpose(V)
end

"""
    jj_factor(G, m, p) -> NamedTuple

`(J_{mr}, J_{pr})`-lossless factorization `G = ΘΠ` of an `(m+r)×(p+r)`
transfer matrix, by Kimura Theorem 6.6 (which subsumes Theorem 6.5 when `G`
is stable, via `X̄ = 0`).  Here `r = size(G.D, 1) - m = size(G.D, 2) - p`.

Returns `(ok, Θ, Π, X, X̄, E, why)`; on failure `ok = false` and `why` names
the condition that failed.
"""
function jj_factor(G::SS, m::Int, p::Int)
    A, B, C, D = G.A, G.B, G.C, G.D
    r = size(D, 1) - m
    r == size(D, 2) - p || error("inconsistent partition")
    J, Jp = Jsig(m, r), Jsig(p, r)
    fail(why) = (
        ok = false, Θ = nothing, Π = nothing, X = nothing,
        X̄ = nothing, E = nothing, why = why,
    )

    # (i)  E with D'JD = E' J' E
    E = signature_sqrt(transpose(D) * J * D, p, r)
    E === nothing && return fail("D'JD has the wrong inertia — (6.18) unsolvable")
    R̃ = transpose(D) * J * D
    S = transpose(C) * J * D
    Q = transpose(C) * J * C

    # (iii) Xbar:  Xbar A' + A Xbar + Xbar C'JC Xbar = 0,  Abar = A + Xbar C'JC stable
    X̄ = ric(transpose(A), -Q, zeros(size(A)))
    X̄ === nothing && return fail("Riccati (6.29) has no stabilizing solution")

    # (ii)  X:  the (6.19) Riccati, Ahat = A + B F stable
    iR = inv(R̃)
    X = ric(A - B * iR * transpose(S), B * iR * transpose(B), Q - S * iR * transpose(S))
    X === nothing && return fail("Riccati (6.19) has no stabilizing solution")

    # (iv)  spectral radius condition
    ρ = maximum(abs.(eigvals(X * X̄)))
    ρ < 1 || return fail(@sprintf("sigma(X Xbar) = %.4f >= 1 — (6.31) violated", ρ))

    F = -iR * (transpose(D) * J * C + transpose(B) * X)       # (6.21)
    Â = A + B * F                                              # (6.20)
    Ā = A + X̄ * Q                                              # (6.30), Q = C'JC
    n = size(A, 1)

    # Theta, (6.32)
    Aθ = [-transpose(Ā)  zeros(n, n); zeros(n, n)  Â]
    Bθ = [I(n)  -X; -X̄  I(n)] \ [transpose(C) * J * D; B]
    Cθ = [-C * X̄   C + D * F]
    Θ = SS(Aθ, Bθ, Cθ, D) * inv(E)

    # Pi, (6.33)
    Π = E * SS(Ā, -(B + X̄ * transpose(C) * J * D), F * inv(I(n) - X̄ * X), Matrix(1.0I, p + r, p + r))

    return (ok = true, Θ = Θ, Π = Π, X = X, X̄ = X̄, E = E, why = "")
end

"""
    dual_jj_factor(G, m, q; strict = true) -> NamedTuple

Dual `(J_{mr}, J_{mq})`-lossless factorization `G = ΩΨ` of an `(m+q)×(m+r)`
transfer matrix, by Kimura Theorem 6.12 (6.50)–(6.66).  `Ω` is unimodular
(`(m+q)` square), `Ψ` is dual (J,J′)-lossless, and
`G J_{mr} G~ = Ω J_{mq} Ω~` (6.9).

Note the asymmetry with the primal [`jj_factor`](@ref): there the two
signatures share `r` (`J_{mr}`, `J_{pr}`); here they share `m` (`J_{mr}`,
`J_{mq}`).  That is why the primal routine cannot be reused by transposition.

`strict = true` enforces Theorem 6.12 as stated, including `Y ≥ 0`, `Ȳ ≥ 0`
and `σ(YȲ) < 1`.  `strict = false` returns the factors whenever (6.50) is
solvable and both Riccati equations have stabilizing solutions, leaving the
verification to the caller.  On the KLMTV plant the semidefiniteness
conditions fail (see the documentation), so `strict = false` is what the
KLMTV analysis uses.

`Ψ` is returned as `Ω⁻¹G` (definitional) rather than via (6.63); the two
disagree on the KLMTV plant and the discrepancy is unresolved, so the
definitional form is used and checked directly.

Returns `(ok, Ω, Ψ, iΩ, Y, Ȳ, E, minY, minȲ, ρ, why)`.
"""
function dual_jj_factor(G::SS, m::Int, q::Int; strict::Bool = true)
    A, B, C, D = G.A, G.B, G.C, G.D
    r = size(D, 2) - m
    size(D, 1) == m + q || error("inconsistent partition")
    J, Jp = Jsig(m, r), Jsig(m, q)
    n = size(A, 1)
    fail(why) = (
        ok = false, Ω = nothing, Ψ = nothing, iΩ = nothing,
        Y = nothing, Ȳ = nothing, E = nothing,
        minY = NaN, minȲ = NaN, ρ = NaN, why = why,
    )

    # (i)  E with D J D' = E J' E'   (6.50).  Same inertia test as the primal,
    #      transposed: signature_sqrt gives F with M = F' Jp F, so E = F'.
    F = signature_sqrt(D * J * transpose(D), m, q)
    F === nothing && return fail("D J D' has the wrong inertia — (6.50) unsolvable")
    E = Matrix(transpose(F))

    R̃ = D * J * transpose(D)
    S = B * J * transpose(D)
    Q0 = B * J * transpose(B)
    iR = inv(R̃)

    # (iii) Ybar:  Ybar A + A' Ybar - Ybar B J B' Ybar = 0,  Abar = A - BJB'Ybar stable  (6.62)
    Ȳ = ric(A, B * J * transpose(B), zeros(n, n))
    Ȳ === nothing && return fail("Riccati (6.62) has no stabilizing solution")

    # (ii)  Y from (6.51), rewritten as a standard ARE in Y with A -> Ahat0'
    Â0 = A - S * iR * C
    Y = ric(transpose(Â0), -transpose(C) * iR * C, S * iR * transpose(S) - Q0)
    Y === nothing && return fail("Riccati (6.51) has no stabilizing solution")

    # (ii)/(iii)/(iv) as literally stated in Theorem 6.12
    minY, minȲ = minimum(eigvals(Symmetric(Y))), minimum(eigvals(Symmetric(Ȳ)))
    ρ = maximum(abs.(eigvals(Y * Ȳ)))
    if strict
        minY > -1.0e-9 || return fail(@sprintf("Y is indefinite (min eig %.3e) — (6.12ii)", minY))
        minȲ > -1.0e-9 || return fail(@sprintf("Ybar is indefinite (min eig %.3e) — (6.12iii)", minȲ))
        ρ < 1 || return fail(@sprintf("sigma(Y Ybar) = %.4f >= 1 — (6.12iv)", ρ))
    end

    L = -(S - Y * transpose(C)) * iR                    # (6.53)
    Ā = A - B * J * transpose(B) * Ȳ                    # (6.62)

    # Omega and its inverse, (6.64) and (6.66)
    Ω = SS(Ā, -((I(n) - Y * Ȳ) \ L), C - D * J * transpose(B) * Ȳ, Matrix(1.0I, m + q, m + q)) * E
    iΩ = inv(E) * SS(
        A + L * C, L, (C - D * J * transpose(B) * Ȳ) / (I(n) - Y * Ȳ),
        Matrix(1.0I, m + q, m + q)
    )

    # Psi definitionally, as Omega^{-1} G (see the docstring).
    Ψ = minreal(iΩ * G)

    return (
        ok = true, Ω = Ω, Ψ = Ψ, iΩ = iΩ, Y = Y, Ȳ = Ȳ, E = E,
        minY = minY, minȲ = minȲ, ρ = ρ, why = "",
    )
end

# --- losslessness tests ------------------------------------------------------

"""
    jjlossless_residual(Θ, m, p; Ωs) -> NamedTuple

Kimura Theorem 4.5 test for `(J_{mr}, J_{pr})`-losslessness, on the minimal
realization of `Θ`.  Returns

  - `dD`    : `‖DᵀJD - J′‖`                                          (4.55)
  - `dLyap` : `‖DᵀJC + BᵀP‖` for `P` solving `PA + AᵀP + CᵀJC = 0`   (4.56)
  - `dFreq` : `max` over the grid `Ωs` of `‖Θ(jω)* J Θ(jω) - J′‖`    (4.50)
  - `minP`  : smallest eigenvalue of `P` (`≥ 0` ⟺ lossless, not merely unitary)
  - `nx`    : number of states after reduction
"""
function jjlossless_residual(Θ0::SS, m::Int, p::Int; Ωs = exp10.(range(-2, 2; length = 200)))
    Θ = minreal(Θ0)
    r = size(Θ.D, 1) - m
    J, Jp = Jsig(m, r), Jsig(p, r)
    dD = norm(transpose(Θ.D) * J * Θ.D - Jp)
    λ = eigvals(Θ.A)
    wellposed = all(abs(λ[i] + λ[j]) > 1.0e-8 for i in eachindex(λ), j in eachindex(λ))
    P = wellposed ? lyap(transpose(Θ.A), transpose(Θ.C) * J * Θ.C) : fill(NaN, size(Θ.A))
    dL = norm(transpose(Θ.D) * J * Θ.C + transpose(Θ.B) * P)
    dF = 0.0
    for Ω in Ωs
        T = evalfr(Θ, im * Ω)
        dF = max(dF, norm(T' * J * T - Jp))
    end
    mP = wellposed ? minimum(real.(eigvals(Symmetric((P + transpose(P)) / 2)))) : NaN
    return (dD = dD, dLyap = dL, dFreq = dF, minP = mP, nx = nx(Θ))
end

"""
    dual_lossless_residual(Ψ, m, q; Ωs) -> NamedTuple

Kimura Lemma 4.11 test for dual `(J_{mr}, J_{mq})`-losslessness:

  - `dD`    : `‖DJ_{mr}Dᵀ - J_{mq}‖`                                  (4.74)
  - `dLyap` : `‖DJBᵀ - CQ‖` with `QAᵀ + AQ - BJBᵀ = 0`                 (4.75)
  - `dFreq`, `dFreqRel` : `max ‖Ψ J_{mr} Ψ~ - J_{mq}‖` over `Ωs`, absolute
    and relative to `max ‖Ψ‖²`                                          (4.71)
  - `minQ`  : smallest eigenvalue of `Q` (`≥ 0` for losslessness)
"""
function dual_lossless_residual(Ψ0::SS, m::Int, q::Int; Ωs = exp10.(range(-2, 2; length = 200)))
    Ψ = minreal(Ψ0)
    r = size(Ψ.D, 2) - m
    J, Jp = Jsig(m, r), Jsig(m, q)
    dD = norm(Ψ.D * J * transpose(Ψ.D) - Jp)
    λ = eigvals(Ψ.A)
    wp = all(abs(λ[i] + λ[j]) > 1.0e-8 for i in eachindex(λ), j in eachindex(λ))
    Q = wp ? lyap(Ψ.A, -Ψ.B * J * transpose(Ψ.B)) : fill(NaN, size(Ψ.A))
    dL = norm(Ψ.D * J * transpose(Ψ.B) - Ψ.C * Q)
    dF = sc = 0.0
    for Ω in Ωs
        T = evalfr(Ψ, im * Ω); Tt = transpose(evalfr(Ψ, -im * Ω))
        dF = max(dF, norm(T * J * Tt - Jp)); sc = max(sc, norm(T)^2)
    end
    mQ = wp ? minimum(real.(eigvals(Symmetric((Q + transpose(Q)) / 2)))) : NaN
    return (dD = dD, dLyap = dL, dFreq = dF, dFreqRel = dF / max(sc, 1), minQ = mQ, nx = nx(Ψ))
end

# --- Potapov pole splitting --------------------------------------------------

"""
    potapov_split(Θ, m, p, inC1) -> (Θ1, Θ2)

Kimura Lemma 4.9: split a (J,J′)-lossless `Θ` into `Θ1` (poles selected by
`inC1(λ) -> Bool`) times `Θ2` (the remaining poles), by reordering the Schur
form of `A` and Schur-complementing the Lyapunov solution `P`.  The Schur
form is assumed to be block diagonal after reordering.
"""
function potapov_split(Θ::SS, m::Int, p::Int, inC1)
    r = size(Θ.D, 1) - m
    J = Jsig(m, r)
    F = schur(Θ.A)
    sel = inC1.(F.values)
    ordschur!(F, sel)
    T = F.Z
    n1 = count(sel)
    A = transpose(T) * Θ.A * T            # block upper-triangular; assume block diagonal
    B = transpose(T) * Θ.B
    C = Θ.C * T
    A1, A2 = A[1:n1, 1:n1], A[(n1 + 1):end, (n1 + 1):end]
    B2 = B[(n1 + 1):end, :]
    C1, C2 = C[:, 1:n1], C[:, (n1 + 1):end]
    P = lyap(transpose(A), transpose(C) * J * C)
    P11, P12 = P[1:n1, 1:n1], P[1:n1, (n1 + 1):end]
    Θ1 = SS(A1, -inv(P11) * transpose(C1) * J, C1, Matrix(1.0I, size(Θ.D, 1), size(Θ.D, 1)))  # (4.69)
    Θ2 = SS(A2, B2, C2 - C1 * inv(P11) * P12, Θ.D)                                            # (4.70)
    return Θ1, Θ2
end
