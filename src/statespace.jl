# Minimal real/complex state-space plumbing.  Pure LinearAlgebra.

"""
    SS(A, B, C, D)

A state-space system `G(s) = D + C (sI - A)⁻¹ B`.  The four matrices are promoted
to a common element type on construction.
"""
struct SS{T}
    A::Matrix{T}
    B::Matrix{T}
    C::Matrix{T}
    D::Matrix{T}
end
SS(A, B, C, D) = SS(promote_mat(A, B, C, D)...)
function promote_mat(xs...)
    T = promote_type(map(eltype, xs)...)
    return map(x -> Matrix{T}(x), xs)
end

"Number of states."
nx(G::SS) = size(G.A, 1)

"""
    evalfr(G, s)

Evaluate the transfer matrix at the complex frequency `s`.  A plain matrix is
treated as a constant system.
"""
evalfr(G::SS, s) = G.D + G.C * ((s * I - G.A) \ G.B)
evalfr(D::AbstractMatrix, s) = D

"Eigenvalues of `A`."
poles(G::SS) = eigvals(G.A)

"""
    G1 * G2

Series connection, `G2` first and then `G1` (Kimura's concatenation rule (2.13)).
Multiplying by a constant matrix on either side scales the input or output.
"""
function Base.:*(G1::SS, G2::SS)
    n1, n2 = nx(G1), nx(G2)
    A = [G1.A  G1.B * G2.C; zeros(n2, n1)  G2.A]
    B = [G1.B * G2.D; G2.B]
    C = [G1.C  G1.D * G2.C]
    D = G1.D * G2.D
    return SS(A, B, C, D)
end
Base.:*(G::SS, M::AbstractMatrix) = SS(G.A, G.B * M, G.C, G.D * M)
Base.:*(M::AbstractMatrix, G::SS) = SS(G.A, G.B, M * G.C, M * G.D)

"Para-Hermitian conjugate `G~(s) = G(-s)ᵀ` (real coefficients)."
tilde(G::SS) = SS(-transpose(G.A), -transpose(G.C), transpose(G.B), transpose(G.D))

"""
    minreal(G; tol = 1e-9) -> SS

Minimal realization by Kalman decomposition: restrict to the controllable
subspace (range of `[B AB … Aⁿ⁻¹B]`), then to the observable one.  Kimura's
(6.32) deliberately produces a non-minimal 2n-state `Θ`, and the Theorem 4.5
Lyapunov equation is singular on the unreduced realization, so this is not
cosmetic.
"""
function minreal(G::SS; tol = 1.0e-9)
    A, B, C, D = G.A, G.B, G.C, G.D
    for pass in 1:2
        n = size(A, 1)
        n == 0 && break
        K = copy(B)                                  # controllability matrix
        M = copy(B)
        for _ in 1:(n - 1)
            M = A * M; K = [K M]
        end
        F = svd(K)
        k = count(>(tol * max(1, F.S[1])), F.S)
        if k < n
            T = F.U[:, 1:k]
            A, B, C = T' * A * T, T' * B, C * T
        end
        # dualize and repeat to remove unobservable modes
        A, B, C = transpose(A), transpose(C), transpose(B)
    end
    return SS(Matrix(A), Matrix(B), Matrix(C), D)
end

"All eigenvalues of `A` have real part `< -tol`."
isstable(A; tol = 1.0e-9) = all(real(λ) < -tol for λ in eigvals(A))
"All eigenvalues of `A` have real part `> tol`."
isantistable(A; tol = 1.0e-9) = all(real(λ) > tol for λ in eigvals(A))

"Signature matrix `J_{mr} = diag(I_m, -I_r)`, Kimura (4.24)."
Jsig(m, r) = Matrix(Diagonal([ones(m); -ones(r)]))
