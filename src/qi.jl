# Quadratic invariance (Rotkowitz & Lall 2006) of pointwise controller
# patterns.
#
# A subspace S of controllers is quadratically invariant under the plant
# G = P22 iff K G K ∈ S for all K ∈ S.  For a pattern that is the same at every
# frequency, write K = Σᵢ cᵢ Eᵢ over a basis of S:
#
#     K G K = Σᵢ cᵢ² (Eᵢ G Eᵢ) + Σ_{i<j} cᵢ cⱼ (Eᵢ G Eⱼ + Eⱼ G Eᵢ),
#
# and since the cᵢ are free, QI holds iff every coefficient matrix lies in S.

"Orthogonal projector onto `span(basis)` in the Frobenius inner product."
function projector(basis::Vector{<:AbstractMatrix})
    V = hcat([vec(ComplexF64.(E)) for E in basis]...)
    Q = Matrix(qr(V).Q)[:, 1:size(V, 2)]
    return M -> reshape(Q * (Q' * vec(ComplexF64.(M))), size(M))
end

"""
    qi_defect(basis, G) -> (defect, witness)

Exact quadratic-invariance test at one frequency.  Returns the largest
relative component of `EᵢGEⱼ + EⱼGEᵢ` (or `EᵢGEᵢ`) lying outside `span(basis)`,
over all basis pairs, and the pair `(i, j)` attaining it.  `0` means QI at
this frequency; `1` means the products lie entirely outside the subspace.
"""
function qi_defect(basis::Vector{<:AbstractMatrix}, G)
    Π = projector(basis)
    n = length(basis)
    worst = 0.0; witness = nothing
    for i in 1:n, j in i:n
        M = i == j ? basis[i] * G * basis[i] :
            basis[i] * G * basis[j] + basis[j] * G * basis[i]
        nm = norm(M)
        nm < 1.0e-12 && continue
        d = norm(M - Π(M)) / nm
        if d > worst
            worst = d; witness = (i, j)
        end
    end
    return worst, witness
end

const _U2 = [[1.0 0; 0 0], [0.0 1; 0 0], [0.0 0; 1 0], [0.0 0; 0 1]]   # 2x2 matrix units
_blk(a, b, c, d) = [a b; c d]

"""
    qi_basis(kind) -> Vector{Matrix}

Bases of controller patterns on a 4×4 matrix of 2×2 quadrature blocks
(`kind = :single` gives a single 2×2 block):

  - `:tied`     — `diag(F, F)`: two passes of one element, tied
  - `:untied`   — `diag(F₁, F₂)`: two independent elements, no backscatter
  - `:antidiag` — `[0 X; Y 0]`: a purely back-reflecting element
  - `:full`     — unconstrained 4×4
  - `:single`   — unconstrained 2×2 (single-pass element)
"""
function qi_basis(kind::Symbol)
    Z2 = zeros(2, 2)
    untied = vcat([_blk(u, Z2, Z2, Z2) for u in _U2], [_blk(Z2, Z2, Z2, u) for u in _U2])
    antidiag = vcat([_blk(Z2, u, Z2, Z2) for u in _U2], [_blk(Z2, Z2, u, Z2) for u in _U2])
    kind === :tied && return [_blk(u, Z2, Z2, u) for u in _U2]
    kind === :untied && return untied
    kind === :antidiag && return antidiag
    kind === :full && return vcat(untied, antidiag)
    kind === :single && return copy(_U2)
    throw(ArgumentError("unknown pattern $kind"))
end
