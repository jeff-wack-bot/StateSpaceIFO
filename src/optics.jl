# Lossless filter cavities in the two-photon quadrature picture, and the
# pole <-> cavity dictionary.

"""
    filter_cavity(δ, Δ) -> SS

Lossless single-sided detuned filter cavity in the two-photon quadrature
picture, with half-bandwidth `δ` and detuning `Δ`.  State `x = (x₁, x₂)` is the
cavity mode:

    ẋ = [-δ  Δ; -Δ  -δ] x + √(2δ) u
    y = -u + √(2δ) x

The transfer matrix is inner (`G*G = I` on the imaginary axis) and equals
`exp(-iα_com) R(α_rot)`; see [`rotation_angle`](@ref).
"""
function filter_cavity(δ, Δ)
    A = [-δ Δ; -Δ -δ]
    B = sqrt(2δ) * Matrix(1.0I, 2, 2)
    return SS(A, B, sqrt(2δ) * Matrix(1.0I, 2, 2), -Matrix(1.0I, 2, 2))
end

"""
    sideband_phases(ξ, δ, Ω) -> (α₊, α₋)

Phase shifts `α± = 2 atan(ξ ± Ω/δ)` imparted by a lossless filter cavity with
offset `ξ = -Δ/δ` on the upper and lower sidebands at sideband frequency `Ω`.
The factor of 2 corrects KLMTV Eq. (88).
"""
sideband_phases(ξ, δ, Ω) = (2atan(ξ + Ω / δ), 2atan(ξ - Ω / δ))

"""
    rotation_angle(ξ, δ, Ω)

Quadrature rotation `α_rot = (α₊ + α₋)/2` of a lossless filter cavity — even
in `Ω`, and zero for a tuned cavity.
"""
rotation_angle(ξ, δ, Ω) = atan(ξ + Ω / δ) + atan(ξ - Ω / δ)

"""
    common_phase(ξ, δ, Ω)

Common (sideband-averaged) phase `α_com = (α₊ - α₋)/2` of a lossless filter
cavity — odd in `Ω`.
"""
common_phase(ξ, δ, Ω) = atan(ξ + Ω / δ) - atan(ξ - Ω / δ)

"""
    cavity_from_pole(p) -> (ξ, δ, Δ)

The unique lossless filter cavity whose rotation phasor `e^{2iα_rot}` is the
Blaschke factor `(s - p̄)/(s - p)` in `s = Ω²`: `ω = √p` on the branch
`Im ω < 0`, half-bandwidth `δ = -Im ω`, detuning `Δ = Re ω`, offset
`ξ = -Δ/δ`.  Any `p ∉ [0, ∞)` is realizable.
"""
function cavity_from_pole(p)
    ω = sqrt(complex(p))
    imag(ω) > 0 && (ω = -ω)      # physical branch: Im ω < 0
    δ = -imag(ω)
    Δ = real(ω)
    return (ξ = -Δ / δ, δ = δ, Δ = Δ)
end

"""
    klmtv_filter_poles(I0 = 1.0) -> (p_I, p_II)

Poles in `s = Ω²` (units `γ = 1`) of the KLMTV variational-readout all-pass
`(s² + s + iΛ⁴)/(s² + s - iΛ⁴)`:
`p_{I,II} = (-1 ± √(1 + 8i I₀/I_SQL))/2`, with `I0 = I₀/I_SQL`.
"""
function klmtv_filter_poles(I0 = 1.0)
    root = sqrt(1 + 8im * I0)
    return ((-1 + root) / 2, (-1 - root) / 2)
end

"""
    klmtv_cavities(I0 = 1.0) -> Vector of (ξ, δ, Δ)

The two KLMTV variational-readout filter cavities (KLMTV Eqs. (89)) at input
power `I0 = I₀/I_SQL`, in units `γ = 1`, ordered (I, II).
"""
klmtv_cavities(I0 = 1.0) = [cavity_from_pole(p) for p in klmtv_filter_poles(I0)]

"Roots of the polynomial `c[1] + c[2] x + … + c[n+1] xⁿ` via the companion matrix."
function roots_of(c)
    c = collect(float.(complex.(c)))
    while length(c) > 1 && abs(c[end]) < 1.0e-14
        c = c[1:(end - 1)]
    end
    n = length(c) - 1
    n < 1 && return ComplexF64[]
    C = zeros(ComplexF64, n, n)
    for i in 2:n
        C[i, i - 1] = 1
    end
    for i in 1:n
        C[i, n] = -c[i] / c[end]
    end
    return eigvals(C)
end

"""
    cavities_from_rational(N, D) -> Vector of (δ, ξ)

Filter cavities for a real-rational readout coordinate `q = N(x)/D(x)`,
`x = Ω²`, with `N`, `D` given as ascending coefficient vectors.  The filter
chain is the Blaschke product whose poles are the roots of `D - iN`; each
pole gives one cavity via [`cavity_from_pole`](@ref).  Sorted by `δ`.
"""
function cavities_from_rational(N, D)
    nn = max(length(N), length(D))
    c = [(i <= length(D) ? D[i] : 0.0) - im * (i <= length(N) ? N[i] : 0.0) for i in 1:nn]
    out = NamedTuple[]
    for p in roots_of(c)
        cav = cavity_from_pole(p)
        push!(out, (δ = cav.δ, ξ = cav.ξ))
    end
    return sort(out; by = c -> c.δ)
end
