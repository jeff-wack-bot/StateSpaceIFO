# qi_test.jl — quadratic invariance of the bidirectional tie.
#
# Run:  julia --project=examples examples/qi_test.jl
#
# ---------------------------------------------------------------------------
# THE GEOMETRY.  Unfolding the signal-extraction cavity, an internal filter
# cavity is traversed twice, in opposite directions:
#
#     y1 = field arriving at the FC from the readout side
#     u1 = field leaving the FC toward the arm            u1 = F y1  (forward)
#     y2 = field arriving at the FC from the arm side
#     u2 = field leaving the FC toward the readout        u2 = F y2  (backward)
#
# so the controller is BLOCK DIAGONAL with EQUAL blocks,
#     K = diag(F, F),
# while the plant closes the loop by connecting each output to the OTHER input,
#     y2 = A u1   (arm reflection: ponderomotive shear + all-pass)
#     y1 = R u2   (SRM reflection)
# so
#     G = [ 0  R ; A  0 ]     — purely BLOCK ANTI-DIAGONAL.
#
# The constraint structure and the plant structure are exact complements.  That
# is the whole result; the code below confirms it and maps out which variations
# of the topology change the answer.
# ---------------------------------------------------------------------------

using StateSpaceIFO
using LinearAlgebra, Printf

# ======================================================================
#  The physical blocks
# ======================================================================

const Z2 = zeros(2, 2)

"Arm reflection: all-pass times the ponderomotive shear, KLMTV Eq. (16)."
function arm(Ω; γ = 1.0, 𝒦₀ = 2.0)
    𝒦 = 𝒦₀ / (Ω^2 * (γ^2 + Ω^2))
    φ = (γ - im * Ω) / (γ + im * Ω)
    return φ * [1.0 0; -𝒦 1]
end

"Signal-extraction mirror reflection (constant, |r| < 1)."
srm(; r = 0.9, ϕ = 0.3) = r * exp(im * ϕ) * Matrix(1.0I, 2, 2)

"Plant G = P22 for the double-passed internal element."
plantG(Ω; withSRM = true, withArm = true) =
    [
    Z2                                 (withSRM ? srm() : Z2);
    (withArm ? arm(Ω) : Z2)            Z2
]

# ======================================================================
#  Run
# ======================================================================

hdr(s) = (println(); println("="^78); println(s); println("="^78))
verdict(d) = d < 1.0e-10 ? "QI      " : (d > 1 - 1.0e-10 ? "NOT QI (maximal)" : "NOT QI")

const ΩS = [0.05, 0.2, 0.7, 1.0, 2.5, 8.0, 30.0]

function run_case(name, basis, Gf; note = "")
    d = 0.0
    for Ω in ΩS
        dd, _ = qi_defect(basis, Gf(Ω))
        d = max(d, dd)
    end
    @printf("  %-44s  %-8.3f  %-17s %s\n", name, d, verdict(d), note)
    return d
end

hdr("QI test: does K G K stay in S for all K in S?")
println("  `defect` = largest relative part of the products falling OUTSIDE S.")
println("  0 => QI (convex in Youla / SLS / IOP);  1 => products entirely outside S.\n")
@printf("  %-44s  %-8s  %-17s %s\n", "case", "defect", "verdict", "")
println("  " * "-"^92)

println("\n  -- positive controls (the test must be able to return QI) --")
run_case(
    "unconstrained S, full plant", qi_basis(:full), Ω -> plantG(Ω);
    note = "trivially QI"
)
run_case(
    "KLMTV out-of-loop readout (G = P22 = 0)", qi_basis(:tied), Ω -> zeros(4, 4);
    note = "out-of-loop readout layer is convex"
)
run_case(
    "block-diagonal S, block-diagonal plant", qi_basis(:untied),
    Ω -> [arm(Ω) Z2; Z2 srm()]; note = "classic decoupled case"
)
run_case(
    "single-pass internal element", qi_basis(:single), Ω -> srm() * arm(Ω);
    note = "no block structure at all"
)

println("\n  -- the bidirectional cases --")
d_tied = run_case(
    "double-passed FC, TIED  diag(F,F)", qi_basis(:tied),
    Ω -> plantG(Ω)
)
d_untied = run_case(
    "double-passed FC, UNTIED diag(F1,F2)", qi_basis(:untied),
    Ω -> plantG(Ω)
)
run_case(
    "...with the SRM removed (R = 0)", qi_basis(:untied),
    Ω -> plantG(Ω; withSRM = false)
)
run_case(
    "...with the arm removed (A = 0, unphysical)", qi_basis(:untied),
    Ω -> plantG(Ω; withArm = false)
)
run_case(
    "purely back-reflecting element  [0 X; Y 0]", qi_basis(:antidiag),
    Ω -> plantG(Ω); note = "a curiosity, not a filter"
)

# ======================================================================
#  Why: the structure, exhibited
# ======================================================================

hdr("Why it fails, and how badly")

let Ω = 1.0, G = plantG(Ω), Π = projector(qi_basis(:tied))
    F = [0.6 + 0.2im 0.1; -0.3 0.8 - 0.4im]        # an arbitrary K in S
    K = [F Z2; Z2 F]
    M = K * G * K
    @printf("  K = diag(F,F) in S;   K G K  has block structure\n")
    @printf(
        "     [ %s  %s ]\n", norm(M[1:2, 1:2]) < 1.0e-12 ? "   0   " : "nonzero",
        norm(M[1:2, 3:4]) < 1.0e-12 ? "   0   " : "nonzero"
    )
    @printf(
        "     [ %s  %s ]\n", norm(M[3:4, 1:2]) < 1.0e-12 ? "   0   " : "nonzero",
        norm(M[3:4, 3:4]) < 1.0e-12 ? "   0   " : "nonzero"
    )
    @printf(
        "\n  ||projection of K G K onto S||  =  %.3e     (relative %.3e)\n",
        norm(Π(M)), norm(Π(M)) / norm(M)
    )
    println(
        """
        So the projection is not merely small — it is ZERO.  S is block diagonal and
        G is block anti-diagonal, so K G K is block anti-diagonal for every K in S and
        meets S only at the origin.  The failure is maximal, not marginal: there is no
        sense in which the constraint is 'nearly' QI, and no perturbation of the
        parameters repairs it."""
    )
end

hdr("What the variations say")
println(
    """
    1. The tie is NOT the culprit.  Untying the two passes (letting the forward and
       backward traversals be different elements) leaves the defect at $(round(d_untied; digits = 3)) —
       identical to the tied case, $(round(d_tied; digits = 3)).  What breaks QI is the *bidirectionality*
       itself: an internal element that does not backscatter is block diagonal,
       while the surrounding interferometer does nothing except couple the two
       directions.  The parameter tie is an extra constraint layered on a set that
       was already non-QI.

    2. Removing the SRM does not help, and neither would removing the arm.  Any one
       coupling path suffices to push K G K off the diagonal.

    3. SINGLE-PASSING the internal element restores QI outright, because then there
       is no block structure to violate: the controller is one block and the plant
       is the composite round trip.

    4. The chain-scattering reading of the same fact: in the unfolded cascade the
       filter appears TWICE in the product,  Theta_FC * Theta_arm * Theta_FC, so
       the closed loop is QUADRATIC in the design variable rather than affine.

    5. Scope.  This applies to a filter-cavity tie, which is a genuine subspace
       (two blocks constrained to be equal).  It does NOT cover a multiplicative
       tie such as an OPA pair with G_f * G_b = 1, which is not a subspace at all —
       QI theory does not apply to it."""
)

println("\ndone.")
