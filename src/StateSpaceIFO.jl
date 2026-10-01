"""
    StateSpaceIFO

State-space and chain-scattering tools for the quantum noise of
gravitational-wave interferometers: a minimal state-space type, Kimura's
(J,J′)-lossless factorizations, the lossless filter-cavity / Blaschke-pole
dictionary, the 4-state KLMTV plant with its optimal-readout quantities, and an
exact quadratic-invariance test.
"""
module StateSpaceIFO

using LinearAlgebra
using Printf

include("statespace.jl")
include("chainscattering.jl")
include("optics.jl")
include("klmtv.jl")
include("qi.jl")

# state space
export SS, nx, evalfr, poles, tilde, minreal, isstable, isantistable, Jsig
# chain scattering
export chain, hm, dhm, ric, jlossless_conjugator, signature_sqrt
export jj_factor, dual_jj_factor, jjlossless_residual, dual_lossless_residual, potapov_split
# filter cavities
export filter_cavity, sideband_phases, rotation_angle, common_phase
export cavity_from_pole, klmtv_filter_poles, klmtv_cavities, roots_of, cavities_from_rational
# KLMTV plant and readout
export klmtv_plant, transfer_matrices, klmtv_K, ponderomotive_gain, ponderomotive_gain_ss
export backaction_gain, arm_phase
export sh_min, optimal_readout, reality_defect, double_angle_phasor, sh_min_real, readout_gap
export excess_ratio, dual_chain_plant, dual_factor_design
# quadratic invariance
export projector, qi_defect, qi_basis

end
