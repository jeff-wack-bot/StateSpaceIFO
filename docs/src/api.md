```@meta
CurrentModule = StateSpaceIFO
```

# API

```@index
```

## State space

```@docs
SS
Base.:*(::SS, ::SS)
nx
evalfr
poles
tilde
minreal
isstable
isantistable
Jsig
```

## Chain scattering

```@docs
chain
hm
dhm
ric
jlossless_conjugator
signature_sqrt
jj_factor
dual_jj_factor
jjlossless_residual
dual_lossless_residual
potapov_split
```

## Filter cavities

```@docs
filter_cavity
sideband_phases
rotation_angle
common_phase
cavity_from_pole
klmtv_filter_poles
klmtv_cavities
roots_of
cavities_from_rational
```

## KLMTV plant and readout

```@docs
klmtv_plant
transfer_matrices
klmtv_K
ponderomotive_gain
ponderomotive_gain_ss
backaction_gain
arm_phase
sh_min
optimal_readout
reality_defect
double_angle_phasor
sh_min_real
readout_gap
excess_ratio
dual_chain_plant
dual_factor_design
```

## Quadratic invariance

```@docs
projector
qi_defect
qi_basis
```

## Module

```@docs
StateSpaceIFO
```
