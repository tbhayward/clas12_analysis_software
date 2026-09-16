# External π0 normalization: Hall-A Dlamini 2021 vs AAOgen/VPK

This directory validates the high-W π0 model actually used by the RGA AAOgen
production against the published Hall-A E12-06-114 π0 structure functions.

## Physics convention

For RGA production input `physics model = 5`, `aao_norad/dsigma.F` switches
π0 events with W > 1.9 GeV to `dvmpw`, the March-2021 Valery Kubarovsky
parameterization in `dvmpx.F`.

Dlamini et al. define

d4σ = Γ/(2π) [ σT + εσL
             + sqrt(2ε(1+ε)) σLT cosφ
             + ε σTT cos2φ
             + h sqrt(2ε(1-ε)) σLT' sinφ ].

This is the same φ convention used by `DVMPX`. Therefore the comparison uses
the underlying `XSIGMA_*` functions, not the electron-level `DVMPX` cross
section and not the transformed `DVMPW` quantities used internally by the
generator.

For this VPK implementation:
- XSIGMA_L = 0 exactly
- XSIGMA_U = XSIGMA_T
- XSIGMA_LT' = 0 exactly

The Python code is a literal transcription of the π0 parameters and formulas
from `aao_norad/dvmpx.F`, so it does not require AAOgen or a Fortran compiler
at runtime.

## Run

From:

    /u/home/thayward/clas12_analysis_software/analysis_scripts/dvcs_cross_section/external_pi0_normalization/

run:

    python3 validate_aao_pi0_model.py

Input expected at:

    import/dlamini2021_halla_pi0_structure_functions.csv

Outputs are written to:

    output/dlamini2021_vpk_point_comparison.csv
    output/dlamini2021_vpk_setting_summary.csv
    output/png/

## First checks

The script compares its calculated epsilon and -tmin against the rounded
values published in Dlamini Table I. Large disagreements here mean the
kinematic convention has been implemented incorrectly and the model/data
ratios should not be used.

The principal normalization diagnostic is the per-setting weighted
Dlamini/VPK ratio for dσU/dt. The t' dependence within each setting should
also be inspected before interpreting a setting-level ratio as a pure
normalization difference.

## Source references

- M. Dlamini et al., Phys. Rev. Lett. 127, 152301 (2021),
  arXiv:2011.11125.
- M. Dlamini, PhD thesis, Ohio University (2018), "Measurement of Hard
  Exclusive Electroproduction of Neutral Meson Cross Section in Hall A of
  JLab with CEBAF at 12 GeV."
- AAOgen production source: `aao_norad/dsigma.F`, `dvmpw.F`, `dvmpx.F`.
