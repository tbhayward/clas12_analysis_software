# RGK Fall 2018 6.535 GeV provisional reduced cross sections

Provisional RGK Fall 2018 extraction using the closure-selected one-iteration IBU result and joint data-and-response count-bootstrap covariance.

Each campaign CSV contains only bins passing that artifact's final validity mask.
The reduced cross-section unit is `nb/(GeV^2 rad)`. Phi is recorded in degrees.
The Q2 and xB coordinate columns are the coordinates used to evaluate the virtual-photon flux.
When a beam energy is supplied, epsilon is evaluated at those same Q2 and xB coordinates.
The -t and phi coordinate columns are geometric bin centers; all four bin edges are included.

The error column is the square root of the stored covariance diagonal. It includes
data counting fluctuations, finite response-simulation counting fluctuations, and
finite radiative-correction simulation uncertainty. It does not include detector,
selection, luminosity, current-efficiency, bin-centering-model, or other pending
campaign systematic covariance. The current-efficiency treatment is
`campaign current-efficiency weighting referenced to 35 nA`. Treat these exports as provisional until the
remaining systematics are finalized.

`rgk_fa18_6535_reduced_cross_section_phi_covariance.npz` stores the full within-cell
phi covariance and validity mask for each sample. Use it instead of independent
diagonal errors when fitting harmonics across phi.

## Samples

- `rgk6535`: RGK Fall 2018, 6.535 GeV, torus+1
  - valid bins: 2032
  - beam energy: 6.535 GeV
  - bins with finite physical epsilon: 2032
  - source SHA-256: `f0e7be8a3aa0fe3b60d2c6b42780f5dedf5d29ed37308945c18d4b79363ebb1c`
  - covariance replicas: 300
