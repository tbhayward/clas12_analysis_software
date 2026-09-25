# RGA Fall 2018 10.604 GeV provisional reduced cross sections

Provisional RGA Fall 2018 extraction using closure-selected polarity-specific IBU iterations and joint data-and-response count-bootstrap covariance.

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
`unit current weights with the same retained-run selections`. Treat these exports as provisional until the
remaining systematics are finalized.

`rga_fa18_10604_reduced_cross_section_phi_covariance.npz` stores the full within-cell
phi covariance and validity mask for each sample. Use it instead of independent
diagonal errors when fitting harmonics across phi.

## Samples

- `torus_plus1`: RGA Fall 2018 torus+1
  - valid bins: 4595
  - beam energy: 10.604 GeV
  - bins with finite physical epsilon: 4595
  - source SHA-256: `ceac5a59f3f742068d452cadf278a7ca66c6bf7958b3799c6c11185684396b08`
  - covariance replicas: 300
- `torus_minus1`: RGA Fall 2018 torus-1
  - valid bins: 4140
  - beam energy: 10.604 GeV
  - bins with finite physical epsilon: 4140
  - source SHA-256: `47a1fcbd1c9937617fe527b60b9632b0ae1b328cf053873d3d6ecac71068b534`
  - covariance replicas: 300
- `combined`: RGA Fall 2018 combined polarities
  - valid bins: 5293
  - beam energy: 10.604 GeV
  - bins with finite physical epsilon: 5293
  - source SHA-256: `b70a6ec08c524c869185a057f3a2d7b5bc477977ac937238a379a0b0dcedd41d`
  - covariance replicas: 300
