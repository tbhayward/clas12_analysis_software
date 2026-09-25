# FA18 provisional reduced cross sections

Generated UTC: 2026-09-24T17:07:53Z

This package contains provisional reduced cross sections for:

- RGA Fall 2018 at 10.604 GeV: torus+1, torus-1, and their BLUE combination.
- RGK Fall 2018 at 6.535 GeV: torus+1.

Each CSV contains only bins passing the corresponding final validity mask.
Virtual-photon epsilon is evaluated at the stored flux Q2 and xB coordinates.

The uncertainty column includes data counting fluctuations, finite response-GEMC counting fluctuations, and finite radiative-GEMC uncertainty. Remaining detector, selection, luminosity, current-efficiency-parameter, bin-centering-model, and other campaign systematic covariance is not included.

The RGA results use unit current weights with the retained-run selections. The RGK result uses the campaign current-efficiency correction referenced to 35 nA.

The NPZ files contain within-cell phi covariance. Use those covariances for harmonic fits rather than treating CSV errors as independent across phi.

The combined RGA sample already contains the torus+1 and torus-1 information. Do not treat the combined sample and its component polarities as independent measurements.

No matching or interpolation between the RGA and RGK kinematic grids has been applied. These files provide inputs for a Rosenbluth sensitivity study; they are not themselves an L/T separation.

See the README and export_summary.json files inside each energy directory for sample-specific provenance.
