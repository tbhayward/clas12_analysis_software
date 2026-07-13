# Pass-1 systematic interpolation and systematic projection plots

## Interpolation

Exact pass-1/pass-2 bin matches are copied unchanged. A pass-2 bin without a
pass-1 match is filled by averaging the four component systematic uncertainties
from the nearest kinematic source below and above it. The search uses bin-center
coordinates in xB, Q2, |t| and phi, scaled by representative analysis-bin
spacings. At an outer phase-space edge where one side does not exist, the two
nearest source bins overall are used.

The point-to-point total imported at this stage is recalculated in quadrature
from the four interpolated component values.

An audit file is written beside the pass-2 CSV:

    <pass2 csv>.interpolated_pass1_systematics.csv

It records the target bin, both source bins, scaled distances and all assigned
values.

## Projection plots

`make_systematic_projection_plots()` is called by `main_systematics` after the
systematic columns have been populated. It makes PNG-only two-panel canvases for
xB, Q2, |t|, phi, electron theta, proton theta and photon theta.

The upper panel shows the arithmetic mean absolute uncertainty in nb/GeV^4.
The lower panel shows the arithmetic mean of 100*s_i/|sigma_i|. Curves are made
for pi0 subtraction, acceptance, Frad, Fbin, exclusivity cuts, fiducial cuts and
the total point-to-point uncertainty.

Outputs are written to:

    output/systematics/point_to_point_projections/

A CSV containing the plotted averages accompanies every PNG.
