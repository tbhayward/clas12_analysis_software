# Manley–Bhattacharya exclusive pi+ predictions

The CSV and `predictions.tex` are preserved exactly as supplied by Brandon Manley and Shohini
Bhattacharya on 2026-11-19. Source hashes are in
`provenance.json`. The source document is retained for the theory definitions;
it is not an input to the analysis-note build.

The CSV contains 16 -t-prime points per xB group (0.02–1.20 GeV^2), for
(xB,Q2/GeV^2) = (0.2,2.2), (0.3,2.6), (0.4,3.1), (0.5,4.3).
Its seven `A_*` columns are already structure-function ratios normalized to
F_UU,T + epsilon F_UU,L, with depolarization factors divided out. They use Trento
angles and target polarization along the virtual photon. Do not multiply by
additional depolarization factors or change signs when reading these columns.

`_lo` and `_hi` are nonnegative asymmetric **error lengths**, not endpoints.
The one-sigma fit bands are [central - _lo, central + _hi]. As described in the
supplied note, these bands come from the baseline fit to CLAS6 pi0 unpolarized
cross sections; they are not an additional experimental systematic uncertainty.

`plot_kroll_comparisons.py` maps A_UU_cosphi -> u1, A_UU_cos2phi -> u2,
A_LU_sinphi -> lu1, A_UL_sinphi -> ul1,
A_UL_sin2phi -> ul2, A_LL -> ll0, A_LL_cosphi -> ll1. It keeps the supplied
kinematics and t-prime points. It does not average the curves over experimental
bins or extrapolate this CSV to zero. Each curve is compared to the matching
experimental xB group. The two theories have slightly different supplied Q2
values; joint-plot legends show both. Existing measurement/frame definitions
are retained without a target-axis conversion.

From the RGC_enpi+ directory, regenerate all 14 PNGs using the regular results:

```bash
python3 asymmetry_extraction/plot_kroll_comparisons.py --input-csv asymmetry_extraction/output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv --output-dir asymmetry_extraction/output/kroll_comparisons
```

The five existing Kroll PNG names stay unchanged in the output root. Five
Manley–Bhattacharya PNGs (four individual groups plus all-xB) are written in
`bhattacharya_manley/`. Four individual-xB joint plots are written in
`both_theories/`; there is no all-xB joint plot. In joint plots, Kroll is black
solid and Manley–Bhattacharya is purple dashed with its one-sigma shaded band.
The original data colors/markers and systematic boxes are retained. New-theory
canvases extend the standard vertical range only when necessary to contain the
supplied theory band. All outputs remain PNG; `--final` suppresses the existing
PRELIMINARY watermark. `--prediction-csv` accepts an alternate supplied CSV.

The updated pi+, pi0 and pi- CSVs are stored here. Only pi+ is compared to the
RGC exclusive pi+ measurements; pi0/pi- tables are retained for future use.
Manley–Bhattacharya and joint canvases now have seven panels: UU on the top row,
single-spin on the middle row and double-spin on the bottom row. Joint plots
show only Manley–Bhattacharya predictions in the UU panels, because no UU GK
input was supplied. GK-only canvases retain their original five-panel layout.
