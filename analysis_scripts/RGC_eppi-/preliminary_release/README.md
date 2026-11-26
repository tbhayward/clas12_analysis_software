# Preliminary π− release

The release reports six −t′ bins per xB for LU sinφ and three for each UL ratio.
Both UU cosine ratios stay fixed zero. All five polarized amplitudes float in
the nominal fits. LL and LLcosφ are nuisance amplitudes and are not plotted.
No BM-fixed result is substituted for a free-fit measurement.

## Farm command

After committing/pushing development and pulling on the farm, run inside tmux:

```bash
python3 analysis_scripts/RGC_eppi-/preliminary_release/run_release.py --source /volatile/clas12/thayward/eppi-/loose/full_20261124/analysis --workers 12
```

This reuses the completed nominal six-bin fits and `three_bin_LLcos_study/free`
combined fits. It adds 36 three-bin, individual-RGC-period baseline fits. It does
not repeat channel selection, dilution replicas, or BM-fixed combined fits.
Inputs, cut identities, selected-file identities, metadata and fit validity are
checked. Invalid period fits are recorded rather than silently plotted.

Full output stays at `SOURCE/preliminary_release/`. The compact review folder
and archive appear in the home checkout at `analysis_scripts/RGC_eppi-/`:

- `review_preliminary_release/`
- `review_preliminary_release.tar.gz`

Copy the archive to the desktop and extract it. It includes all release figures,
48 measurement rows, bin kinematics, period comparisons, fit JSONs, selection
and dilution metadata/figures, generated LaTeX tables, and Overleaf figure assets.
It excludes selected ROOTs, calibration caches and replica ensembles.
The output includes no LL or LLcosφ result plots. Original full-study diagnostic
plots remain in their existing source directories, separate from this release.

`--skip-periods` is for a clearly marked local plotting preview only. The final
Overleaf builder requires the full command's period-check completion flag.
`--study PATH`, `--output PATH`, `--review-output PATH`, `--theory PATH` and
`--dpi` can override defaults. The default figure resolution is 220 dpi. Output
figures are PNG only. The requested combined theory-comparison PDFs are built by
the separate theory-comparison routine.

## Outputs

- `measurements.csv`: 24 LU and 12+12 UL measurements in the theory-plot schema.
- `bin_kinematics.csv`: pooled measured means and population standard deviations.
- `release_manifest.json`: input signatures and selected sample counts.
- `release_fits/`, `period_fits/`: exact combined and individual-period fit results.
- `period_results.csv`, `period_fit_summary.json`: period-check results and validity.
- `figures/data_*`: three-observable all-xB, individual-xB, overlay and 2x2 plots.
- `figures/BM_total_*`, `figures/BM_stages_*`: full and cumulative BM comparisons.
- `figures/periods_*`: fine LU and three-bin UL period comparisons.
- `overleaf_additions/`: figure assets and generated result/count/dilution tables.
- `support/`: current selection fits, frozen windows, kinematics and dilution products.

Only statistical uncertainties are shown. The UL errors include first-order
dilution propagation. Point-to-point systematic values are recorded as zero in
the normalized plotting schema to denote unassigned errors, not an assessment
that the systematics vanish. Photon-axis BM target predictions are overlaid
directly on lab-axis measurements as requested, without frame transformation.

## Finalize the Overleaf source locally

Use the exact supplied original π− ZIP as the template, after the full release
archive is available:

```bash
python3 analysis_scripts/RGC_eppi-/preliminary_release/build_overleaf.py --template ORIGINAL_NOTE.zip --release REVIEW_RELEASE_DIRECTORY --output OUTPUT_DIRECTORY
```

The builder preserves unrelated sections and bibliography, replaces current
selection/kinematic figures, updates the abstract, dual binning, selected counts,
dilution statistics, three-observable results, generated numerical tables and
BM comparisons, and checks literal figure references. Compile and visually
review the whole document before delivering the final ZIP. This scripted merge
is preparation for that final editorial/render review, not a substitute for it.

## Refresh the shared pion comparison collection

The comparison runner accepts `--pim-release REVIEW_RELEASE_DIRECTORY` and keeps
the original pi+ source values/systematics unchanged. It replaces all pi-minus
measurement inputs with this release table. The regenerated collection includes
LU, UL sinφ and UL sin2φ charged-pion comparisons, individual π− UL panels, and
13-page shared PDFs containing the original seven π+ panels, three π− panels,
the LU charge/π0 comparison and the two UL charge comparisons.
