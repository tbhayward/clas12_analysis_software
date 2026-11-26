# Current predictions: cumulative BM models, 2026-11-25

The main shared plots now use Brandon’s `chiral_split_asym_photon_pred.zip`:

- Orange: longitudinal/pole contributions only (`step1`).
- Red: longitudinal/pole plus chiral-odd contributions (`step1_2`).
- Blue: full BM model including chiral-even transverse contributions (`step1_2_3`).

Each curve is a complete prediction with its own numerator and denominator. They are not additive contributions. The full stage exactly matches all 21 earlier species/observable central-value grids. Shaded theory bands use the supplied lower/upper error widths. These are theory fit uncertainties, not experimental systematic uncertainties.

Predictions remain virtual-photon referenced. UL/LL data retain laboratory-axis labels. Beam-referenced predictions have not yet been supplied. GK remains light gray where available, without annotations on plots where it is unavailable.

`figures/shared_plots/bm/` and `figures/shared_plots/bm+gk/` each contain seven π+ cumulative comparisons, π− LU, and the charge-comparison plot with π0, plus a nine-page PDF. Run the usual `run_comparisons.py` entry point to regenerate them. The PDF filename is retained for compatibility.

Older common-denominator CFF-product plots remain in the detailed collection. The previous shared set is preserved in `archive/fixed_denominator_comparison_20261125/`. Their colored pieces must not be interpreted as the new cumulative models.

The reference material below describes those earlier product-decomposition views.

---

# Pion theory comparisons

Open **index.html** in the Desktop bundle. The first four figures are the π+ UL sinφ and UL sin2φ CFF comparisons, with an all-products and a leading-products version of each.

## Where to find the figures

| Folder | Contents |
|---|---|
| `figures/pi_plus/featured_stories/` | The main π+ UL sinφ and UL sin2φ 2×2 figures |
| `figures/pi_plus/total_comparisons/` | π+ data with updated BM totals |
| `figures/pi_plus/cff_breakdowns/` | π+ CFF comparisons for all seven ratios |
| `figures/pi_plus/BM_and_GK/` | π+ comparisons retaining both theory totals |
| `figures/pi_plus/alternative_error_display/` | π+ statistical and point-to-point uncertainties combined in quadrature |
| `figures/pi_plus/reference_original_layouts/` | Remade versions of the original π+ plotting layouts |
| `figures/pi_minus/` | The corresponding π− total and CFF comparisons |
| `figures/charge_comparison/` | π+ and π− LU together, plus the optional π0 theory and pointwise-residual versions |

Within `cff_breakdowns`, **all_products** shows every supplied CFF product and the supplied `other` term. **leading_products** shows the same five leading products across all four xB panels, with every omitted contribution summed into a gray remainder. Both reproduce the same total to CSV rounding precision.

`four_xB_panels` contains focused 2×2 figures for a single ratio. `seven_ratio_canvases` contains the larger canvases for a single xB range. Filenames identify the channel, ratio, decomposition and binning.

## The two main π+ stories

The four featured files are:

- `pi_plus_UL_sin_phi_all_CFF_products_2x2.png`
- `pi_plus_UL_sin_phi_leading_CFF_products_2x2.png`
- `pi_plus_UL_sin_2phi_all_CFF_products_2x2.png`
- `pi_plus_UL_sin_2phi_leading_CFF_products_2x2.png`

The simpler leading-products plots use Brandon's familiar blue, orange, green, red and purple palette, with a thick black total. For UL sinφ, the red curve is H_T × Etilde_L(pole), orange is Etilde_L(pole) × calE, green is Ebar_T × Htilde_L, blue is F_1,1 × Htilde_L and purple is H_T × Etilde_L(non-pole). For UL sin2φ, F_1,1 × G_1,4 is blue, the pole × calE term is orange, H_T × Etilde_L(pole) is red, F_1,1 × F_1,4 is purple and H_T × calE is green.

Colors remain fixed across xB and between the two pion channels for a given observable. Additional products receive separate solid colors rather than the previous global dash scheme. The labels refer to the new individual-product CSVs. They are not relabeled as the older grouped A/P curves.

## Data and uncertainty conventions

π+ uses the final 24-bin result table. Its default figures show statistical error bars and separate point-to-point systematic boxes. Correlated scale uncertainties are excluded.

π− LU uses the six original −t′ bins. Its UL/LL ratios use the two-bin refit within each xB interval, 0.05–0.65 and 0.65–1.25 GeV². Its fixed UU inputs are labeled and are not plotted as measurements. π− has statistical errors only.

The supplied CFF tables contain central values without updated theory bands. The old bands are not reused. All supplied products plus `other` close to the supplied totals within five-decimal rounding precision. No signs or depolarization factors are changed. The π0 tables are preserved for reference, without inventing π0 data.

UL/LL data retain their lab target-axis definition. No target-axis conversion or experimental-bin folding has been applied to the supplied theory. The faint red PRELIMINARY watermark is 90% transparent. Figures are PNG. The shared collection also includes an explicitly requested nine-page PDF.

## Regenerate or use updated extraction results

In the Desktop bundle, run this from its parent directory:

```bash
python3 /Users/tbhayward/Desktop/pion_theory_comparisons/run_comparisons.py
```

From the analysis repository root:

```bash
python3 analysis_scripts/pion_theory_comparisons/run_comparisons.py
```

The repository command writes to `analysis_scripts/pion_theory_comparisons/output/`, containing the same `figures/`, `data/diagnostics/` and `index.html` organization.

To use new π− fits after reprocessing, provide both the fine and coarse extraction directories. The latter must include `mean_tprime.json`:

```bash
python3 analysis_scripts/pion_theory_comparisons/run_comparisons.py --pim-fine analysis_scripts/RGC_eppi-/asymmetry_extraction/output_fixed_UU --pim-coarse analysis_scripts/RGC_eppi-/asymmetry_extraction/output_two_tprime_bins --output analysis_scripts/pion_theory_comparisons/output_new_campaign
```

Use `--pip-results /path/to/structure_function_ratios.csv` for a replacement π+ table. Use `--theory /path/to/theory` for replacement `pip`, `pim` and `pi0` CSV subdirectories. The measurement reader checks that the π− fine and coarse fits use the same selected sample and polarization metadata.

For a quick style preview:

```bash
python3 analysis_scripts/pion_theory_comparisons/plot_comparisons.py --quick --output analysis_scripts/pion_theory_comparisons/preview
```

The default watermark alpha is 0.10 and resolution is 300 dpi. `--channels pip` or `--channels pim` limits the channel canvases. `--skip-legacy` skips the reference versions of the old π+ layouts. The existing π+ `plot_kroll_comparisons.py` entry point still defaults to the new CFF totals.

## Supporting files

The Desktop bundle keeps the theory CSVs in `data/theory`, normalized measured results in `data/measurements`, and numerical audits in `data/diagnostics`. Plotting code is under `scripts`.

The repository keeps its source snapshots in `inputs` and `measurements`. Its plotting routine produces the same organized figure tree. The archived prediction tables and their uncertainties remain available.

`data/diagnostics/cff_colors.json` defines the color key. `manifest.json` lists every generated figure. The other audit tables record exact input hashes, decomposition closure and differences from the previously stored model predictions. Those earlier predictions may predate the immediate baseline in the email.

`data/model_UU_inputs` contains optional 24-bin UU inputs for a later fit variation. These are model values interpolated at the plotted mean −t′, not bin-folded predictions. They have not been applied to the current fits.

The old Desktop layout is archived outside this cleaned bundle. No fits or numerical results were changed, and no Git commit was created.

## Axis and line-width update

Fixed ranges are UU cosφ: −0.4 to 0.8, UU cos2φ: −0.8 to 0.6, LU sinφ: −0.2 to 0.3 for both channels, UL sinφ: −0.3 to 0.3, UL sin2φ: −0.4 to 0.1, LL: −0.2 to 1.0, and LL cosφ: −0.4 to 0.8. Major ticks include the endpoints, spaced by 0.1 for sine observables and 0.2 for the others. CFF curves were reduced slightly in thickness. These are display changes only.

## Chiral-sector comparison

The first two featured π+ figures now show three known sectors, with unclassified `other` kept separate. Blue contains products involving at least one added chiral-even transverse F/G CFF, including the calE/calEt partners defined by the supplied F11/G14 kernels and interference with older sectors. Red contains products involving H_T or Ebar_T but none of those added CFFs. Orange contains only the longstanding Htilde_L, Etilde_L and pion-pole CFFs. Ordinary chiral-even longitudinal GPDs are not counted as new. Gray is the supplied `other` column, whose components are unknown.

Every subtotal uses the full supplied BM denominator. None is presented as a GK prediction or as a recomputed old-model ratio. A genuine amplitudes-off comparison would also require recalculating the denominator and separating the contents of `other`. The individual-product sector plots use the same colors. Exact assignments are in `data/diagnostics/chiral_sector_classification.csv`.

## All-seven sector plots and CFF contribution ranking

All seven π+ ratios now have featured 2×2 chiral-sector canvases. Equivalent π− canvases remain under its CFF breakdown folder. Legend wording is BM Total, Chiral-even transverse CFFs + interf., Chiral-odd CFFs + interf., and Longitudinal / pion pole products. The unclassified curve and the footer are hidden in every sector plot. Its values remain in the total and numerical audit.

The contribution-size study is in `figures/cff_sensitivity` and `data/diagnostics/cff_sensitivity`. It integrates the absolute net contribution of explicitly listed F14 or G11 products over 0.05–1.20 GeV², sums the four xB settings equally, and divides by the corresponding integral of the absolute total. A separate gross-product share is also tabulated. This avoids singular ratios where the total crosses zero. It is not a complete CFF sensitivity calculation because products hidden in other and the denominator response are unavailable. Entries not separately listed are not asserted to have zero sensitivity. The optional error-weighted comparison is a descriptive pointwise proxy, not a significance.

Shared collections are in `figures/shared_plots/bm/` and `figures/shared_plots/bm+gk/`. Each contains seven π+ plots, π− LU sinφ, and an automatically rebuilt nine-page PDF. GK curves are light gray. The supplied GK table covers only π+ LU, UL and LL ratios. The other three plots retain BM curves without a GK curve or missing-prediction annotation.

The final page of each shared PDF is the π+/π− LU comparison with the π0 theory curve, also copied into each shared folder as a PNG.

Additional theory-comparison display: `figures/shared_plots/bm+gk_combined_errors/` contains the same nine-page set with one data error bar equal to sqrt(stat² + point-to-point systematic²), without separate systematic boxes. The charge-comparison final page uses this convention too. Correlated scale errors are not included. LU sinφ uses −0.2 to 0.3 with 0.1 tick spacing.

## π− conference release update

Use `--pim-release PATH` on `run_comparisons.py` to refresh from the validated
preliminary release folder. It uses six LU bins and three bins for each UL ratio
per xB, preserving the pi+ source table. No pi-minus LL measurements are included.
Shared collections now contain 13 PNG pages and `pion_chiral_sector_comparisons_13_pages.pdf`:
seven pi+ observables, pi-minus LU and both UL observables, the LU charged-pion/π0
comparison, and two charged-pion UL comparisons. π− uses no GK overlays.
The π− UL sin2φ panels use [-1,1] so the new positive measurements remain visible.
