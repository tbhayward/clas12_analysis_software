# Corrected bin volumes and bin centering

The original CSVs and existing test directories are preserved. No ROOT event
files are read by the numerical rerunner.

## Convention

Sangbaek Lee's final April 2025 pass-1 note, Section 5.4 (printed p.126),
and the embedded paper Eq. (10) (PDF p.166) define both imported factors
as denominators:

    Frad = average radiative / average Born
    Fbin = average Born / Born(reference)
    sigma = corrected_yield / (luminosity * bin_volume * Frad * Fbin)

The former production formula multiplied the unchanged imported values. That
was a convention error. The extraction now divides by both and writes a
`cross section correction convention` column (`pass1-divide`). Imported
numbers remain unchanged. Recomputed Fbin uses the same average/reference
convention. Relative factor uncertainties still enter in quadrature.

For a fixed yield, volume and normalization, correcting an OLD multiplication
result changes sigma by `1/(Frad*Fbin)^2`. Replacing Fbin in an already correct
division result changes sigma by `Fbin_old/Fbin_new`. These are distinct steps.
Previous new-factor previews used reciprocal Fbin values and must be regenerated
from the original extraction; do not feed those previews into the repair.

The numerical rerunner refuses to guess the convention for unmarked CSVs.
Declare `--input-correction-convention legacy-multiply` for the original saved
nominal test. Newly generated CSVs carry a marker and are detected automatically.
A conflicting explicit convention is rejected, preventing accidental double repair.
All standalone results hold calibration and normalization fixed; a future full
run recomputes downstream stages in the normal way.

## Numerical changes

- `kinematic_bin_integration.h` replaces the 10x10x10 volume count with the
  exact allowed t length and converged adaptive integration in xB and Q².
  The existing generated-domain cuts are retained: W threshold from the saved
  configuration, 0.19<y<0.85, and t_min. Detector thresholds remain in acceptance.
- `external_scripts/kinematic_domain.py` implements the same domain for the
  standalone numerical stages. C++/Python agreement was checked on all 242
  energy/kinematic boxes in the supplied nominal run; the largest difference
  in allowed fraction was about 6e-12. Custom W thresholds are also tested.
- `km15_bin_center_batch.py` integrates only inside the allowed boundaries.
  It splits changing integration limits and maps the lower t boundary to
  smooth square-root behavior. It never drops a bad model value and then
  silently renormalizes over surviving points.
- Quadrature starts at order 3, compares with order 4, and refines only rows
  that need it. Both successive-factor changes and integrated-volume checks
  must meet the requested relative criterion (default 0.002). This is a
  convergence criterion, not a rigorous bound on unknown physics uncertainty.
- Each worker imports KM15 once. Phi-independent CFF/Fourier coefficients at
  identical physical inputs are reused; initial cross-phi reuses are checked
  against the original model. The propagators and trigonometric dependence
  are still evaluated at each phi. All caches are bounded. The model parameters
  remain fixed within a worker. `KM15_DISABLE_CFF_CACHE=1` restores original
  coefficient calls for verification. Real-model cached/scalar checks matched
  exactly.
- Volume and factor fields retain more numerical precision. No empirical
  normalization or detector efficiency is refit by the standalone rerunner.

The full supplied CSV check has 2,720 populated factor energy-rows. 2,714
converge at the requested criterion; six Sp19 rows near phi=0/360 remain
unconverged. Their saved bin indices are 642, 663, 1136, 1152, 2258 and 2269.
They require further model/integration investigation. The code does not invent
an additional generator cut or label their last estimates as valid factors.

The default is strict: a failed centering check prevents publication of the
combined CSV. The accurate volume-only CSV is still saved by the rerunner.
For an exploratory plot, `--on-unconverged omit` explicitly leaves those Fbin
and cross-section energy-rows blank in the new combined CSV. The corresponding
original and volume-only values remain available. Failed rows and convergence
information are written to numerical_diagnostics.csv. This preview is not a
complete replacement for production results.

## Reference coordinates in plots

The final-by-t plotter supports:

    --reference-kinematics imported
    --reference-kinematics means
    --reference-kinematics midpoints

The default is `imported`, matching the nominal extraction's imported Fbin.
Cross sections use the donor table's common xB/Q²/|t| means and each phi bin's
reference phi for both point positions and exact BH/KM15 comparisons. Dense
curves keep those common xB/Q²/|t| coordinates fixed while sampling phi.
The table is discovered near the input run as `imports/all_bin_v3.csv`, or can
be supplied using `--centering-reference-csv`. Matching uses the same bin index
as the production importer, with geometry and Fbin consistency checks. Missing
or inconsistent donors are errors; there is no silent nearest-bin fallback.
A merged test whose donor rows do not share a common reference needs an
explicit `--reference-kinematics means` approximate comparison.

BSA does not receive the cross-section Fbin correction and retains its measured
means in the default mode. `means` explicitly selects the earlier per-phi
measured-mean diagnostic for both observables, including independently
mean-centered CSVs; `midpoints` retains the legacy geometric-center comparison.
The measured-kinematics diagnostic profiles continue to show the actual data
means rather than imported coordinates.

The plotter writes `theory_reference_points.csv` alongside its figures so every
reference coordinate and mode can be inspected. Cache keys include the actual
coordinates and mode. Imported-reference XS predictions cannot be confused
with earlier measured-mean predictions; compatible BSA predictions are reused.
Exact reference-point pulls are independent of display interpolation.

All numerical rerunner comparisons now use measured means, including the
convention-only baseline and the volume-only control. Pass-1 used xB/Q²/|t|
means integrated over phi; this helper uses each individual row's measured
means, so importing pass-1 factors still does not exactly recenter to pass-2 means. Published Pass-1 overlays are disabled in these diagnostic plots.

## Build and validate

From dvcs_cross_section in the farm ROOT/GCC/Python environment:

```sh
make -j10 dvcs_analysis && python3 test_numerical_corrections.py --real-km15
```

The tests compile the standalone C++ volume routine, compare Python and C++,
check allowed quadrature points, constant/nonconstant-model centering, denominator
correction/volume scaling, legacy conversion and idempotence, default imported XS references and explicit mean references, missing-factor handling, and exact reference-point
model routing. The optional real-model check requires gepard in that Python.
The modified ROOT-dependent C++ units also passed local syntax checks; full
farm linking remains part of the build above.

## Rerun on the existing nominal extraction

This one-line preview command creates a NEW directory, produces both numerical
controls and runs final-by-t plotting. Existing directories are refused.

```sh
python3 run_numerical_corrections.py --input-csv test_runs/proton_t_default_001/output/csvs/dvcs_pass2_analysis.csv --output-dir test_runs/proton_t_default_001/numerical_corrections_002 --workers 10 --input-correction-convention legacy-multiply --on-unconverged omit --plot
```

Omit `--on-unconverged omit` if all requested rows must converge before the
combined CSV can be written. Omit `--plot` to calculate the numerical tables
without model plotting.

Outputs:

- `original.csv`: byte-for-byte snapshot of the input.
- `convention_fixed.csv`: division by imported factors, original volumes retained.
- `volumes_only.csv`: division by imported factors and improved volumes.
- `new_fbin_and_volumes.csv`: improved volumes and converged new factors;
  optional preview omissions are explicit blanks.
- `comparison.csv`: old/new factors, volumes and resulting cross-section ratios.
- `numerical_diagnostics.csv`: convergence status, quadrature order, numerical
  differences, model-at-mean and model-bin-average values.
- `numerical_corrections.json`: input hash, settings and omitted-row list.
- `numerical_comparison.png`: factor and volume-change distributions.
- `plots_convention_fixed/`, `plots_volumes_only/` and `plots_new_fbin_and_volumes/`: final-by-t plots.

After repairing the input convention, the frozen-input cross-section update is exactly

    sigma_new / sigma_old = (Fbin_old / Fbin_new) * (V_old / V_new)

Values and statistical errors are rescaled together. Raw/MC/signal yields,
acceptance, luminosities, Frad, BSA and normalization settings are preserved.
Systematics are not reevaluated. Existing systematic fields cannot be used as
an updated uncertainty budget for the new centering model.

The normal full extraction automatically uses the improved C++ volumes after
recompilation. It continues to import pass-1 Fbin by default; selecting a new
production centering prescription is a separate physics decision. Do not run
the normal cross-section importer over the numerical preview: it deliberately
restores the pass-1 factors.

The rerunner rejects overlapping phi input partitions, including the earlier
unrepaired merged test. Numerical corrections cannot repair its event binning.

For the fastest correction-only check, without new KM15 integrations:

```sh
python3 run_numerical_corrections.py --input-csv test_runs/proton_t_default_001/output/csvs/dvcs_pass2_analysis.csv --output-dir test_runs/proton_t_default_001/convention_fix_001 --input-correction-convention legacy-multiply --convention-only --workers 10 --plot
```

This preserves source bytes, yields, volumes, factors and normalization. It
rescales experimental central values and their saved errors, including normed
columns, using the denominator convention. Systematic budgets are not rederived.
The comparison-table ratios refer to the convention-fixed baseline. `--plot`
is the only part of this fast command that needs model evaluation.

## Independent pass-1 convention check and fresh baseline

The main note also explicitly divides in Eq. (3.54), printed page 50 (PDF
page 51); Eqs. (3.56)-(3.58) define the ratios. Fig. 5.21 (printed p.127,
PDF p.128) has an incorrect inset kinematic label: the caption and bin
number 88 identify xB=.204-.268, Q²=1.912-2.510, |t|=.25-.40. All 57
central vector marker values (19 each for Frad, Fbin, and Frad*Fbin) match
the raw imported CSV columns within 2.9e-10, whereas reciprocals do not.
Production uses the two individual factor columns, not the combined column.

To rebuild, preserve available final-by-t model caches, delete the three
named test directories, run a fresh nominal proton-t extraction without
systematics, and plot it:

```sh
make -j10 dvcs_analysis && python3 rerun_nominal_baseline.py --reset-test-runs --plot
```

This explicitly deletes `test_runs/merged_t_001`, `test_runs/photon_t_001`,
and `test_runs/proton_t_default_001` (all contents). It preserves all found
model cache entries in `test_runs/theory_cache_preserved` first, then seeds
the fresh nominal plotting cache. Only matching existing cache keys are
reused; midpoint curves cannot replace mean-reference predictions.
Unrelated directories are retained. Do not run two reset jobs concurrently;
a lock rejects another reset launcher. A non-reset invocation refuses an
existing nominal destination. No previous cross sections are rescaled.

The nominal analysis retains imported pass-1 Frad and Fbin and divides by
both, uses the improved volumes, and defaults to imported-reference
cross-section theory comparisons (measured means for BSA). It does not select the exploratory new centering prescription.
The final output is in `test_runs/proton_t_default_001/output`.

Pull panels now default to [-10,10], configurable by `--pull-limit 10`.
Outliers are shown as red boundary triangles, with count and maximum
absolute pull annotated. Raw pull values, distributions and statistics
are retained without clipping. Single-cell canvases reserve enough title
space. `python3 test_baseline_refresh.py` verifies this and checks the
reset/cache workflow on disposable fixtures without ROOT inputs.

Percent-difference canvases use fixed y limits of -100% to +100%. Full values
remain in their summaries. KM15 is drawn as a solid line and BH as a dashed
line, joining the cached reference-point predictions without changing pulls.

## Dense displayed theory curves

`--phi-dense 24` controls the number of uniformly spaced model evaluations
used to draw KM15 and BH, including both 0 and 360 degrees. This applies to
all reference modes. In the default imported mode, XS curves use the common
coordinates from the Fbin donor table. In explicit means mode (also the BSA
reference in default mode),
displayed curves hold xB, Q² and |t| fixed at the cell's raw-count-weighted
measured means across phi. Weights sum ep->epg unpolarized raw yields over
the three detector topologies and the appropriate energy's run periods.
These are visual reference curves; exact row means are retained separately
for every pull. If count weights are absent, the plotter reports an explicit
fallback to the unweighted average of available row means.

The displayed arrays have separate cache metadata including their reference
kinematics and sampling count. Existing exact row predictions are preserved
when adding display arrays. Changing `--phi-dense` recomputes the display
samples while retaining cached exact row predictions. Distinct requests are
batched and deduplicated within each cell. No ROOT extraction is needed.

`python3 test_dense_theory_curves.py` checks the sample count, fixed reference,
exact pull predictions, cache reuse, and changes in display sampling density.

`python3 test_imported_theory_references.py` checks donor matching, unchanged
cross-section values/errors, imported phi positions, XS/BSA reference separation,
exact pulls, dense curves, cache isolation/reuse and invalid-input rejection.
No C++ rebuild or ROOT extraction is needed for this plotting change.
