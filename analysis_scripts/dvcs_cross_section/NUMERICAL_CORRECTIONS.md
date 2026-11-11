# Corrected bin volumes and bin centering

The original CSVs and existing test directories are preserved. No ROOT event
files are read by the numerical rerunner.

## Convention

The extraction currently uses

    sigma = corrected_yield * Frad * Fbin / (luminosity * bin_volume)
    stored Fbin = model(at measured means) / model average over allowed bin

A paper may instead define the reciprocal and divide by it. The code retains
its existing multiplicative convention. Changing the integration does not
silently invert either imported or recomputed factors.

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

The final-by-t plotter now supports

    --reference-kinematics midpoints
    --reference-kinematics means

The legacy default stays `midpoints` for existing comparisons. In `means`
mode, data use measured mean phi and theory is evaluated at each row's exact
measured xB/Q²/|t|/phi means. Theory is drawn as model points, because the
reference kinematics vary from one phi row to another. Pulls use those exact
predictions; no interpolation across uncomputed mean points is substituted.
Cache keys include the actual reference coordinates and mode, so old midpoint
curves cannot be reused for this comparison accidentally.

The numerical rerunner uses midpoint comparisons for the volume-only control,
which retains the imported factors, and measured-mean comparisons for the newly
centered result. Published Pass-1 overlays are disabled in these diagnostic plots.

## Build and validate

From dvcs_cross_section in the farm ROOT/GCC/Python environment:

```sh
make -j10 dvcs_analysis && python3 test_numerical_corrections.py --real-km15
```

The tests compile the standalone C++ volume routine, compare Python and C++,
check allowed quadrature points, constant-model centering, multiplicative
correction/volume scaling, missing-factor handling, and exact reference-point
model routing. The optional real-model check requires gepard in that Python.
The modified ROOT-dependent C++ units also passed local syntax checks; full
farm linking remains part of the build above.

## Rerun on the existing nominal extraction

This one-line preview command creates a NEW directory, produces both numerical
controls and runs final-by-t plotting. Existing directories are refused.

```sh
python3 run_numerical_corrections.py --input-csv test_runs/proton_t_default_001/output/csvs/dvcs_pass2_analysis.csv --output-dir test_runs/proton_t_default_001/numerical_corrections_001 --workers 10 --on-unconverged omit --plot
```

Omit `--on-unconverged omit` if all requested rows must converge before the
combined CSV can be written. Omit `--plot` to calculate the numerical tables
without model plotting.

Outputs:

- `original.csv`: byte-for-byte snapshot of the input.
- `volumes_only.csv`: improved volumes with imported Fbin retained.
- `new_fbin_and_volumes.csv`: improved volumes and converged new factors;
  optional preview omissions are explicit blanks.
- `comparison.csv`: old/new factors, volumes and resulting cross-section ratios.
- `numerical_diagnostics.csv`: convergence status, quadrature order, numerical
  differences, model-at-mean and model-bin-average values.
- `numerical_corrections.json`: input hash, settings and omitted-row list.
- `numerical_comparison.png`: factor and volume-change distributions.
- `plots_volumes_only/` and `plots_new_fbin_and_volumes/`: final-by-t plots.

The frozen-input cross-section update is exactly

    sigma_new / sigma_old = (Fbin_new / Fbin_old) * (V_old / V_new)

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
