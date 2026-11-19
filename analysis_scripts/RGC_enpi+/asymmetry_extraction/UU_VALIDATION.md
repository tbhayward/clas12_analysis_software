# Unpolarized-modulation validation

This diagnostic imports `extract_structure_function_ratios.py` and uses its nominal
seven-parameter, three-period likelihood and minimization. It writes a separate
`asymmetry_extraction/output/unpolarized_validation/` directory. It does not alter
the production fitter, nominal results, dilution factors, or channel-selection cuts.
All generated figures are PNG.

## Run on ifarm

From the parent `RGC_enpi+` analysis directory:

```bash
python3 asymmetry_extraction/validate_unpolarized_modulations.py --workers 8 --replicas 100
```

This runs all three studies in all 24 bins:

1. Two one-dimensional profile likelihoods per bin for
   `u1 = F_UU^cos(phi)/F_UU` and `u2 = F_UU^cos(2phi)/F_UU`, reoptimizing the other
   coefficients and all three dilution nuisances. It compares the profile intervals
   with the ordinary Hesse errors. Each constrained point uses a warm start and an
   independent nominal start. A minimum below the nominal fit is reported as a
   failure requiring inspection, rather than silently redefining the reference.
2. A joint u1/u2 profile grid per bin, profiling all other parameters.
3. Conditional injection/recovery at the fitted inputs and under the four input
   variations described below, retaining the actual selected event counts and
   measured kinematics in each period.

The first two items together are selected by `--studies profiles`; the injection
studies are selected by `--studies closure`. Both are enabled by default.

The script requires the same Python dependencies as the production extraction.
It uses the production default run-info, cut and dilution files, and corrected ROOT
inputs. Valid selected-event staging is reused; if staging is missing, it is rebuilt
from those existing ROOT inputs. No PID skim, momentum calibration, channel fit,
or dilution-factor extraction is rerun. If existing staging is stale, the provenance
check stops: use `--rebuild-events` to rebuild selected-event staging from the
existing ROOT files. Input overrides are supported via `--input PERIOD=FILE`,
`--run-info-csv`, `--cut-json`, and `--dilution-json`.

For a smaller first pass:

```bash
python3 asymmetry_extraction/validate_unpolarized_modulations.py --bins 2 8 14 20 6 24 --workers 8 --replicas 100 --output-dir asymmetry_extraction/output/unpolarized_validation_pilot
```

To run only profiles, with a separate output directory:

```bash
python3 asymmetry_extraction/validate_unpolarized_modulations.py --studies profiles --workers 8 --output-dir asymmetry_extraction/output/unpolarized_profiles
```

Use the same command to resume an interrupted run. Checkpoints are saved after
each one-dimensional profile, each contour row, each group of ten replicas, and
each completed scenario. The identity includes source code, ROOT input identities,
run information, cut and dilution hashes, runtime versions, and study settings.
Changed inputs or settings require a new output directory. Changing only the
number of workers is allowed. No stale fits are accepted as a new configuration.

## Injection scenarios

- `nominal`: all seven coefficients at the newly refitted nominal solution.
- `zero_uu`: u1=u2=0, with the polarized coefficients unchanged.
- `opposite_uu`: signs of both nominal UU coefficients reversed, with the polarized
  coefficients unchanged.
- `half_polarized`: nominal UU coefficients, with all five polarized coefficients
  multiplied by 0.5.
- `zero_polarized`: nominal UU coefficients, with all polarized coefficients zero.

Four scenarios are refitted with the production fitter for every replica. The
exactly zero-polarization scenario instead verifies that the production likelihood
is flat in UU at zero polarization, on both the original events and generated
samples. This is a structurally unidentifiable limit: Gaussian UU uncertainties,
pulls, and coverage are deliberately left undefined. It is not counted as an
ordinary failed fit. It should demonstrate the loss of information, not manufacture
a UU measurement from a numerical minimizer.

Run/helicity labels are regenerated **jointly**, with probabilities proportional to
`Q(run,h) * g(x;run,h,period)`. The actual run charges, run target polarizations,
period beam polarizations, and period dilution factors are retained. Helicity is
not assigned 50/50, and the original run label is not frozen. Kinematics and period
assignments remain fixed. The sampler uses the production cross-section factor.
It checks positivity across all allowed run/helicity states using the extrema of
the target polarization. A nonphysical input scenario is explicitly recorded and
skipped, never clipped or silently shrunk. A nonphysical nominal injection is a
bin failure. Skips and failures are listed in `run_status.json`.

The true dilution values are the fitted nominal nuisances. By default each replica
also draws new Gaussian auxiliary dilution measurements around those true values,
using the production statistical dilution errors, and profiles the nuisances in
its fit. Gaussian draws are not truncated/redrawn; inadmissible draws are recorded
as failures. This supports repeated-data-and-auxiliary-measurement coverage tests.
`--auxiliary-mode fixed_truth` instead fixes each auxiliary center at its truth;
that is a conditional diagnostic, not the same unconditional coverage claim.

Toy fits use the ordinary production multistart procedure and are not seeded with
the injected truth. Fits with invalid/inaccurate/nonpositive covariance or parameters
at limits are reported explicitly and excluded from the ordinary pull summaries.
The summary also reports coverage with **all attempts** in the denominator so
failed fits cannot disappear from the assessment.

## Inspect these outputs first

Under `asymmetry_extraction/output/unpolarized_validation/`:

- `tables/nominal_fits.csv`: fresh nominal solutions and fit-quality flags.
- `tables/unpolarized_correlations.csv`: u1/u2 correlations with the other physics
  coefficients and dilution nuisances.
- `plots/profiles_bin_XX.png` and `tables/profile_intervals.csv`: real-data profiles,
  68.3%/95% crossings, and profile/Hesse error ratios. Missing crossings remain
  blank. Failed scan points are not bridged when finding crossings.
- `plots/joint_bin_XX.png`: profiled u1/u2 contours. An open contour can mean the
  scan is too narrow, a parameter bound is reached, or fits fail; inspect the
  checkpoint grid rather than interpreting the plot alone.
- `plots/closure_summary_u1.png`, `plots/closure_summary_u2.png`, and
  `tables/closure_summary.csv`: bias, empirical spread, pull mean/RMS, coverage,
  uncertainty changes, and all failure counts for each bin and scenario. The
  `sigma_ratio_to_nominal` column compares uncertainty magnitudes with the nominal
  injection scenario.
- `tables/toy_results.csv`: individual replica estimates, errors and pulls.
- `tables/zero_polarization_identifiability.json` and the bin checkpoints: explicit
  flatness checks in the zero-signal case.
- `manifest.json`: inputs, settings, code identities and interpretation.
- `checkpoints/bin_XX.json`: complete scan points, fitted covariances, injected
  truths, auxiliary measurements and reproducible seeds.
- `run_status.json`: bin failures and nonphysical injection scenarios.

Profile thresholds are asymptotic: 1 and 3.84 in 1D, 2.30 and 5.99 in 2D. Toy
coverage here tests the quoted **Hesse** intervals (estimate +/- 1 sigma and
+/- 1.96 sigma). It does not calibrate the profile intervals. With 100 usable toys,
the nominal binomial sampling uncertainties on 68.3% and 95% coverage are roughly
0.047 and 0.022 per bin/scenario. Those fluctuations are not evidence of failure.
Markedly non-Gaussian profiles or poor Hesse coverage would motivate a follow-up
profile-interval coverage study rather than declaring the coefficient validated.

These tests validate the conditional estimator under its assumed model at the
accepted kinematics. They do not establish robustness to omitted harmonics,
kinematic dependence within a bin, or spin-dependent detector acceptance. Those
remain separate questions before deciding to release the UU coefficients.

## Implementation checks

```bash
python3 -m unittest discover -s asymmetry_extraction -p 'test_validate_unpolarized_modulations.py' -v
```

The tests compare generated run/helicity frequencies with explicit probabilities,
compare the production NLL with an explicit state-by-state normalization, check
exact zero-polarization cancellation, reject nonphysical injections, verify the
auxiliary draw distribution, prevent crossings across failed points, expose failed
fits in coverage denominators, and exercise production-fit recovery with a
nonuniform accepted phi distribution. These are implementation tests, not the
results of the full 24-bin physics validation.
