# Exact production-count pi0 response comparison

```sh
python3 test_production_pi0_response.py && python3 study_production_pi0_response.py --production-run "$PWD/test_runs/current_mc_fixed_001/proton_t_default_001" --response-followup "$PWD/test_runs/particle_response_add_pt_followup_001" --workers 10 --output "$PWD/test_runs/production_pi0_response_001"
```

Use the actual completed baseline run path if yours differs. Requires the baseline
CSV, saved combined/global cuts and saved current-response calibration JSON.
Response parameters are hash-checked against the completed add-pT follow-up.
It uses the full-statistics production loader, not sampled diagnostic caches.

The driver copies/instruments total_counts.cpp in its output backend only; it
never changes the production counter or binary. New headers and the applied source
remain reviewable in that directory. No make/rebuild of dvcs_analysis is needed.
The standalone backend compiles automatically. Ten workers total by default.

Only pi0 reconstructed MC and pi0-as-DVCS leakage MC are counted. DATA, generated
MC, DVCS MC, theory and systematics are not recalculated. Baseline DATA columns
are retained exactly. Each relevant MC channel is scanned once per scenario.
Nominal must pass before morph is run; no extra control weights are used.

The nominal gate compares raw and current-corrected MC triples for all topologies
and period totals in saved valid production rows. Terminal invalid rows are not
used to validate physics yields. A mismatch aborts with nominal_mismatches.csv
and nominal_gate.json. This can expose changed input files, processing code,
selection/calibration provenance or count conventions; do not bypass it.

Both scenarios use original production MC cut windows, full globals and the exact
production xB,Q²,|t|,phi matcher/current-weight code. The morphed scenario changes
four-vectors and derived coordinates before all those steps. Stored phi offsets
are preserved while its physical Trento-angle change is applied. The original
period current model is used without requiring a newer calibration schema; its
existing angle-dependent terms are still evaluated. Individual control photons
are recovered, independently varied, then the pi0 is rebuilt. Failed recovery
aborts; no composite-photon approximation or silent event loss is allowed.

Deterministic per-tree/entry random numbers apply the saved fitted parameter
values to full MC. These are new full-statistics noise realizations, not a replay
of sampled optimizer-cache noise. Parameters and a seed are recorded. This study
does not redetermine cuts, retrain responses or estimate their uncertainty.

Open REPORT.html after success. Files:

* nominal_gate.json: gate status and number of checked yield cells.
* contamination_comparison.csv: production topology-summed fractions per original
  4D row, numerator/denominator MC values, changes in percentage points, and
  subtraction-only signal-yield ratios.
* topology_mc_comparison.csv: all three topology MC yields and changes separately.
  These are diagnostics; production contamination itself is topology-summed.
* summary.csv and plots/PERIOD.png: DATA-yield-weighted contamination versus |t|,
  aggregating per-row predicted background counts, not ratios of aggregate yields.
* nominal.csv / morphed.csv: full CSV copies with their reconstructed MC counts.
  The stored contamination columns remain inherited baseline values; use the
  explicitly recalculated contamination_comparison.csv for morphed fractions.
* nominal.log / morphed.log: actual production loader files and counter diagnostics.

Signal-yield ratios quantify subtraction only, not acceptance or cross sections.
The matched full-statistics farm gate remains essential even when native unit
checks pass. Tests with --native additionally compile the actual counter backend
and exercise response identity, phi range, missing-energy changes and rejection
of control events lacking individual-photon recovery inputs.
