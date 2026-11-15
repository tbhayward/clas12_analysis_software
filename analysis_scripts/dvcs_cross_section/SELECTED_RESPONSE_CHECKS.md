# Frozen selected-response checks

Run on the farm where the completed selection-first scratch exports still exist:

```sh
python3 test_selected_response_checks.py && python3 study_selected_response_checks.py --response-run "$PWD/test_runs/selected_particle_response_001" --joint-run "$PWD/test_runs/particle_response_joint_control_001" --output "$PWD/test_runs/selected_response_checks_001"
```

No compilation, ROOT reading, response fit, or cut redetermination is performed. The
script reads memory-mapped exports and saves only small tables and plots. It runs
serially to keep memory and disk pressure low. Export checksums and the fitted-run
scoring helpers are verified. A fresh output directory is required. The starting
joint run supplies its saved sampling probabilities; original MC cache binaries
are not opened. Missing farm scratch exports are reported before output creation.

Open REPORT.html. Plots are grouped by period, with three checks:

* `mc_response_only_t_migration.png`: same initial MC windows before/after response.
  `mc_response_and_new_cuts_t_migration.png`: initial nominal versus final modified
  selection. Separate versions cover leakage and two-photon MC. Rows are original
  |t| bins and columns modified bins, each row normalized to 100%. Outside the
  analysis |t| range and selection failure are explicit. Both-failing events are
  excluded. `t_migration_summary.csv` separates bin feed-in/out from selection
  entry/exit; counts have no production current/proton weights. Nominal/candidate
  exports come from the same ordered cache, with unchanged topology and half
  checked; this is not a join of unrelated samples.
* `final_window_scores.png`: nominal and modified MC evaluated with common final
  DATA windows, omitting the variable being displayed. Actual MC half 1 is shown;
  both halves are tabulated. The nominal support and variance are held fixed.
  The fitted background nuisance is fixed for both scenarios; control/signal area
  normalizations are recalculated from the selected counts, without fitting.
  Scores are descriptive shape diagnostics, not p-values or independent DATA
  validation. Missing quantities and momentum balance are the headline variables.
  Other angle diagnostics remain in the CSV, without driving these checks.
* `*_mass_kinematics.png`: median and central 68% half-width of the reported
  two-photon mass in FD photon momentum/angle quartiles defined by DATA and shared
  by both MC scenarios. Full common final DATA selection is applied. Both photon
  roles are pooled with identical treatment. FD–FD events contribute twice;
  entries are correlated, so no independent-event precision is claimed. FT–FT has
  no FD photons and is not drawn. Low-support cells remain NaN. Reconstructed
  coordinate migrations, generator populations and the upstream mass window can
  influence these comparisons. Unresolved photon coordinates are excluded.

These checks do not modify production analysis or convert raw counts into final
cross-section corrections. Transfer only the new output directory for assessment;
keep large arrays on scratch.
