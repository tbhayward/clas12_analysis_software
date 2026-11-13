# Independent control response/weighting validation

```sh
python3 test_pi0_control_validation.py && python3 study_pi0_control_validation.py --control-study "$PWD/test_runs/pi0_control_transfer_001" --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/pi0_control_validation_001"
```

No DATA ROOT reads. Reuses verified DATA/add-pT extended exports on scratch.
Nominal MC alone is reevaluated from the verified two-photon cache. The worker
compiles once against saved production globals; at most 10 workers run across
periods. Existing response fits, trees, reports and production results are untouched.

The stable exported event-half bit separates training and test samples. It is
unchanged between nominal and corrected MC. No test event enters training scales
or weights. Each period is fitted separately. Both DATA and MC receive the same
saved DATA control windows. The test compares:

1. Unmodified control MC, scaled by training DATA/MC counts.
2. Fixed add-pT control MC, independently scaled by training DATA/MC counts.
3. Fixed add-pT MC with training-only coarse |t|×Q² weights, shared across all
   control topologies, and its training-only weighted total scale.

The fixed grid is |t|=[.11,.25,.4,.6,1], Q²=[1,2.5,4,12]. Cells need at least 50
DATA and MC training events by default (--min-train-count). Otherwise use a
training |t|-marginal factor if supported, finally unity. Out-of-grid coordinates
receive unity. There is no topology/photon-momentum weighting or test-driven
choice of grid. training_weights.csv records cells, counts, factors and fallback.

Open REPORT.html, heldout_totals.csv, then the all_control_topologies and FD_FD
plots. These plots contain only independent test events. A ratio of one is not
imposed; the totals may differ. Comparing nominal/add-pT locates discrepancies
that predate the particle correction. Purple tests whether a simple shared
kinematic spectrum adjustment generalizes, including to axes not used to weight.

Statistical bars condition on training scales/weights and omit training-fit,
response and model uncertainty. Comparisons share DATA and MC events and are
correlated. Sparse/zero DATA points are not precise exclusions. Do not choose a
model by summing these overlapping projection bins as independent chi-square.

Improved held-out control agreement supports this empirical control description,
not a proven generator correction or missing-second-photon efficiency. Reconstructed
coordinates can absorb detector/selection effects. Control DATA may include other
backgrounds. No weights are applied to one-photon leakage or production results;
that requires validating common-coordinate transfer and exposures separately.

Processing provenance comes from saved run inputs. Reading a file on another Git
branch uses `git show BRANCH:path`, which does not change the working checkout.
The current branch is not evidence for the historical processing version.
