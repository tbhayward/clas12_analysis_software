# Photon ordering investigation

Source reviewed directly with `git show dvcs_processing_branch:PATH`; no branch switch. The inspected branch commit is `5a8b14f1a8a57f638079b99a0ea409fbd35a3c18`. This is source evidence consistent with the observed trees, not a claim that the branch tip identifies the exact historical processing revision for every file.

## Confirmed observations

The existing `cd_proton_angle_response_001/photon_order_audit.csv` shows photon 1 harder in approximately 99% of both-FD MC control events but only 0.01–0.02% of DATA control events, in all five periods. These are conditional ordering fractions, not yield or cross-section ratios.

The processing-branch DATA script `processing_scripts/processing_exclusive_pi0.groovy`, lines 269–290, loops over both ordered photon permutations and overwrites photon-specific detector tags, electron opening angles and azimuths together. For exactly two photons, the final accepted permutation is `(1,0)`.

The rebuilt MC path in `processing_scripts/processing_mc_rec_exclusive_pi0.groovy`, `choosePhotonPair` (lines 78–99), loops `j=i+1`, retains one mass-selected pair, and passes that analyzer to `pairedKinematics`. The latter writes detector, opening-angle and azimuth fields from the same pair and photon role. Thus a two-photon pair is labelled `(0,1)`. Lower index need not mathematically imply higher energy, but the measured ordering fractions show that it usually does in these samples.

The DATA script sets the control `detector2` label from `detector_gamma1` (line 311). Rebuilt MC control readers explicitly use `detector_gamma1`. Consequently, mixed FD/FT pairs can be grouped differently even when both photon coordinate/tag associations are correct. Both-FD pairs retain the same detector grouping under permutation.

## Checks of the response implementation

`recover_photons` receives each photon's azimuth and electron opening angle in matching order and stores the solved momentum/direction in that role. Its detector tag remains in the corresponding slot. `vary_control` selects each photon's response using its own tag and rebuilds the summed four-vector.

`test_photon_order_audit.py` creates a known physical photon pair, recovers it in both orders, and checks that recovered coordinates exchange consistently. Swapping both coordinates and detector tags preserves the morphed sum for distinct FD/FT biases; swapping tags alone changes it and is detected. Existing response tests also cover independent photon fluctuations. The observed ordering reversal therefore is not, by itself, evidence for a wrong detector response or a broken composite four-vector.

A separate limitation exists for multi-photon candidates: the DATA script retains late pair metadata but analyzes composite PID 111 occurrence zero; the rebuilt MC sidecar chooses the pair closest to the nominal pi0 mass while also analyzing composite occurrence zero. Pair/composite agreement is not guaranteed by those selection rules. Current recovery checks constrain momentum/mass consistency and flag unresolved pairs, but these alone do not identify every possible source-pair ambiguity. This investigation does not change processing or drop unresolved events.

## Next quantitative audit

```sh
git pull && python3 test_photon_order_audit.py && python3 study_photon_order_audit.py --frozen-study "$PWD/test_runs/frozen_response_kinematics_001" --output "$PWD/test_runs/photon_order_audit_001"
```

Uses existing checked 28-column frozen-study arrays on scratch. No ROOT reads, fits, production writes or regenerated caches. All periods are analyzed for DATA, nominal MC and frozen-response MC. Transfer the small output directory and open `REPORT.html`.

- `ordering.csv`: both photon tags, role-energy order, counts, ties and unresolved coverage; full globals and existing common DATA windows.
- `topology_regrouping.csv`: saved control topology versus the topology obtained using the harder photon's own detector tag consistently. FD–FT (FD proton, FT harder photon) is explicit even though it is absent from the current three-group scheme.
- `plots/PERIOD.png`: ordering fractions and row-normalized regrouping matrices for the three samples.

These are descriptive regroupings of already retained events. Global cuts and windows are not recalculated using the harder-photon grouping; absent events cannot be recovered. Neither matrix fractions nor DATA/MC event counts measure acceptance changes. A common grouping definition must be chosen explicitly, and its selection and transfer must then be reevaluated, before using it in fits or contamination calculations. Harder-photon grouping is a diagnostic candidate, not an automatic recommendation to change the control selection or its 2/0.4 GeV thresholds.
