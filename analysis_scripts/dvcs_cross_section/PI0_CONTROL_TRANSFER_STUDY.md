# Differential two-photon control check

Run from dvcs_cross_section:

```sh
python3 test_pi0_control_transfer.py && python3 study_pi0_control_transfer.py --common-study "$PWD/test_runs/particle_response_common_cuts_001" --workers 10 --scratch-root /scratch/thayward/dvcs_particle_response --output "$PWD/test_runs/pi0_control_transfer_001"
```

Uses all five periods separately, fixed add-pT response, saved global selection and
saved DATA control windows on both DATA and MC. The fixed control scale is the
one recorded in the common-cut study. It is not refitted per bin or topology.

Open REPORT.html, check anchor_check.csv first, then all_control_topologies.png
and FD_FD.png for Fa18 inb. Each point is a disjoint kinematic interval of xB, Q²,
|t|, electron momentum/angle or proton momentum/angle. It shows DATA divided by
the fixed-scale two-photon MC prediction. The total across all topologies should
reproduce the previous anchor (within counting rounding); individual bins need
not equal one. Errors are conditional on that fixed scale and omit its shared
model uncertainty. MC sampling probabilities are included.

The source analysis_fitter.java currently applies p>2 GeV for num_photon==0,
and p>0.4 GeV otherwise. dvcs_fitter.java applies p>2 GeV. These source thresholds
do not prove the versions used to produce existing ROOT trees. The diagnostic
inherits upstream photon thresholds and cannot reconstruct removed photons.
It does not treat the composite pi0 as an individual leakage photon or impose a
second-photon threshold on it. Any future differential reweighting must use
coordinates common to control and leakage, not pi0 momentum equated to photon
momentum. Control backgrounds must also be assessed before extracting weights.

One ROOT read of two-photon DATA per period, existing two-photon MC cache, no
DVCS/leakage MC reread, no parameter fits, no production updates. Extended arrays
stay on scratch; transfer only the report directory. The standalone worker is
compiled automatically against the saved global cuts. At most 10 total workers.
The previous optional 20-column export is reused without modifying its format.
A future replay can use --reuse-exports PREVIOUS_OUTPUT with the same inputs.

This is a control validation, not an automatic corrected background calculation.
No leakage weights are applied before the differential control, exposure and
selection assumptions have been validated. A matching control alone does not
validate the probability of losing a second photon, especially near threshold.
