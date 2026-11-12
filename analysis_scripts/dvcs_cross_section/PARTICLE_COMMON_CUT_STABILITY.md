# Particle-response common-cut stability

This diagnostic freezes the existing responses and applies **identical windows** to
DATA, DVCS MC and pi0 leakage MC. It does not derive new windows, refit particle
parameters, or apply the production cut finder's additional histogram morph.
Production output and previously completed studies are unchanged.

Run from the dvcs_cross_section directory (all five source periods, all three
topologies and six |t| bins from 0.11 to 1 GeV² by default):

```sh
python3 test_particle_common_cut_stability.py && python3 study_particle_common_cut_stability.py --followup-study "$PWD/test_runs/particle_response_add_pt_followup_001" --closure-study "$PWD/test_runs/particle_response_global_closure_001" --output "$PWD/test_runs/particle_response_common_cuts_001"
```

Optional `--periods fa18_inb` and `--topologies CD_FT FD_FD` restrict the study.
No compilation or ROOT read is performed, and no new large scratch cache is made.
The pre-cut scratch exports retained from the follow-up are required; a missing
cache fails explicitly before creating output. Their hashes are checked against
the completed closure study, streaming in 8 MiB blocks. Processing is sequential
and memory bounded; BLAS thread counts are one. Events are grouped by topology and
|t| before plotting/scanning. Threshold counts use one sort and binary searches,
not one event loop per threshold. Scans do not need a worker option.

Open REPORT.html. Each period/topology/|t| cell has two plots:

* `_shapes.png`: full-global pre-exclusivity shapes, then common-window N−1
  shapes (other saved DATA windows retained). These first two rows independently
  match the **summed MC** area to DATA while preserving signal/background relative
  weights. The final row uses fixed selected-anchor normalization, so its area
  need not equal DATA. Histogram under/overflow participates in area normalization.
* `_stability.png`: vary a common pTmiss upper bound or a symmetric Mx2 window
  about the DATA-window center, with other cuts fixed. The rows use saved DATA
  windows versus the add-pT candidate's refitted DATA windows. Each window family
  is applied unchanged to every response and every one-photon sample. The dotted
  line is the reference threshold; if that variable has no saved cut, defaults
  are 0.20 GeV and 0.02 GeV². The reference is normalized to one by definition.

The three responses are unmodified MC, the original-training fitted response,
and the add-pT fitted response. All include consistently modified leakage and
two-photon MC. The two-photon normalization is taken once per period using the
**same saved DATA control windows on control DATA and MC**. It fixes the leakage
scale using the saved exposure ratio and response's transfer multiplier. The
DVCS scale then matches the selected one-photon yield once per topology. No
observable-dependent mixture fit is performed. Invalid control/residual values
are flagged instead of clipped.

For each threshold c the diagnostic is R(c) = (D(c) − alpha B(c)) / S(c), with
MC sampling probabilities accounted for. It plots R(c)/R(reference) within the
same |t| bin and response. `cut_stability.csv` contains the underlying unit-weight
counts, subtraction, absolute signal scale, relative scale and uncertainty.
`normalizations.csv` records the selected anchors and flags.

Statistical errors use the actual overlap of the threshold and reference samples
and propagate the same control-count normalization into both. They are conditional
on fitted particle responses, cut windows, exposure/transfer multiplier and the
Poisson counting approximation. They exclude model and fit uncertainties. A
reference point has zero conditional error because its ratio is identically one;
this does not mean its absolute normalization is known exactly.

A flat curve supports selection stability, not absolute acceptance validation.
The exports contain |t| but not xB, Q² or track polar angles. Differing kinematic
populations, correlated cuts or background transfer can therefore also cause
drift. These are diagnostic unit-weight counts, not final corrected cross sections.
A data-fitted 95% marginal window is not proof of 95% joint signal efficiency.
Do not use a tighter cut or freely adjusted background per plot to conceal a drift.
