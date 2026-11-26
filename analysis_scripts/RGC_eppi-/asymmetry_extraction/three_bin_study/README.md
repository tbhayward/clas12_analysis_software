# Three-bin polarized extraction and LL cosφ constraint study

Run after the relaxed campaign's nominal selection, dilution and 24-bin fits.
From the farm repository root, inside tmux:

```bash
python3 analysis_scripts/RGC_eppi-/asymmetry_extraction/three_bin_study/run_three_bin_study.py --source /volatile/clas12/thayward/eppi-/loose/full_20261124/analysis --workers 12
```

The original fine LU extraction remains unchanged. This study merges adjacent
pairs of selected-event bins into three −t′ bins per xB: [0.05,0.45),
[0.45,0.85), [0.85,1.25) GeV². It fits all five polarized ratios in 12 bins.
It recomputes the shared-period ND3 dilution from summed target counts with
2,000 Poisson replicas per bin. Charges remain unchanged when summing bins.
Each saved target-count vector is checked against the selected ROOT files.

Three combined RGB/RGC fits are run per bin:

- `free`: both UU ratios fixed zero; all five polarized ratios free.
- `step1_2`: same, but LL cosφ fixed to the supplied BM stage-2 ratio.
- `step1_2_3`: same, but LL cosφ fixed to the supplied full-BM ratio.

Both UL ratios remain free. No UL constraint is applied. LU in these merged fits
is a cross-check; the nominal LU result remains at six −t′ bins per xB.

The two bundled theory files are the supplied π− photon-axis predictions. At the
user's request they are applied directly to the lab-axis extraction without a
frame conversion. Values are linearly interpolated at the combined six-period
mean −t′ in each bin, using the four supplied xB/Q² settings. No extrapolation,
event folding, theory uncertainty or fitting of theory parameters is performed.
The actual fixed values and source checksums are saved. This is a conditional
sensitivity study and does not replace the free baseline automatically.

A fixed physical LL cosφ ratio implies beta_LLcos = dilution × fixed ratio in
the likelihood. The free-amplitude response to dilution is obtained by implicit
differentiation of the score equations. That response and the direct division
by dilution are included in first-order covariance propagation. The fixed ratio
itself has zero conditional error. UU and LU can receive an indirect dilution
contribution through the LLcos constraint. Three starts, optimizer/covariance
checks, positivity diagnostics and exact convex-hull state checks are retained.

Full output: `SOURCE/three_bin_LLcos_study/`. Includes 12 dilution fits, 36
polarized fits, comparison CSV, conditional error ratios, sensitivity shifts,
correlations, 10 PNG figures and the original fine LU reference. Shift / baseline
error is not a significance because fits share the same events.

A compact review directory and archive are written automatically beside the
analysis scripts as `review_three_bin_LLcos` and `review_three_bin_LLcos.tar.gz`.
ROOT files, NPY caches and replica ensembles are omitted from that review copy.
Override with `--review-output`. The source data stay on volatile.
Checkpoints resume only when their input and code signatures match.

`--prepare-only` checks saved counts/metadata/mean-kinematics consistency and
writes the exact rebinning and theory inputs without reading selected ROOTs or
running fits. `--df-replicas` overrides dilution trial count for explicit tests.
A production run should retain the default 2,000.

Validation:

```bash
python3 -m unittest discover -s analysis_scripts/RGC_eppi-/asymmetry_extraction/three_bin_study -p test_study.py -v
```
