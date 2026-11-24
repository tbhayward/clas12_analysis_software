# Selected kinematic distributions

Matches the eight plot products in `RGC_enpi+/kinematic_distributions/plot_kinematic_distributions_review.py` and the pi+ analysis note Section 3.7 (Figures 31–37). The default output contains separate combined RGB D2 and combined RGC ND3 sets, plus one normalized RGB/RGC comparison. No control targets are mixed into these plots.

Run from the repository root **after nominal channel selection finishes**:

```bash
python3 analysis_scripts/RGC_eppi-/kinematic_distributions/plot_kinematic_distributions.py --workers 6
```

All outputs are saved to `kinematic_distributions/output/` beside the script. The default input registry is `channel_selection/output/selected_inputs.json`. Cut-set identifiers, dataset identity and analysis-bin annotations are verified. Selection is not reapplied. Uses pion `phi2` (radians converted to 0–360 degrees) and `p2_p`. All angles stored in the original files remain unchanged.

Eight PNG pages per sample: 1D xB/Q2/−t′/phi/W/y, the three core correlations, three phi correlations, phi projections by xB, Q2 distributions in all 24 bins, W distributions in all 24 bins, Q2–W correlations by xB, and pion momentum in all 24 bins. Like the pi+ code, the Q2/W pages have six −t′ panels with four xB curves, and the pion-momentum page has four xB panels with six −t′ curves. Individual curves are normalized to the entries in their displayed range. These are raw selected combination shapes, not luminosity-normalized yields or acceptance corrections. 2D colors show unweighted counts on logarithmic scales.

Add `--per-period` for all six original period samples as well:

```bash
python3 analysis_scripts/RGC_eppi-/kinematic_distributions/plot_kinematic_distributions.py --workers 6 --per-period
```

The script scans only eight necessary branches, once per selected file, using six independent workers and bounded chunks. It caches histograms for replotting:

```bash
python3 analysis_scripts/RGC_eppi-/kinematic_distributions/plot_kinematic_distributions.py --plots-only
```

`--inputs FILE` and `--output PATH` override the registry and output directory. Q2 is plotted over 1–12 GeV² and W over 2–4.5 GeV, with common edges for RGB and RGC. Out-of-range counts are recorded in `kinematic_summary.json`. The portable pion-momentum CSV uses the actual pi+ code's 0.50–10.50 GeV range in 0.25-GeV bins. `bin_kinematic_means.csv` supplies per-period and combined-bin means and standard deviations. Figures are PNG only. `index.html` links the output.

`09_proton_momentum_FD_CD.png` adds the detected-proton momentum distribution with all/FD/CD curves and separate FD and CD momentum-versus-polar-angle maps. FD=detector1 1 and CD=detector1 2. Track p1_theta is converted to degrees for these plots. The three momentum curves share the all-proton normalization so the FD/CD proportions are visible. Momentum covers 0–5 GeV and polar angle 0–180 degrees, with tails recorded in the summary. The default therefore produces 19 PNG pages (nine for RGB, nine for RGC and one comparison). With `--per-period` there are 73 PNG pages.

Analysis-bin overlays use an explicit six-color palette, with fixed colors for each xB or −t′ interval across RGB and RGC. After channel-selection version 2 completes, rerun normally to read the new selected samples. --plots-only rejects a histogram cache whose cut-set identifier differs from the current selected-input registry.
