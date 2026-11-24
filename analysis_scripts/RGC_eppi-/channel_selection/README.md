# Nominal preliminary exclusive pi− channel selection

This is the shared channel-selection workflow for dilution-factor and asymmetry analyses. The earlier `diagnostics/` and `exclusivity/` studies are retained as development history. Their caches and cut files are not automatically used by this workflow.

## Run

From `~/clas12_analysis_software` on one farm node:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --workers 6
```

The command builds six calibration caches, fits five groups, freezes a versioned cut set, and writes **19 selected ROOT files**, one per original period/target. All outputs are saved under `analysis_scripts/RGC_eppi-/channel_selection/output/` beside the nominal scripts. Inputs are the existing merged files in `/volatile`, configured by `../imports/datasets.json`.

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --status
```

Audit the finished release:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --audit
```

Python dependencies: NumPy, SciPy, uproot 5 and Matplotlib. ROOT itself is not required. No PDF figures are produced.

## Calibration groups and original samples

| Fit group | Calibration sample | Windows applied to |
|---|---|---|
| RGB spring 2019 | Spring-2019 D2 | Spring-2019 D2 |
| RGB fall 2019 | Fall-2019 D2 | Fall-2019 D2 |
| RGB spring 2020 | Spring-2020 D2 | Spring-2020 D2 |
| RGC Su22 | Su22 ND3 | Su22 ND3, C, CH2, ET, He |
| RGC Fa22+Sp23 | Pooled Fa22 and Sp23 ND3 | Fa22 and Sp23 ND3, C, CH2, ET, He, and Sp23 CD2 |

Pooling determines common windows only. It does not merge original ROOT datasets, run numbers, target labels, polarization metadata, charge normalizations or acceptance groups. Fits use unweighted candidate counts integrated over helicities. Control targets receive exactly the same windows as ND3 in their period. They are not fitted independently.

All 24 pi+ analysis bins are retained: xB edges [0.10, 0.25, 0.35, 0.45, 0.60] and −t′ edges [0.05, 0.25, 0.45, 0.65, 0.85, 1.05, 1.25] GeV². Intervals include lower edges and exclude upper edges. Bin numbers are 1–24 in xB-major order. The saved signed tprime is negated.

## Nominal cuts

The base selection requires finite necessary kinematics, fiducial_status == 3, an FD pion, an FD or CD proton, Q² ≥ 1 GeV², W ≥ 2 GeV, y ≤ 0.80, and membership in an analysis bin. The processing corrections are already in the saved lab momenta. They are not applied again.

1. Reconstruct the deuteron-reference missing mass squared `(beam + stationary deuteron − electron − proton − pion)²`. Fit the spectator-proton peak with Gaussian plus a nonnegative linear background. Select mean ± two Gaussian widths. Histograms have 0.010 GeV² bins. The initial fit range 0.70–1.30 GeV² expands as needed within 0.40–1.50 GeV². The peak mean is restricted to 0.78–1.05 GeV² and width to 0.01–0.20 GeV².
2. After the nominal mass cut, fit the full spatial missing momentum `|beam − electron − proton − pion|`. This is not transverse missing momentum. Use an asymmetric gamma signal plus a higher-mode gamma tail. Both components are positive, vanish at zero and decrease at large momentum. Momentum histograms use 0.005 GeV bins and the fitted interval is 0–0.60 GeV.
3. The momentum window is the fitted signal's **central 95.45% probability interval**, given by its 2.275% and 97.725% quantiles. This replaces Gaussian mean ± two widths when the signal model becomes asymmetric. It has the same probability coverage for an ideal Gaussian. This is model coverage, not a purity estimate.
4. Plot theta_measured − theta_predicted in degrees after both cuts. The predicted pion uses a stationary neutron. The angular Gaussian is a visual guide, not an additional cut. Original tree angles retain their original units.

The momentum fit parameterization is `Ns*Gamma(p; k, scale=m/(k−1)) + Nb*Gamma(p; kb, scale=(m+d)/(kb−1))`. Signal mode m is 0.02–0.20 GeV, k is 2–15 and d is 0.04–0.50 GeV. The nominal tail shape kb=3 is compared with kb=2. Component yields are nonnegative. Fits integrate the model over histogram bins and minimize the binned Poisson deviance. The fitted tail is an empirical description and does not identify a unique physical background process.

Mass fits compare linear and nonnegative concave-down quadratic backgrounds. Both stages compare a second histogram with twice the bin width. Shifts larger than two nominal bins or 10% of the nominal window width are flagged for review. Poor Pearson chi-square/ndf > 5 is reported separately. These checks are conditional on the nominal upstream mass cut, rather than propagated full systematic variations.

## Reference assistance in weak Su22 bins

All bins must have a documented cut before a complete nominal release can be published. Independent fits still require 100 entries, convergence, a peak away from parameter bounds and a deviance improvement of at least 9 over the background-only model. This improvement is a practical resolution screen, not a calibrated significance.

The four RGB/Fa22+Sp23 groups are fitted first. Su22 is then fitted independently. If its same-bin mass or momentum peak is unresolved:

1. Retry the independent fit with the corresponding Fa22+Sp23 peak as the optimizer starting point.
2. If still unresolved, transfer the reference signal shape. For mass, shift the reference mean by the median Su22-minus-reference mean difference and multiply its width by the median Su22/reference width ratio measured in reliable Su22 bins. For momentum, scale the reference mode by the median Su22/reference mode ratio and retain its shape parameter.
3. Prefer donors with usable, unflagged fits in both groups. If fewer than six exist, use all mutually usable bins. If six donors are still unavailable, retain the unshifted reference and flag that limitation.
4. With the signal shape held fixed, fit the Su22 signal normalization and background/tail parameters. The window is explicitly reference-assisted even if the fitted Su22 signal amplitude is compatible with zero. It is not an independent Su22 peak measurement.

Reference-assisted fits save the reference group and bin, donor bins, median correction, donor scatter, original failed fit and unadjusted reference window for later sensitivity studies. Their selected events carry `exclusivity_reference_assisted = 1`. The same assisted window is applied to every Su22 target. In the supplied full data, only the **mass cuts in Su22 bins 6 and 18** need this assistance. All five groups then have all 24 bins available. Remaining model-sensitivity flags are retained as preliminary-release provenance.

A missing usable reference cannot be replaced silently. A full nominal cut set is not published if a bin remains unavailable after this procedure. Partial smoke-test fits can be saved for diagnostics but cannot be loaded as nominal selections.

## Outputs and provenance

- `nominal_cuts.json`: frozen numerical windows, all five groups, bin edges, version `eppi_minus_preliminary_v1`, code fingerprint, source calibration identities, assisted-bin provenance and a SHA-256 cut-set identifier.
- `fits/GROUP/`: mass and momentum fit results, cut and stability tables, per-member cut flows, angular guides and **12 PNG pages per group** (60 in total).
- `selected/selected_DATASET.root`: 19 separate selected `PhysicsEvents` TTrees. Every original scalar branch is preserved with its original values and type. Added branches are `analysis_bin`, `exclusivity_group`, `Mx2_deuteron`, `p_miss`, `delta_theta_deg` and `exclusivity_reference_assisted`.
- ROOT metadata objects `selection_cutset_id` and `selection_manifest`: exact cut identity, original source identity, period/target label, fit group, selected counts, original schema and full/partial status.
- `datasets/DATASET.json` and `.csv`: per-bin counts, per-run combination counts and provenance. Input/selected counts are candidate combinations, not unique events.
- `selected_inputs.json`: the downstream input registry, including relative selected-file paths, group identities, reference-assisted bins and branch conventions. Relative paths remain valid when the output folder is copied.
- `calibration_cache/`: bounded, double-precision derived candidate chunks. Calibration ROOT files are scanned once. Group fits, cut changes and reference checks reuse these chunks. The initial nominal workflow builds its own cache and does not depend on diagnostic caches.
- `index.html`: figures and dataset tables.

The cut-set checksum is stored in every selected ROOT file. The audit checks dataset coverage, branch names/types, entry counts, group labels, reference-assistance labels and that every exported candidate satisfies the frozen windows. Selected files are written through temporary files and published after validation. Matching completed files are reused on restart.

For export speed, the saved missing energy and neutron missing mass give `p_miss² = Emiss² − Mx2` and `Mx2_deuteron = Mx2 + 2*(Md−Mn)*Emiss + (Md−Mn)²`. These identities are checked against independent four-vector reconstruction for each source chunk. Angular residuals are evaluated only for retained candidates during full-branch export. Missing-energy shortcuts are not used if the branch is absent.

## Use for dilution factors and asymmetries

Use the selected ROOT files or import the same selection class when scanning original files. Do not reimplement numerical cuts in a new script. To import from a sibling analysis script, add the channel_selection directory to Python's module path, then:

```python
from selection import ChannelSelection, load_selected_inputs
registry = load_selected_inputs(channel_selection_dir / 'output' / 'selected_inputs.json')
selection = ChannelSelection(registry['cut_file'])
# Already selected files:
selected_path = registry['files']['RGC_spring2023_ND3']['path']
# Or apply exactly the same cuts to a chunk of original scalar arrays:
mask, added_branches = selection.evaluate(arrays, 'RGC_spring2023_ND3')
```

`load_selected_inputs` resolves paths relative to the registry and checks ROOT cut-set identifiers, entry counts and dataset identities. `evaluate` returns a full-length mask plus annotation arrays aligned with selected rows. It does not alter or weight any original variable.

The production azimuth remains `phi2`. The original `helicity`, `beam_pol`, `target_pol`, beam energy and depolarization factors remain in the selected files. RGB target polarization is exactly zero. Carbon, CH2, CD2, ET and He also have no target polarization. ND3 retains its original signed, run-dependent values. The exporter rejects nonfinite downstream metadata or invalid selected beam helicities, and rejects a nonzero target polarization in an unpolarized dataset.

**Do not derive luminosity or accepted charge from selected candidate counts.** Keep QA-selected run-based charge and polarization normalization independently of event-level exclusivity selections. Runs with zero selected events can still carry charge. The registry preserves original period/target separation for that purpose. This workflow does not yet implement the CD2/CH2 transfer, dilution calculation or asymmetry likelihood. It supplies their common event selection and input contract.

## Restart and options

Apply existing frozen cuts without refitting:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --workers 6 --apply-only
```

Recalibrate without exporting selected TTrees:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --workers 6 --fit-only
```

Force fresh calibration ROOT scans:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --workers 6 --rebuild-cache
```

A limited calibration smoke test must use a separate destination:

```bash
python3 analysis_scripts/RGC_eppi-/channel_selection/run_channel_selection.py --workers 6 --fit-only --max-entries 200000 --output analysis_scripts/RGC_eppi-/channel_selection/output/smoke_test
```

Use `--input-root PATH`, `--datasets FILE` or `--output PATH` to override defaults. No previous diagnostic cache or source ROOT file is deleted automatically.

## Validation

```bash
python3 -m unittest discover -s analysis_scripts/RGC_eppi-/channel_selection -p 'test_nominal_selection.py' -v
```

Tests cover gamma-signal recovery and falling tails, model probability coverage, shared control-target cuts, reference transfer and provenance, frozen-cut tamper checks, full branch/polarization preservation, zero-entry output, resuming exports and agreement of fast missing-energy identities with full four-vectors. An end-to-end six-worker synthetic run creates and audits all 19 selected files. The supplied complete candidate caches were also fitted with the five-group nominal procedure. Full selected production files still require the farm command above.
