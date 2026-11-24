# Exclusive pi− structure-function ratios

The current preliminary default fixes both UU ratios to zero and fits the five spin-dependent ratios. Its output is `output_fixed_UU/`. The earlier seven-free-ratio mode is retained with `--float-uu` and defaults to `output_free_UU/`.

Run from the repository root:

`python3 analysis_scripts/RGC_eppi-/asymmetry_extraction/extract_structure_function_ratios.py --workers 6`

Preflight without fitting: append `--preflight`. For a short end-to-end check, append `--bins 1`. Check progress with the same command and `--status`. A repeat launch resumes checkpoints matching the selected-file identities, metadata, dilution inputs and code. Outputs default to `output/` beside the scripts. No new channel selection or original campaign scan is required.

Inputs can be changed with `--selected-inputs`, `--run-info`, and `--dilution`. Dependencies are numpy, scipy, iminuit, uproot, matplotlib (threadpoolctl is optional), plus the existing channel-selection reader. Each worker uses one numerical-library thread.

## Results

- `combined`: all six periods, with seven common ratios.
- Six separate period fits. RGB fits UU cosφ, UU cos2φ and LU sinφ. RGC fits all seven ratios, allowing comparisons of the same first three and of the four target-dependent ratios.
- `ratios.csv`, `fits/*.json`: central values, symmetric statistical errors, full covariance, correlations, positivity and optimizer diagnostics.
- `figures/nominal_*.png`: seven combined-result canvases.
- `figures/period_UU_LU_*.png`: three six-period comparison canvases.
- `figures/period_UL_LL_*.png`: four RGC-period comparison canvases.
- `input_audit.json`, `fit_summary.json`: run coverage, charge conventions, candidate/unique-event counts, bin populations and fit status.
- `period_likelihood_ratio_checks.json`: approximate omnibus comparisons of 30 independent-period amplitudes with seven shared amplitudes. The groups overlap through shared dilution information and the cosine ratios can be weakly identified, so this diagnostic is not a substitute for closure checks.

## Physics and normalization

The pion production angle is `phi2`, in radians. The model uses DepB/A, DepC/A, DepV/A and DepW/A per event. Stored beam polarizations are checked against the campaign's run-dependent table. The signed elastic target polarization is read from the same run-information data used in processing and checked against the selected ROOT values. RGB has zero target polarization. No cos(theta_gamma) projection or imported transverse-target terms enter the fit.

The conditional probability is Q(r,h)G(x;r,h) / sum_(r',h' in same period) Q(r',h')G(x;r',h'). Charge moments compress the normalization exactly, including run-dependent Pb and signed Pt. Different periods retain different normalization and acceptance groups. Every original pre-exclusivity source run contributes its exposure, including runs with no final exclusive candidate. Missing metadata and inconsistent polarization stop the calculation. Positive states only enter the normalization. The cross-section factor must remain positive for every allowed spin state at each selected event's kinematics.

RGC uses measured helicity charges in the six-column records: run,total,Qplus,Qminus,Pt,Pt_uncertainty. RGC stored helicity retains the raw sign, as does the processing campaign. New RGB records have seven columns: run,total,Qplus,Qminus,Qunassigned,Pt,Pt_uncertainty. Their unavailable measured helicity columns remain zero. With the user's explicit CEBAF assumption, the likelihood assigns Qplus=Qminus=total/2 for RGB. Consequently the RGB run-dependent helicity sign reversal does not change the exposure normalization. This is an assumption about exposure, not measured helicity-resolved charge. No empirical balancing using candidate counts is performed. Target polarization and uncertainty are zero for RGB.

## Dilution uncertainty

The nominal shared 24-bin dilution and symmetric replica standard deviations come from `../dilution_factor/output/joint_period_method1.json`. A common factor applies to all three RGC periods. The likelihood first fits the products beta_target=Df*ratio_target. This avoids a numerical scale degeneracy. The physics ratios and their covariance are then transformed at the central Df. A common Gaussian dilution nuisance adds sigma_Df²*(d ratios/d Df)*(d ratios/d Df)^T. This is exactly the local-Hessian covariance of the Gaussian-constrained nuisance formulation under this parameterization, without extra amplitude constraints. The target terms receive correlated uncertainty. UU and LU do not receive a direct dilution correction. Individual-period comparisons share the same dilution uncertainty and must not be treated as independent for subsequent significance calculations.

## Identifiability and limitations

The spin-independent UU terms cancel when all spin modulations are zero. RGB conditional helicity dependence constrains their effect on the beam asymmetry denominator, not an acceptance-free helicity-summed cross section. Three initial values are tried, Minuit validity and accurate covariance are required, and positivity and covariance conditioning are recorded. Invalid fits remain flagged and are omitted from nominal plotted points. A numerically valid result is not sufficient evidence that a weak UU ratio is physically well constrained. Inspect its correlations and multistart minima before using it. No harmonics are fixed to external pion-channel values.

The initial covariance treats retained proton–pion combinations as independent observations, as in the earlier pipeline. The audit exposes repeated run/event IDs. Event-level clustering, detector-state acceptance changes within a period, tensor polarization, background-dependent differences between RGB and RGC, and polarization scale uncertainties require separate studies. Dilution uncertainty is propagated at first order, consistently with symmetric preliminary reporting.

Tests: `python3 -m unittest discover -s analysis_scripts/RGC_eppi-/asymmetry_extraction/tests -v`

The first full local extraction completed all 168 fit tasks. Eighteen nominal bins passed the Minuit covariance requirements. Nominal bins 3, 6, 11, 12, 16 and 23 reached the cross-section positivity boundary and remain invalid rather than being assigned misleading symmetric uncertainties. Individual-period invalid fits are also retained in JSON and shown in the fit-status canvas. These are diagnostic results, not a completed physics validation.

## Production plots and methodical audit

`python3 analysis_scripts/RGC_eppi-/asymmetry_extraction/plot_production_results.py --selected-inputs analysis_scripts/RGC_eppi-/channel_selection/output/selected_inputs.json`

Future extraction runs create this production set automatically. It follows the pi+ `plot_final_results.py` seven-panel layout, with statistical error bars only. All UU and LL panels use [-1,1], all LU and UL panels [-0.4,0.4], and all horizontal axes [0,1.30] GeV². Mean selected-event tprime values define the plotted positions. Four separate xB canvases and an overlay are accompanied by four six-period UU/LU comparison canvases and four RGC UL/LL comparison canvases. No systematic bands are drawn. Edge triangles disclose intervals extending outside the fixed axes. Omitted fits, clipping and plot positions are recorded as CSV files. Old figures are retained, with the new canonical set under `figures/production/`.

`python3 analysis_scripts/RGC_eppi-/asymmetry_extraction/audit_fits.py`

The audit verifies input checksums and reconstructs all saved likelihood values. It records analytic scores and Hessian eigenvalues, active positivity constraints, the worst state/event and its harmonic decomposition, and the distance to positivity in the local Gaussian covariance. This distance is a diagnostic, not a coverage statement. Constrained-boundary stationarity is checked with nonnegative constraint multipliers. Both-UU-zero, one-UU-zero and leave-RGB-Fa19-out fits are diagnostic only. They never modify the nominal fits or the production ratios. They test parameter-identification and sample leverage, not hardware helicity calibration. See `audit/bin_audit.csv`, `audit/bin_audit.json`, `audit/audit_report.txt`, `audit/fixed_UU_zero_diagnostic.json`, and `audit/isolation_diagnostics.json`.

## RGB-only isolation

`python3 analysis_scripts/RGC_eppi-/asymmetry_extraction/rgb_only_check.py --workers 6`

This writes `output/rgb_only/`, fitting only LU and the two UU ratios to the three combined RGB periods and to each RGB period separately. All target-polarization terms are absent and no dilution uncertainty enters these fits. Separate diagnostic variants test the pi+ ±1.5 coefficient limits and its observed-factor/observed-probability support checks. The latter are a comparison only: positivity of every other allowed state is not implied by checking one observed run-state probability. None of these variants changes the nominal extraction.

The initial combined RGB run has 672,750 selected combinations, with 4,555 to 111,885 per two-dimensional bin. It gives 21/24 valid covariance estimates, with bins 10, 12 and 15 at the positivity boundary. The bound/support comparison leaves the same three failures. Individual periods give 19/24, 18/24 and 18/24 valid fits. Convergence alone is not validation of the UU estimates.

For balanced helicity charges, the conditional RGB probability measures the beam-asymmetry numerator divided by the UU angular factor. It does not obtain acceptance-free UU information from the helicity-summed phi distribution. The UU terms disappear entirely when LU is zero. `expected_information.csv` illustrates the conditional Fisher information at the measured kinematics for injected UU values (0.15,-0.15) and two LU magnitudes, 0.03 and 0.20. These are illustrative assumed truths, not measured amplitudes or closure results.

## Fixed UU inputs and later model variations

The default command fixes UU_cos_phi=UU_cos_2phi=0 in every bin and allows LU, UL sinphi, UL sin2phi, LL and LL cosphi to float in the combined RGB+RGC fit. Individual RGB period fits have only LU free. Individual RGC periods have all five spin ratios free. The common dilution uncertainty is propagated as before.

To supply nonzero inputs, use `--uu-inputs path/to/UU_model.csv`. Copy `UU_inputs_template.csv` and fill the 24 rows. Columns are exactly `bin,UU_cos_phi,UU_cos_2phi`, in the same xB-major ordering as the analysis (six tprime bins per xB interval). Values must be structure-function ratios, not raw structure functions or coefficients already multiplied by depolarization factors. Every bin must occur exactly once, with finite numbers. The same imported ratios apply to all six periods. They are held exactly fixed, with no statistical nuisance or model uncertainty assigned. They are suitable for a separate systematic variation; the input does not itself define that systematic uncertainty.

Imported runs default to `output_fixed_UU_<input-file-stem>/`; `--output` overrides this. `UU_settings.json`, the input audit and every fit JSON retain the fixed values. The import contents and SHA256 enter the checkpoint signature. CSV rows identify these values as `fixed_input`. Production panels use crosses labelled as fixed inputs rather than zero-error measurements. The period-comparison likelihood-ratio diagnostic uses 13 extra free parameters in fixed-UU mode, versus 23 for the earlier fully free mode.

A fresh local fit with zero UU values passed all 168 convergence/covariance checks. This resolves the earlier optimizer failures but does not replace the remaining physics and systematic validation.
