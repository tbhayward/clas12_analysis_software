# π− fixed-UU model-input study

This is a paired refit of the original selected e′pπ− sample, without changing exclusivity, charge normalization, polarizations, dilution factors, or likelihood constraints. All six RGB/RGC periods enter the combined fit. RGB target polarization remains zero.

Three scenarios fix the UU cosφ and cos2φ ratios to zero, Brandon’s `step1_2` predictions, or his `step1_2_3` predictions. The five polarized ratios float in every case. Theory values are interpolated linearly in −t′ at the measured mean of the combined sample for each fit bin. The four supplied theory xB/Q² settings are used directly. No extrapolation and no theory-input uncertainty propagation are performed.

The script refits all 24 original bins and all eight merged bins. LU comparisons use the original six −t′ bins. The principal UL/LL comparisons use the lower/upper three-bin mergers. Existing coarse dilution factors and their statistical uncertainties are reused, after checking their cut-set and metadata provenance. No dilution refit is needed because selected events and target yields are unchanged.

`likelihood.py` is the validated likelihood with an exact convex-hull reduction of the allowed spin states used for positivity checking. Charge-weighted spin normalization still uses every original state. The original zero-UU fit is rerun and compared with saved fine/coarse results in `baseline_reproduction.json`.

Outputs are written to `output/` by default. `comparison.csv` records values, errors, fixed UU settings, shifts, convergence, and positivity diagnostics. `shift_summary.csv` summarizes shifts in units of the zero-scenario statistical uncertainty. These are correlated refits, so that quantity is a sensitivity scale, not a significance. `UU_inputs/` contains the actual fixed values used. The figures show absolute results and changes from zero-UU for fine LU and coarse UL/LL.

Dependencies: Python, NumPy, SciPy, iminuit, uproot, Matplotlib. threadpoolctl is optional. Set `PYTHONPATH` to any local dependency installation if needed.

Example (single line, substitute your analysis and input locations):

```bash
python3 run_study.py --analysis-root /path/to/RGC_eppi- --selected-inputs /path/to/channel_selection/output/selected_inputs.json --coarse-dilution /path/to/output_two_tprime_bins/coarse_dilution.json --theory /path/to/pion_theory_comparisons/data/theory --baseline-fine /path/to/output_fixed_UU --baseline-coarse /path/to/output_two_tprime_bins --workers 3
```

Every fit is checkpointed with an input and code signature. Existing fits are reused only when the signature matches. Nominal extraction results are never overwritten.
