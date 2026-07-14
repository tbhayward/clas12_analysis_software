# Final systematic-column update

## CSV schema

Removed the redundant `BSA, sigma, ...` columns. The `BSA, counts, ...` cells already store `(value, stat, sys)` tuples.

The BSA central-value tuple columns now appear first, followed by the BSA systematic columns for `10.6 GeV`, `Fa18`, `Sp18` and `Sp19 Inb`:

- `BSA, counts, <group>, combination sys`
- `BSA, counts, <group>, beam polarization sys`
- `BSA, counts, <group>, total scale sys`

The beam-polarization component is an absolute BSA uncertainty equal to `0.04*abs(A_LU)`. The total BSA scale uncertainty is the quadrature sum of that component and the absolute BSA combination systematic.

For each combined cross-section target, the schema now includes:

- `..., combination sys`
- `..., target thickness and charge sys`
- `..., total scale sys`

The target-thickness/absolute-charge component is stored as the constant fractional uncertainty `0.0476`. The total cross-section scale systematic is the row-dependent fractional quadrature sum of `combination sys` and `0.0476`.

## Production and plotting changes

- `initialize_pass2_csv.cpp`: reorganized the schema and removed redundant BSA sigma columns.
- `combination_systematics.cpp`: fills the new component and total-scale columns after the final theta-p-dependent combination values are assigned.
- `main_systematics.cpp`: creates and validates the complete new schema.
- `external_scripts/plot_pass2_vs_pass1_model_comparison.py`: uses the CSV point-to-point total for vertical pass-2 errors and the CSV total scale systematic for scale boxes and ratios. BSA panels now show pass-2 total-scale boxes.
- `external_scripts/integrated_cross_section_study.py`: uses the CSV point-to-point and total-scale systematics instead of the fixed 10% estimate.
- `external_scripts/clas6_cross_check.py`: uses the CSV point-to-point and total-scale systematics instead of the fixed pass-2 estimate.
- `external_scripts/hall_a_cross_check.py`: uses the CSV point-to-point and total-scale systematics instead of the fixed 15% estimate.

The pass-1 prescriptions in the comparison/cross-check scripts are unchanged.
