# Total systematic and median projection update

## Point-to-point total

`pass1_systematics_import.cpp` imports only the four legacy pass-1 components:

- pi0 subtraction
- acceptance
- Frad
- Fbin

It no longer imports or reconstructs the pass-1 four-component total as the final pass-2 total. After the four imported values are assigned, the pass-2 total is recalculated bin by bin as the quadrature sum of all six available components, including exclusivity and fiducial cuts.

If any required component is unavailable or invalid, the total is left blank rather than silently treating it as zero.

## Projection plots

Only `systematic_projection_plots.cpp` changes from mean to median aggregation. The underlying CSV systematic values and the bin-by-bin quadrature total are unchanged.

For each projected kinematic group, the plots now show:

- median absolute systematic value
- median of the per-bin relative values, `100*s/abs(cross section)`

The axes use log scale with lower limits of `1e-5 nb/GeV^4` and `1e-2 %`. Extra intermediate `2x10^n` labels are disabled.
