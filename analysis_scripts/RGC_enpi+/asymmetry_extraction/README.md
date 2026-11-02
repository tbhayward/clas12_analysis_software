# RGC helicity-bilinear companion extraction

This directory contains two intentionally separate tools.

## 1. Seven-combination production companion

`extract_helicity_bilinear_combinations.py` consumes the production
`structure_function_ratios.csv` and detailed JSON from
`extract_structure_function_ratios.py`. It applies an exact, invertible linear
change of basis from

`(u1, u2, lu1, ul1, ul2, ll0, ll1)`

to seven normalized amplitude bilinears. It transforms the complete Minuit
statistical covariance according to

`V_H = M V_F M^T`.

This gives the same likelihood minimum and information as fitting the seven
new variables directly, while avoiding a duplicate copy of the 15,000-line
production analysis.

Run from the original `asymmetry_extraction` directory, for example:

```bash
python /path/to/extract_helicity_bilinear_combinations.py \
  --input-csv output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv \
  --input-json output/asymmetry_extraction/nominal/json/structure_function_ratios.json \
  --output-dir output/helicity_bilinear_extraction
```

The output contains a table, one exact statistical covariance and correlation
matrix per bin, a machine-readable definition file, and four seven-panel
figures.

The production CSV stores assigned systematic magnitudes rather than complete
signed systematic covariance matrices. Systematic columns produced here are
therefore explicitly named `*_diagonal`. For exact systematic propagation,
transform each signed alternative-fit vector with the same matrix before
forming its assigned magnitude.

## 2. Quick publication-table plotter

`plot_published_five_helicity_bilinears.py` hard-codes Tables 25, 28, 31, and
34 of the supplied analysis note. It creates four figures containing the five
bilinears available from the intended published observables alone. It does not
use the unpolarized cos(phi) or cos(2phi) nuisance results.

```bash
python plot_published_five_helicity_bilinears.py
```

Because the analysis-note tables do not publish the observable covariance,
the quick script propagates statistical and point-to-point systematic errors
diagonally and labels them accordingly. The production companion should be
used for final numerical work.

## Amplitude-index convention

Every label uses the PDF convention `T_{mu,lambda}^{0 sigma}`:

- `mu`: virtual-photon helicity (`-1`, `0`, or `+1`);
- `lambda`: initial-proton helicity (`-1/2` or `+1/2`);
- upper `0`: spin-zero pion helicity;
- `sigma`: recoil-neutron helicity, summed over `-1/2` and `+1/2`.

The plotted quantities are normalized bilinears, not individual complex
helicity amplitudes.
