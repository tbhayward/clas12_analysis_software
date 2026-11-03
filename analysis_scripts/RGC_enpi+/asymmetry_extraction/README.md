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
matrix per bin, a machine-readable definition file, four single-`x_B`
seven-panel figures, and one all-`x_B` overlay figure.

The production CSV stores assigned systematic magnitudes rather than complete
signed systematic covariance matrices. Systematic columns produced here are
therefore explicitly named `*_diagonal`. For exact systematic propagation,
transform each signed alternative-fit vector with the same matrix before
forming its assigned magnitude.

## 2. Quick publication-table plotter

`plot_published_five_helicity_bilinears.py` hard-codes Tables 25, 28, 31, and
34 of the supplied analysis note. It creates four single-`x_B` figures and one
all-`x_B` overlay containing the five combinations available from the intended
published observables alone. It does not use the unpolarized cos(phi) or
cos(2phi) nuisance results.

```bash
python plot_published_five_helicity_bilinears.py
```

Because the analysis-note tables do not publish the observable covariance,
the quick script propagates statistical and point-to-point systematic errors
diagonally and labels them accordingly. The production companion should be
used for final numerical work.

Both plotters place the systematic rectangle around the central value, put the
full amplitude expression on each y axis, and start every `-t'` axis at zero.
The y ranges are calculated once from all `x_B` bins and then reused in every
canvas. Physically paired panels also share a range (the two initial-proton
helicities for the same interference), while the strength-difference panel and
unlike interference types retain independent ranges.

## Amplitude-index convention

Every label uses the PDF convention `T_{mu,lambda}^{0 sigma}`:

- `mu`: virtual-photon helicity (`-1`, `0`, or `+1`);
- `lambda`: initial-proton helicity (`-1/2` or `+1/2`);
- upper `0`: spin-zero pion helicity;
- `sigma`: recoil-neutron helicity, summed over `-1/2` and `+1/2`.

The plotted quantities are normalized bilinears, not individual complex
helicity amplitudes.

## 3. Restricted no-nucleon-flip amplitude test

`fit_no_nucleon_flip_amplitudes.py` asks a different question from the two
plotters above: can three parity-related nonflip complex amplitudes reproduce
all seven ratios in each kinematic bin?

The model keeps

```text
L = T_{0,+}^{0,+}
P = T_{+1,+}^{0,+}
M = T_{-1,+}^{0,+}
```

and sets amplitudes in which the recoil neutron has the opposite helicity to
the initial proton to zero. Parity supplies the corresponding amplitudes for
initial-proton helicity `-`. The common phase and common scale are unobservable
in ratios, so the three complex amplitudes contain four fit parameters. The
seven measured ratios therefore normally provide three goodness-of-fit degrees
of freedom in each bin.

The denominator used by the model is

```text
D = 0.5 (|P|^2 + |M|^2) + epsilon |L|^2,
```

which is `F_UU/C` in the restricted model. The common kinematic factor `C`
cancels. By default the script uses the full seven-parameter statistical
covariance from the production JSON. The optional diagonal-systematics mode is
only an approximation because the production table does not contain a signed
systematic covariance matrix.

Run it beside the production extraction with

```bash
python /path/to/fit_no_nucleon_flip_amplitudes.py \
  --input-csv output/asymmetry_extraction/nominal/tables/structure_function_ratios.csv \
  --input-json output/asymmetry_extraction/nominal/json/structure_function_ratios.json \
  --output-dir output/no_nucleon_flip_amplitude_test
```

It writes a per-bin fit summary, a long table of measured and predicted ratios,
four measured-versus-model canvases, a goodness-of-fit figure, a residual-pull
heatmap and a representative normalized amplitude solution. The last plot uses
an arbitrary common phase and scale convention and must not be read as an
absolute amplitude measurement.

A satisfactory fit means that helicity-flip amplitudes are not required by
these seven ratios at the present precision. It does not demonstrate that the
flip amplitudes do not exist. A poor fit rejects this restricted model, subject
to the same frame, convention and covariance checks that apply to the inputs.
