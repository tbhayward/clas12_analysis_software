#!/usr/bin/env python3
"""
plot_ImH_fit_results.py

Load fitted GPD-H parameters from a text file and compare the original (VGG) vs. fitted
Im H(ξ,t) ansatz. Produces a 2×2 panel for four ξ values, each showing Im H vs. −t
for three t‐values (0.1, 0.4, 0.7 GeV²), with solid lines for the original parameters
and dashed lines for the fit. Saves to output/plots/ImH_dependence_<TIMESTAMP>.pdf.

Usage:
    python plot_ImH_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""

import sys
import os
import re
import numpy as np
import matplotlib.pyplot as plt

def parse_fit_results(fname):
    """
    Parse fit_results_<TIMESTAMP>.txt, returning:
      - params: dict of fitted parameters
      - chi2_ndf: float
      - timestamp: str
    """
    # extract TIMESTAMP from filename
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt', fname)
    if not m:
        raise ValueError("Filename must contain fit_results_<YYYYMMDD_HHMMSS>.txt")
    timestamp = m.group(1)

    params = {}
    chi2 = None
    ndf  = None

    with open(fname) as f:
        for line in f:
            line = line.strip()
            # look for lines like " alpha0     = 0.398334 ± 0.0697363"
            mparam = re.match(r'(\w+)\s*=\s*([-\d\.eE]+)\s*±\s*([-\d\.eE]+)', line)
            if mparam:
                key, val, err = mparam.groups()
                params[key] = float(val)
            # look for "χ²/ndof    = 2480.42/2973 = 0.834315"
            mchi = re.match(r'χ²/ndof\s*=\s*([-\d\.eE]+)/(\d+)\s*=\s*([-\d\.eE]+)', line)
            if mchi:
                chi2, ndf, perndf = mchi.groups()
                chi2 = float(chi2)
                ndf  = int(ndf)
    if chi2 is None:
        raise ValueError("Could not parse χ²/ndof from file")
    return params, chi2/ndf, timestamp

def ImH(xi, t, p):
    """
    Im H(ξ,t) ansatz:
      r = 0.9
      α = α0 + α1 * t
      n = n_val
      b = b_val
      Mm2 = Mm2_val
      P = P_val
      pref = π·5/9·n·r / (1+ξ)
      ImH = renormImag * pref * (2ξ/(1+ξ))^(−α)
                       * ((1−ξ)/(1+ξ))^b
                       * [1 − ((1−ξ)/(1+ξ))·t/Mm2]^(−P)
                       * 2
    """
    r   = 0.9
    α   = p['alpha0'] + p['alpha1'] * t
    n   = p['n_val']
    b   = p['b_val']
    Mm2 = p['Mm2_val']
    P   = p['P_val']
    pref = np.pi * 5.0/9.0 * n * r / (1 + xi)
    xfac = (2*xi/(1+xi))**(-α)
    yfac = ((1 - xi)/(1+xi))**( b )
    tfac = (1 - ((1 - xi)/(1+xi))*t/Mm2 )**(-P)
    return p['renormImag'] * pref * xfac * yfac * tfac * 2.0

def main():
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)

    fitfile = sys.argv[1]
    params_fit, chi2ndf, timestamp = parse_fit_results(fitfile)

    # original VGG parameters
    params_orig = {
        'renormImag': 1.0,
        'alpha0':     0.43,
        'alpha1':     0.85,
        'n_val':      1.35,
        'b_val':      0.4,
        'Mm2_val':    0.64,
        'P_val':      1.0
    }

    # four ξ values (e.g. ξ ≃ x_B/(2−x_B) for x_B=0.1,0.2,0.3,0.4 → ξ≈0.053,0.111,0.176,0.242)
    xi_vals = [0.05, 0.10, 0.20, 0.30]

    # three -t values (GeV²)
    t_vals = [0.1, 0.4, 0.7]
    colors = ['tab:blue', 'tab:orange', 'tab:green']

    # ensure output directory
    outdir = f"output/plots"
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"ImH_dependence_{timestamp}.pdf")

    # prepare a fine t-array for smooth curves
    t_plot = np.linspace(0.05, 1.0, 200)

    # Set up LaTeX style
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(8, 8))
    fig, axes = plt.subplots(2, 2, sharex=True, sharey=True)

    for ax, xi in zip(axes.flat, xi_vals):
        for t, c in zip(t_vals, colors):
            # original (solid)
            ax.plot(t_plot, ImH(xi, t_plot, params_orig),
                    linestyle='-', color=c,
                    label=fr'orig, $-\,t={t:.1f}$')
            # fitted (dashed)
            ax.plot(t_plot, ImH(xi, t_plot, params_fit),
                    linestyle='--', color=c,
                    label=fr'fit,  $-\,t={t:.1f}$')
        ax.set_title(rf"$\xi = {xi:.2f}$")
        ax.grid(True)
        ax.set_xlim(t_plot.min(), t_plot.max())

    # labels & legend
    fig.text(0.5, 0.04, r"$-\,t\,$[GeV$^2$]", ha='center', va='center')
    fig.text(0.06, 0.5, r"$\mathrm{Im}\,H(\xi,t)$", ha='center', va='center', rotation='vertical')

    # one shared legend
    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize='small', ncol=1,
               title=r'Original vs. Fit')

    fig.tight_layout(rect=[0, 0.03, 1, 0.95])
    fig.suptitle(r"\textbf{Comparison of Original vs. Fitted ImH}", y=0.98)

    plt.savefig(outpath)
    print(f"[+] Saved ImH dependence plot to {outpath}")

if __name__ == "__main__":
    main()