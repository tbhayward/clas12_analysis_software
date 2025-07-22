#!/usr/bin/env python3
"""
plot_ImH_fit_results.py

Load fitted GPD-H parameters from a text file and compare the original (VGG) vs. fitted
Im H(ξ,t) ansatz. Produces a 2×2 panel for four ξ values, each showing Im H vs. −t
for three t‐values (0.1,0.4,0.7 GeV²), with solid lines for the original parameters
and dashed lines for the fit. Saves to output/plots/ImH_dependence_<TIMESTAMP>.pdf.

Usage:
    python plot_ImH_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""

import sys, os, re
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
    chi2_ndf = None

    with open(fname, encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # match lines like " alpha0     = 0.398334 ± 0.0697363"
            mparam = re.match(r'(\w+)\s*=\s*([-\d\.eE]+)\s*±\s*([-\d\.eE]+)', line)
            if mparam:
                key, val, err = mparam.groups()
                params[key] = float(val)

            # match either Unicode χ²/ndof or ASCII chi2/ndof
            mchi = re.search(
                r'(?:(?:χ²)|(?:chi2))\s*/\s*ndof\s*=\s*([-\d\.eE]+)\s*/\s*(\d+)\s*=\s*([-\d\.eE]+)',
                line, flags=re.IGNORECASE
            )
            if mchi:
                num, denom, perndf = mchi.groups()
                chi2_ndf = float(perndf)  # this is χ²/ndof already

    if chi2_ndf is None:
        raise ValueError("Could not parse χ²/ndof (or chi2/ndof) from file")

    return params, chi2_ndf, timestamp

def ImH(xi, t, p):
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

    # four ξ values (common DVCS range)
    xi_vals = [0.05, 0.10, 0.20, 0.30]
    t_vals  = [0.1, 0.4, 0.7]
    colors  = ['tab:blue', 'tab:orange', 'tab:green']

    outdir = "output/plots"
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"ImH_dependence_{timestamp}.pdf")

    t_plot = np.linspace(0.05, 1.0, 200)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig, axes = plt.subplots(2, 2, figsize=(8,8), sharex=True, sharey=True)

    for ax, xi in zip(axes.flat, xi_vals):
        for t, c in zip(t_vals, colors):
            ax.plot(t_plot, ImH(xi, t_plot, params_orig),
                    '-', color=c,
                    label=fr'orig, $-\,t={t:.1f}\,\mathrm{{GeV}}^2$')
            ax.plot(t_plot, ImH(xi, t_plot, params_fit),
                    '--', color=c,
                    label=fr'fit,  $-\,t={t:.1f}\,\mathrm{{GeV}}^2$')
        ax.set_title(rf"$\xi={xi:.2f}$")
        ax.grid(True)

    fig.text(0.5, 0.04, r"$-\,t\ \mathrm{[GeV^2]}$", ha='center')
    fig.text(0.06, 0.5, r"$\mathrm{Im}\,H(\xi,t)$", va='center', rotation='vertical')

    handles, labels = axes[0,0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right', fontsize='small',
               title=r'\textbf{Original vs.\ Fit}')

    fig.suptitle(r"\textbf{Comparison of Original vs.\ Fitted ImH}", y=0.98)
    fig.tight_layout(rect=[0,0.03,1,0.95])

    plt.savefig(outpath)
    print(f"[+] Saved ImH dependence plot to {outpath}")
    print(f"[+] Parsed χ²/ndof = {chi2ndf:.3f}")

if __name__ == "__main__":
    main()