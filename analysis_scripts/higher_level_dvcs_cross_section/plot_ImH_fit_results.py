#!/usr/bin/env python3
"""
plot_ImH_fit_results.py

Load fitted GPD‐H parameters from a fit_results_<TIMESTAMP>.txt file
and compare the original (VGG) vs. fitted Im H(ξ,t) ansatz.

Makes a 2×2 panel for four ξ values, each plotting Im H vs. −t at
three fixed values (0.1, 0.4, 0.7 GeV²). Solid lines = original VGG;
dashed lines = your fit.

Usage:
    python plot_ImH_fit_results.py output/fit_results/fit_results_<TIMESTAMP>.txt
"""

import os, sys, re
import numpy as np
import matplotlib.pyplot as plt

def load_fit(fname):
    """
    Parse the new style fit_results file:
      # values: renormImag alpha0 alpha1 n_val b_val Mm2_val P_val renormReal
      <eight floats>
      # errors:
      <eight floats>
      # chi2 ndof chi2/ndof
      <chi2> <ndof> <chi2/ndof>
    Returns:
      params: dict mapping name→value
      chi2ndf: float
      timestamp: str
    """
    # extract timestamp
    m = re.search(r'fit_results_(\d{8}_\d{6})\.txt$', fname)
    if not m:
        raise RuntimeError(f"Couldn't extract timestamp from '{fname}'")
    timestamp = m.group(1)

    keys = None
    vals = None
    chi2ndf = None

    with open(fname, 'r', encoding='utf-8') as f:
        lines = [l.strip() for l in f if l.strip()]

    for i, line in enumerate(lines):
        if line.startswith("# values:"):
            # pull out the parameter names
            keys = line.split(":",1)[1].split()
            # next line are the values
            vals = list(map(float, lines[i+1].split()))
        if line.startswith("# chi2"):
            parts = lines[i+1].split()
            # parts = [chi2, ndof, chi2/ndof]
            chi2ndf = float(parts[2])

    if keys is None or vals is None:
        raise RuntimeError("Failed to parse '# values:' block in fit file")
    if chi2ndf is None:
        raise RuntimeError("Failed to parse '# chi2' block in fit file")

    if len(keys) != len(vals):
        raise RuntimeError("Number of fit keys != number of fit values")

    params = dict(zip(keys, vals))
    return params, chi2ndf, timestamp

# ──────────────────────────────────────────────────────────────────────────────
# ImH ansatz
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

    # 1) parse fit results
    params_fit, chi2ndf, timestamp = load_fit(fitfile)

    # 2) original VGG defaults
    params_orig = {
        'renormImag': 1.0,
        'alpha0':     0.43,
        'alpha1':     0.85,
        'n_val':      1.35,
        'b_val':      0.4,
        'Mm2_val':    0.64,
        'P_val':      1.0
    }

    # 3) choose ξ panels and t‐points
    xi_vals = [0.05, 0.10, 0.20, 0.30]       # reasonable DVCS skewness
    t_vals  = [0.1, 0.4, 0.7]                # fixed -t [GeV^2]
    colors  = ['tab:blue', 'tab:orange', 'tab:green']

    # 4) prepare output dir
    outdir = "output/plots"
    os.makedirs(outdir, exist_ok=True)
    outname = os.path.join(outdir,
        f"ImH_dependence_{timestamp}.pdf"
    )

    # 5) set style & LaTeX labels
    try:
        plt.style.use('ggplot')
    except OSError:
        pass
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=12)

    # 6) draw 2×2 grid
    fig, axes = plt.subplots(2, 2, figsize=(10,8), sharey=True)
    axes = axes.flatten()

    for ax, xi in zip(axes, xi_vals):
        for t, c in zip(t_vals, colors):
            # plot original
            ax.plot(
                t_vals,
                [ImH(xi, tt, params_orig) for tt in t_vals],
                '-', lw=2, color=c,
                label=rf"orig, $-t={t:.1f}\,\mathrm{{GeV}}^2$"
            )
            # plot fitted
            ax.plot(
                t_vals,
                [ImH(xi, tt, params_fit) for tt in t_vals],
                '--', lw=2, color=c,
                label=rf"fit,  $-t={t:.1f}\,\mathrm{{GeV}}^2$"
            )
        ax.set_title(rf"$\xi = {xi:.2f}$")
        ax.set_xlabel(r"$-t\,[\mathrm{GeV}^2]$")
        ax.grid(True)

    axes[0].set_ylabel(r"$\mathrm{Im}\,H(\xi,t)$")
    axes[2].set_ylabel(r"$\mathrm{Im}\,H(\xi,t)$")

    # single legend
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(
        handles, labels,
        loc='upper right', fontsize='small',
        title=r'\textbf{Original vs. Fit}'
    )

    fig.tight_layout(rect=[0,0,1,0.95])
    fig.suptitle(r"\textbf{Im$H$ Dependence: Original vs.\ Fitted}", y=0.99)

    # 7) save
    fig.savefig(outname)
    print(f"[+] Saved ImH dependence plot to {outname}")
    print(f"[+] Parsed χ²/ndof = {chi2ndf:.3f}")

if __name__ == "__main__":
    main()