#!/usr/bin/env python3
"""
Reproduce the RGC positive-hadron beta-vs-p PID diagnostic before/after cuts.

For each RGC running period, make a 1x2 figure:
  left  : positive FD hadron candidates before the requested PID cuts
  right : after |chi2pid| < 3.5 and 0.5 < p < 5.0 GeV

Input ROOT trees are the calibration trees listed below.  The script searches
for the TTree automatically and accepts the common branch aliases used by the
RGC calibration output.

Run:
    python plot_beta_vs_p_pid.py

Output:
    output/beta_vs_p_pid_su22.png
    output/beta_vs_p_pid_fa22.png
    output/beta_vs_p_pid_sp23.png
"""

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import uproot

INPUTS = {
    "Su22": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/"
                 "rgc_su22_inb_NH3_epi+X_calibration.root"),
    "Fa22": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/"
                 "rgc_fa22_inb_NH3_epi+X_calibration.root"),
    "Sp23": Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/calibration/"
                 "rgc_sp23_inb_NH3_epi+X_calibration.root"),
}

OUTDIR = Path("output")

# Match the original plot's positive-hadron species.
POSITIVE_HADRONS = (211, 321, 2212)

# FD particles in CLAS12 have |status| in [2000, 4000).
FD_STATUS_MIN = 2000
FD_STATUS_MAX = 4000

P_BINS = np.linspace(0.0, 6.0, 241)
BETA_BINS = np.linspace(0.2, 1.2, 201)


def find_tree(root_file):
    """Return the first TTree in the ROOT file."""
    for key, obj in root_file.items(recursive=True):
        if isinstance(obj, uproot.behaviors.TTree.TTree):
            return obj
    raise RuntimeError("No TTree found in ROOT file.")


def branch_name(tree, *aliases):
    names = {str(k).split(";")[0] for k in tree.keys()}
    for name in aliases:
        if name in names:
            return name
    raise KeyError(
        f"None of {aliases} found. Available branches include: "
        + ", ".join(sorted(names)[:50])
    )


def make_period_plot(period, filename):
    print(f"[{period}] reading {filename}")
    with uproot.open(filename) as root_file:
        tree = find_tree(root_file)

        b_pid = branch_name(tree, "particle_pid", "pid")
        b_p = branch_name(tree, "p", "particle_p")
        b_beta = branch_name(tree, "particle_beta", "beta")
        b_chi2 = branch_name(tree, "particle_chi2pid", "chi2pid")
        b_status = branch_name(tree, "particle_status", "status")

        arrays = tree.arrays(
            [b_pid, b_p, b_beta, b_chi2, b_status],
            library="np",
        )

    pid = np.asarray(arrays[b_pid])
    p = np.asarray(arrays[b_p], dtype=float)
    beta = np.asarray(arrays[b_beta], dtype=float)
    chi2pid = np.asarray(arrays[b_chi2], dtype=float)
    status = np.asarray(arrays[b_status])

    finite = (
        np.isfinite(p)
        & np.isfinite(beta)
        & np.isfinite(chi2pid)
    )
    positive_hadron = np.isin(pid, POSITIVE_HADRONS)
    abs_status = np.abs(status)
    fd = (abs_status >= FD_STATUS_MIN) & (abs_status < FD_STATUS_MAX)

    before = finite & positive_hadron & fd
    after = (
        before
        & (np.abs(chi2pid) < 3.5)
        & (p > 0.5)
        & (p < 5.0)
    )

    print(
        f"[{period}] positive FD hadrons: before={before.sum():,}, "
        f"after={after.sum():,} "
        f"({100.0 * after.sum() / max(before.sum(), 1):.2f}% retained)"
    )

    fig, axes = plt.subplots(
        1, 2, figsize=(12, 5), sharex=False, sharey=False,
        constrained_layout=True,
    )

    # Use a common logarithmic normalization so the two panels are directly
    # visually comparable.
    h_before, _, _ = np.histogram2d(
        p[before], beta[before], bins=(P_BINS, BETA_BINS)
    )
    h_after, _, _ = np.histogram2d(
        p[after], beta[after], bins=(P_BINS, BETA_BINS)
    )
    vmax = max(float(h_before.max()), float(h_after.max()), 1.0)

    for ax, mask, title in (
        (axes[0], before, "Before PID cuts"),
        (axes[1], after, r"After $|\chi^2_{\rm PID}|<3.5$, $0.5<p<5.0$ GeV"),
    ):
        h = ax.hist2d(
            p[mask],
            beta[mask],
            bins=(P_BINS, BETA_BINS),
            norm=LogNorm(vmin=1.0, vmax=vmax),
            cmap="turbo",
        )
        ax.set_xlabel(r"$p$ (GeV)")
        ax.set_title(title)
        fig.colorbar(h[3], ax=ax, label="Counts")
    # endfor

    # Keep the first panel as the broad diagnostic view and zoom the
    # post-selection panel onto the accepted pion band.
    axes[0].set_xlim(0.0, 6.0)
    axes[0].set_ylim(0.2, 1.2)
    axes[1].set_xlim(0.5, 5.0)
    axes[1].set_ylim(0.8, 1.1)

    axes[0].set_ylabel(r"$\beta$")
    fig.suptitle(
        f"{period} FD positive hadrons "
        r"($\pi^+$, $K^+$, $p$)"
    )

    outfile = OUTDIR / f"beta_vs_p_pid_{period.lower()}.png"
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"[{period}] wrote {outfile}")


def main():
    OUTDIR.mkdir(parents=True, exist_ok=True)

    for period, filename in INPUTS.items():
        if not filename.is_file():
            raise FileNotFoundError(f"Missing input ROOT file: {filename}")
        make_period_plot(period, filename)
    # endfor


if __name__ == "__main__":
    main()
