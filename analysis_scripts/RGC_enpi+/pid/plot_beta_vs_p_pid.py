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

BEAM_ENERGY = {
    "Su22": 10.5473,
    "Fa22": 10.5563,
    "Sp23": 10.5593,
}
M_PROTON = 0.9382720813
M_PION = 0.13957039

# A single broad neutron-peak window for this diagnostic.  The final analysis
# has bin/period-dependent mu +/- 2 sigma cuts; averaging those values gives
# approximately mu ~ 0.88 GeV^2 and sigma ~ 0.07 GeV^2, motivating the rounded
# common window 0.74 < Mx2 < 1.02 GeV^2.
MX2_MIN = 0.74
MX2_MAX = 1.02

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

        b_run = branch_name(tree, "config_run", "runnum", "run")
        b_event = branch_name(tree, "config_event", "evnum", "event")
        b_pid = branch_name(tree, "particle_pid", "pid")
        b_px = branch_name(tree, "particle_px", "px")
        b_py = branch_name(tree, "particle_py", "py")
        b_pz = branch_name(tree, "particle_pz", "pz")
        b_p = branch_name(tree, "p", "particle_p")
        b_beta = branch_name(tree, "particle_beta", "beta")
        b_chi2 = branch_name(tree, "particle_chi2pid", "chi2pid")
        b_status = branch_name(tree, "particle_status", "status")

        arrays = tree.arrays(
            [b_run, b_event, b_pid, b_px, b_py, b_pz,
             b_p, b_beta, b_chi2, b_status],
            library="np",
        )

    pid = np.asarray(arrays[b_pid])
    p = np.asarray(arrays[b_p], dtype=float)
    beta = np.asarray(arrays[b_beta], dtype=float)
    chi2pid = np.asarray(arrays[b_chi2], dtype=float)
    status = np.asarray(arrays[b_status])
    run = np.asarray(arrays[b_run], dtype=np.int64)
    event = np.asarray(arrays[b_event], dtype=np.int64)
    px = np.asarray(arrays[b_px], dtype=float)
    py = np.asarray(arrays[b_py], dtype=float)
    pz = np.asarray(arrays[b_pz], dtype=float)

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

    # Match every positive-hadron row to the scattered electron from the same
    # (run,event).  The calibration tree contains one row per reconstructed
    # particle, so the event identifiers provide an exact event-level join.
    # Treat every positive-hadron candidate under the pion mass hypothesis when
    # forming ep -> e' pi+ X missing mass; this avoids using the assigned
    # hadron PID itself to manufacture the exclusivity discrimination.
    key = (run.astype(np.uint64) << np.uint64(32)) | (
        event.astype(np.uint64) & np.uint64(0xffffffff)
    )
    electron_mask = (pid == 11) & np.isfinite(px) & np.isfinite(py) & np.isfinite(pz)
    e_idx = np.flatnonzero(electron_mask)
    e_keys = key[e_idx]
    order = np.argsort(e_keys)
    e_keys_sorted = e_keys[order]
    e_idx_sorted = e_idx[order]

    h_idx = np.flatnonzero(after)
    pos = np.searchsorted(e_keys_sorted, key[h_idx])
    matched = pos < len(e_keys_sorted)
    safe_pos = np.minimum(pos, max(len(e_keys_sorted) - 1, 0))
    if len(e_keys_sorted):
        matched &= e_keys_sorted[safe_pos] == key[h_idx]
    else:
        matched[:] = False
    # endif

    matched_h = h_idx[matched]
    matched_e = e_idx_sorted[pos[matched]]

    Ebeam = BEAM_ENERGY[period]
    epx, epy, epz = px[matched_e], py[matched_e], pz[matched_e]
    Ee = np.sqrt(epx**2 + epy**2 + epz**2)  # electron mass negligible here
    qx, qy, qz = -epx, -epy, Ebeam - epz
    nu = Ebeam - Ee
    q2_four = nu**2 - (qx**2 + qy**2 + qz**2)
    Q2 = -q2_four
    W2 = M_PROTON**2 + 2.0 * M_PROTON * nu - Q2
    W = np.sqrt(np.maximum(W2, 0.0))
    y = nu / Ebeam

    hpx, hpy, hpz = px[matched_h], py[matched_h], pz[matched_h]
    hp2 = hpx**2 + hpy**2 + hpz**2
    Eh_pi = np.sqrt(hp2 + M_PION**2)
    mx_E = M_PROTON + nu - Eh_pi
    mx_px = qx - hpx
    mx_py = qy - hpy
    mx_pz = qz - hpz
    Mx2 = mx_E**2 - (mx_px**2 + mx_py**2 + mx_pz**2)

    exclusive_local = (
        (W > 2.0)
        & (y < 0.8)
        & (Mx2 > MX2_MIN)
        & (Mx2 < MX2_MAX)
    )
    exclusive_idx = matched_h[exclusive_local]
    exclusive = np.zeros(len(pid), dtype=bool)
    exclusive[exclusive_idx] = True

    print(
        f"[{period}] positive FD hadrons: before={before.sum():,}, "
        f"PID={after.sum():,}, matched-to-e={len(matched_h):,}, "
        f"PID+DIS+Mx2={exclusive.sum():,}"
    )

    fig, axes = plt.subplots(
        1, 3, figsize=(17, 5), sharex=False, sharey=False,
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
    h_exclusive, _, _ = np.histogram2d(
        p[exclusive], beta[exclusive], bins=(P_BINS, BETA_BINS)
    )
    vmax = max(
        float(h_before.max()), float(h_after.max()),
        float(h_exclusive.max()), 1.0
    )

    for ax, mask, title in (
        (axes[0], before, "Before PID cuts"),
        (axes[1], after, r"After $|\chi^2_{\rm PID}|<3.5$, $0.5<p<5.0$ GeV"),
        (
            axes[2], exclusive,
            rf"PID + $W>2$, $y<0.8$, ${MX2_MIN:.2f}<M_X^2<{MX2_MAX:.2f}$ GeV$^2$"
        ),
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
    axes[2].set_xlim(0.5, 5.0)
    axes[2].set_ylim(0.8, 1.1)

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
