#!/usr/bin/env python3
"""Publication-style RGC e n pi+ kinematic-distribution plots for reviewer response.

Reads the nominal NH3 paper-version ROOT trees for Su22/Fa22/Sp23 and applies
exactly the period- and (xB,-t')-dependent nominal +/-2 sigma Mx2 windows from
channel_selection's final_carbon_assisted_mx2_cuts.json.

Outputs:
  01_kinematic_1d.png          compact 2x3 one-dimensional overview
  02_core_correlations.png     three non-redundant core 2D correlations
  03_phi_correlations.png      Q2, xB, and -t' versus phi with identical phi axis
  04_phi_projections.png       phi projections for all selected events and by xB row
  kinematic_summary.json       counts/ranges and exact inputs used

No asymmetry-extraction quantities are fitted or modified.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Iterable

import matplotlib.pyplot as plt
import numpy as np
import uproot
from matplotlib.colors import LogNorm

PERIODS = ("su22", "fa22", "sp23")
XB_BINS = ((0.10, 0.25), (0.25, 0.35), (0.35, 0.45), (0.45, 0.60))
TP_BINS = ((0.05, 0.25), (0.25, 0.45), (0.45, 0.65),
           (0.65, 0.85), (0.85, 1.05), (1.05, 1.25))

PAPER_DIR = Path("/work/clas12/thayward/CLAS12_exclusive/enpi+/data/pass2/data/paper_versions")
DEFAULT_INPUTS = {
    p: PAPER_DIR / f"rgc_{p}_inb_NH3_epi+_mom_corrections.root" for p in PERIODS
}
DEFAULT_CUTS = Path(
    "../channel_selection/output/channel_selection_mx2_fit_stability/"
    "final_carbon_assisted_cuts/tables/final_carbon_assisted_mx2_cuts.json"
)
DEFAULT_OUT = Path("output/kinematic_distributions_review")
TREE = "PhysicsEvents"
CHUNK = "250 MB"

ALIASES = {
    "xB": ("xB", "x", "xb", "x_b"),
    "Q2": ("Q2", "q2"),
    "tprime": ("tprime", "t_prime", "tp", "tPrime"),
    "Mx2": ("Mx2", "mx2", "Mx2_epi", "Mx2_epip", "missing_mass_squared", "missing_mass2"),
    "phi": ("phi", "phi1", "phi_h", "trento_phi"),
    "W": ("W", "w"),
    "y": ("y", "inelasticity"),
}


def resolve(tree: uproot.behaviors.TTree.TTree, logical: str, required: bool = True) -> str | None:
    names = set(tree.keys())
    for candidate in ALIASES[logical]:
        if candidate in names:
            return candidate
    if required:
        raise KeyError(f"Could not resolve {logical}; tried {ALIASES[logical]}")
    return None


def phi_degrees(values: np.ndarray, branch: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return values
    # Same practical unit test used by the extraction: radian-valued phi has
    # magnitude O(pi), while degree-valued phi is O(180).
    if np.nanmax(np.abs(finite)) <= 2.0 * np.pi + 0.25:
        values = np.degrees(values)
    # Put both [-180,180) and [0,360) conventions onto [0,360).
    return np.mod(values, 360.0)


def bin_index(xb: np.ndarray, mtp: np.ndarray) -> np.ndarray:
    out = np.full(xb.shape, -1, dtype=np.int16)
    for ix, (xl, xh) in enumerate(XB_BINS):
        xm = (xb >= xl) & (xb < xh)
        for it, (tl, th) in enumerate(TP_BINS):
            tm = (mtp >= tl) & (mtp < th)
            out[xm & tm] = ix * len(TP_BINS) + it + 1
    return out


def load_nominal_windows(path: Path) -> dict[str, dict[int, tuple[float, float]]]:
    payload = json.loads(path.read_text())
    if float(payload.get("nominal_sigma_multiple", np.nan)) != 2.0:
        raise RuntimeError(f"Expected nominal 2-sigma windows in {path}")
    windows: dict[str, dict[int, tuple[float, float]]] = {}
    for period in PERIODS:
        windows[period] = {}
        for row in payload["periods"][period]:
            lo, hi = row["nominal"]
            windows[period][int(row["bin_number"])] = (float(lo), float(hi))
        if len(windows[period]) != 24:
            raise RuntimeError(f"Expected 24 nominal windows for {period}")
    return windows


def collect(period: str, path: Path, windows: dict[int, tuple[float, float]]) -> dict[str, np.ndarray]:
    print(f"[load] {period}: {path}", flush=True)
    with uproot.open(path) as f:
        tree = f[TREE]
        branches = {k: resolve(tree, k, required=(k not in {"W", "y"})) for k in ALIASES}
        expressions = [v for v in branches.values() if v is not None]

        acc = {k: [] for k in ("xB", "Q2", "minus_tprime", "Mx2", "phi", "W", "y")}
        seen = selected = 0
        for ichunk, arrays in enumerate(tree.iterate(expressions=expressions, step_size=CHUNK, library="np"), 1):
            xb = np.asarray(arrays[branches["xB"]], dtype=float)
            q2 = np.asarray(arrays[branches["Q2"]], dtype=float)
            mtp = -np.asarray(arrays[branches["tprime"]], dtype=float)
            mx2 = np.asarray(arrays[branches["Mx2"]], dtype=float)
            phi = phi_degrees(arrays[branches["phi"]], branches["phi"])
            w = np.asarray(arrays[branches["W"]], dtype=float) if branches["W"] else np.full(xb.shape, np.nan)
            y = np.asarray(arrays[branches["y"]], dtype=float) if branches["y"] else np.full(xb.shape, np.nan)

            b = bin_index(xb, mtp)
            keep = (b > 0) & np.isfinite(xb) & np.isfinite(q2) & np.isfinite(mtp) & np.isfinite(mx2) & np.isfinite(phi)
            # Apply the exact period/bin-specific nominal Mx2 window.
            for bn in range(1, 25):
                lo, hi = windows[bn]
                m = b == bn
                keep[m] &= (mx2[m] >= lo) & (mx2[m] <= hi)

            # Make the standard DIS requirements explicit if W/y are available.
            keep &= q2 > 1.0
            if branches["W"]:
                keep &= w > 2.0
            if branches["y"]:
                keep &= y < 0.80

            seen += xb.size
            selected += int(np.count_nonzero(keep))
            vals = {"xB": xb, "Q2": q2, "minus_tprime": mtp, "Mx2": mx2,
                    "phi": phi, "W": w, "y": y}
            for key, value in vals.items():
                acc[key].append(value[keep])
            print(f"[load] {period}: chunk {ichunk}; seen={seen:,}; selected={selected:,}", flush=True)

    result = {k: np.concatenate(v) if v else np.empty(0) for k, v in acc.items()}
    result["period"] = np.full(result["xB"].size, period)
    print(f"[load] {period}: DONE; selected={result['xB'].size:,}", flush=True)
    return result


def finite(a: np.ndarray) -> np.ndarray:
    return a[np.isfinite(a)]


def hist2d(ax, x, y, bins, ranges, xlabel, ylabel):
    good = np.isfinite(x) & np.isfinite(y)
    h, xe, ye = np.histogram2d(x[good], y[good], bins=bins, range=ranges)
    positive = h[h > 0]
    norm = LogNorm(vmin=max(1.0, float(positive.min())) if positive.size else 1.0,
                   vmax=max(2.0, float(h.max())))
    mesh = ax.pcolormesh(xe, ye, h.T, shading="auto", norm=norm)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    return mesh


def save(fig, path: Path):
    fig.savefig(path, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"[plot] wrote {path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cuts", type=Path, default=DEFAULT_CUTS)
    ap.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    for p in PERIODS:
        ap.add_argument(f"--{p}", type=Path, default=DEFAULT_INPUTS[p])
    args = ap.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    windows = load_nominal_windows(args.cuts)
    datasets = [collect(p, getattr(args, p), windows[p]) for p in PERIODS]
    data = {}
    for key in datasets[0]:
        data[key] = np.concatenate([d[key] for d in datasets])
    n = data["xB"].size
    print(f"[summary] combined selected sample: {n:,} events", flush=True)

    # 1) Compact 1D overview. No duplicate x/y-swapped panels.
    one_d = [
        ("xB", np.linspace(0.10, 0.60, 51), r"$x_B$"),
        ("Q2", np.linspace(1.0, max(8.0, np.nanpercentile(data["Q2"], 99.7)), 55), r"$Q^2$ (GeV$^2$)"),
        ("minus_tprime", np.linspace(0.05, 1.25, 49), r"$-t'$ (GeV$^2$)"),
        ("phi", np.linspace(0.0, 360.0, 49), r"$\phi$ (deg)"),
    ]
    if np.isfinite(data["W"]).any():
        one_d.append(("W", np.linspace(2.0, max(4.5, np.nanpercentile(data["W"], 99.7)), 50), r"$W$ (GeV)"))
    if np.isfinite(data["y"]).any():
        one_d.append(("y", np.linspace(0.0, 0.80, 49), r"$y$"))
    fig, axes = plt.subplots(2, 3, figsize=(11.5, 7.0))
    for ax, item in zip(axes.flat, one_d):
        key, bins, label = item
        ax.hist(finite(data[key]), bins=bins, histtype="step", linewidth=1.5)
        ax.set_xlabel(label); ax.set_ylabel("Events")
        ax.tick_params(direction="in", top=True, right=True)
    for ax in axes.flat[len(one_d):]: ax.axis("off")
    fig.suptitle(r"RGC $e n \pi^+$ kinematics after nominal exclusivity selection")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, args.output_dir / "01_kinematic_1d.png")

    # 2) Three genuinely distinct core correlations only.
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.1))
    specs = [
        (data["xB"], data["Q2"], (50, 55), ((0.10, 0.60), (1.0, max(8.0, np.nanpercentile(data["Q2"], 99.7)))), r"$x_B$", r"$Q^2$ (GeV$^2$)"),
        (data["xB"], data["minus_tprime"], (50, 48), ((0.10, 0.60), (0.05, 1.25)), r"$x_B$", r"$-t'$ (GeV$^2$)"),
        (data["Q2"], data["minus_tprime"], (55, 48), ((1.0, max(8.0, np.nanpercentile(data["Q2"], 99.7))), (0.05, 1.25)), r"$Q^2$ (GeV$^2$)", r"$-t'$ (GeV$^2$)"),
    ]
    for ax, spec in zip(axes, specs):
        mesh = hist2d(ax, *spec)
        fig.colorbar(mesh, ax=ax, label="Events")
    fig.suptitle("Non-redundant kinematic correlations")
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    save(fig, args.output_dir / "02_core_correlations.png")

    # 3) Direct answer to the reviewer's phi concern: all three observables use
    # the same selected events, the same phi range, and the same phi bin edges.
    fig, axes = plt.subplots(3, 1, figsize=(9.0, 10.0), sharex=True)
    q2max = max(8.0, np.nanpercentile(data["Q2"], 99.7))
    phi_edges = np.linspace(0.0, 360.0, 49)
    ys = [(data["Q2"], np.linspace(1.0, q2max, 55), r"$Q^2$ (GeV$^2$)"),
          (data["xB"], np.linspace(0.10, 0.60, 51), r"$x_B$"),
          (data["minus_tprime"], np.linspace(0.05, 1.25, 49), r"$-t'$ (GeV$^2$)")]
    for ax, (yy, yedges, ylabel) in zip(axes, ys):
        good = np.isfinite(data["phi"]) & np.isfinite(yy)
        h, xe, ye = np.histogram2d(data["phi"][good], yy[good], bins=(phi_edges, yedges))
        pos = h[h > 0]
        mesh = ax.pcolormesh(xe, ye, h.T, shading="auto",
                             norm=LogNorm(vmin=max(1.0, float(pos.min())) if pos.size else 1.0,
                                          vmax=max(2.0, float(h.max()))))
        ax.set_ylabel(ylabel)
        fig.colorbar(mesh, ax=ax, label="Events")
    axes[-1].set_xlabel(r"$\phi$ (deg)")
    axes[-1].set_xlim(0, 360)
    fig.suptitle(r"Common-event-sample correlations with $\phi$")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    save(fig, args.output_dir / "03_phi_correlations.png")

    # 4) Phi projections make acceptance structure at 0/360 immediately visible
    # and show whether it changes primarily with xB rather than implying mutually
    # inconsistent 2D distributions.
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.5), sharex=True, sharey=True)
    for ix, (ax, (xl, xh)) in enumerate(zip(axes.flat, XB_BINS)):
        m = (data["xB"] >= xl) & (data["xB"] < xh)
        counts, edges = np.histogram(data["phi"][m], bins=phi_edges)
        centers = 0.5 * (edges[:-1] + edges[1:])
        # Normalize within xB row so shape, rather than very different statistics,
        # is what the reviewer compares.
        density = counts / counts.sum() if counts.sum() else counts.astype(float)
        ax.step(centers, density, where="mid", linewidth=1.5)
        ax.set_title(fr"${xl:.2f} \leq x_B < {xh:.2f}$  ($N={counts.sum():,}$)")
        ax.set_xlim(0, 360); ax.tick_params(direction="in", top=True, right=True)
    for ax in axes[-1, :]: ax.set_xlabel(r"$\phi$ (deg)")
    for ax in axes[:, 0]: ax.set_ylabel("Fraction of events / bin")
    fig.suptitle(r"$\phi$ acceptance shape in the four analysis $x_B$ intervals")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    save(fig, args.output_dir / "04_phi_projections.png")

    summary = {
        "cuts_json": str(args.cuts.resolve()),
        "nominal_sigma_multiple": 2.0,
        "total_selected_events": int(n),
        "selected_by_period": {p: int(np.count_nonzero(data["period"] == p)) for p in PERIODS},
        "inputs": {p: str(getattr(args, p).resolve()) for p in PERIODS},
        "explicit_dis_cuts": {"Q2_min_gev2": 1.0, "W_min_gev_if_branch_present": 2.0, "y_max_if_branch_present": 0.80},
        "analysis_binning": {"xB": XB_BINS, "minus_tprime_gev2": TP_BINS},
    }
    (args.output_dir / "kinematic_summary.json").write_text(json.dumps(summary, indent=2))
    print(f"[summary] wrote {args.output_dir / 'kinematic_summary.json'}", flush=True)


if __name__ == "__main__":
    main()
