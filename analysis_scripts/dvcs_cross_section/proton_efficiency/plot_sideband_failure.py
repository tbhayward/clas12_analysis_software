#!/usr/bin/env python3
"""Plot the low-momentum proton sideband failure from saved analysis outputs."""

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--input-dir", type=Path, default=Path("output/proton_efficiency/fa18_inb"),
    help="Directory containing mass_shape_v1 and integrated_fit_v2_dvcs")
parser.add_argument("--mass-cache", type=Path,
                    help="Override mass_shape_v1/9_mass_shape_histograms.npz")
parser.add_argument("--cells-csv", type=Path,
                    help="Override integrated_fit_v2_dvcs/7_theta_p_cells.csv")
parser.add_argument("--output", type=Path,
                    help="Output PNG path")
args = parser.parse_args()
mass_cache = args.mass_cache or args.input_dir/"mass_shape_v1/9_mass_shape_histograms.npz"
cells_csv = args.cells_csv or args.input_dir/"integrated_fit_v2_dvcs/7_theta_p_cells.csv"
output = args.output or args.input_dir/"sideband_failure/1_low_momentum_sideband_failure.png"
for source in (mass_cache, cells_csv):
    if not source.is_file():
        parser.error(f"Missing input: {source}")
MASS = np.load(mass_cache)
with cells_csv.open() as stream:
    rows = list(csv.DictReader(stream))
exact = {angle: next((row for row in rows
                      if float(row["theta_low_deg"]) == angle
                      and float(row["p_low_GeV"]) == .2), None)
         for angle in (50, 60)}
if any(row is None for row in exact.values()):
    parser.error("Missing θ=50–60° or 60–70°, p=0.2–0.4 GeV cell")

edges = MASS["m2_edges"]
centers = (edges[:-1] + edges[1:]) / 2
mp2 = .9382720813**2
left = (mp2 - .45, mp2 - .15)
peak = (mp2 - .15, mp2 + .15)
right = (mp2 + .15, mp2 + .45)
def fraction(values, limits):
    if values.sum() == 0:
        return float("nan")
    return values[(centers >= limits[0]) & (centers < limits[1])].sum()/values.sum()

fig, axes = plt.subplots(2, 2, figsize=(13, 9),
                         gridspec_kw={"height_ratios": [2, 1]},
                         constrained_layout=True)
for column, theta in enumerate((50, 60)):
    ax, note = axes[:, column]
    idx = theta // 10
    data = MASS["data_all"][:, idx, 0]
    mc = MASS["mc_all"][:, idx, 0]
    near = MASS["mc_truth_valid_near_aligned"][:, idx, 0]
    ax.axvspan(*left, color="#ef9a9a", alpha=.28, label="Left sideband")
    ax.axvspan(*peak, color="#90caf9", alpha=.23, label="Peak window")
    ax.axvspan(*right, color="#ef9a9a", alpha=.18, label="Right sideband")
    for values, name, color, style in (
        (data, "Data: all tags", "#1f4e79", "-"),
        (mc, "MC: all tags", "#c05a14", "--"),
        (near, "MC: truth-near proton tags", "#15803d", "-"),
    ):
        if values.sum():
            ax.step(centers, values/values.sum(), where="mid", linewidth=1.8,
                    color=color, linestyle=style, label=name)
    ax.set(xlim=(0, 1.4),
           xlabel="Reconstructed missing mass squared (GeV²)",
           ylabel="Fraction per 0.01 GeV²",
           title=f"Predicted θ={theta}–{theta+10}°, p=0.2–0.5 GeV")
    ax.grid(alpha=.16)
    ax.legend(fontsize=8, loc="upper right")
    ax.text(.03, .96,
            f"Truth-near MC in left sideband: {fraction(near,left):.0%}\n"
            f"Truth-near MC in peak: {fraction(near,peak):.0%}",
            transform=ax.transAxes, va="top", fontsize=9, color="#14532d",
            bbox=dict(boxstyle="round,pad=.3", facecolor="white",
                      edgecolor="#15803d"))

    row = exact[theta]
    p, l, r = (int(row[key]) for key in ("mc_peak", "mc_left", "mc_right"))
    estimate = p - .5*(l+r)
    note.axis("off")
    note.set_title(f"Exact sideband cell: θ={theta}–{theta+10}°, p=0.2–0.4 GeV",
                   fontsize=10.5)
    note.text(.04, .82,
              f"MC peak: {p:,}     left: {l:,}     right: {r:,}\n"
              f"Peak − (left + right)/2 = {estimate:,.1f}",
              transform=note.transAxes, fontsize=12, va="top",
              color="#b91c1c" if estimate <= 0 else "#111827")
    note.text(.04, .34,
              "Negative result: background subtraction is invalid here."
              if estimate <= 0 else "Near-zero result: denominator is unstable.",
              transform=note.transAxes, fontsize=10, va="top",
              color="#b91c1c")
fig.suptitle("High-angle low-momentum sideband failure", fontsize=15)
output.parent.mkdir(parents=True, exist_ok=True)
fig.savefig(output, dpi=180)
print(f"[INPUT] mass_shape={mass_cache} cells={cells_csv}")
print(f"[OUTPUT] {output}")
