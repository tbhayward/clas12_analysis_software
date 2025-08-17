# -*- coding: utf-8 -*-
# Markers-only errorbar plot for A_{UL}^{sin^n φ}; saves PDF to output/enpi+/
# Includes "#endfor" / "#endif" comments per your preference.

from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# --------------------------
# Data (triples: [xB_mean (unused for x), asym, err])
# --------------------------
a1raw = np.array([
    [0.240240319, -0.061361791, 0.012107184],
    [0.240240319, -0.075138978, 0.015014592],
    [0.240240319, -0.065822478, 0.012052638],
    [0.240240319, -0.041226811, 0.014264659],
    [0.240240319, -0.073887825, 0.015092830],
], dtype=float)

a2raw = np.array([
    [0.240240319, -0.294339122, 0.028625956],
    [0.240240319, -0.298589432, 0.028476098],
    [0.240240319, -0.296541576, 0.028384420],
    [0.240240319, -0.224695943, 0.026973792],
    [0.240240319, -0.296004943, 0.028522073],
], dtype=float)

# --------------------------
# Remap x → 1..N, extract y and ey
# --------------------------
N1, N2 = len(a1raw), len(a2raw)
if N1 != N2:
    raise ValueError(f"Dataset size mismatch: n1={N1}, n2={N2}")  # endif
# endif

x = np.arange(1, N1 + 1, dtype=int)
y1, e1 = a1raw[:, 1], a1raw[:, 2]
y2, e2 = a2raw[:, 1], a2raw[:, 2]

# --------------------------
# Figure and axes
# --------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=False)
ax.set_xlim(0.5, len(x) + 0.5)
ax.set_ylim(-0.4, 0.4)
ax.set_xticks(x)
ax.set_xlabel(r"Test $\#$")
# Interpret \sinn\phi as \sin^{n}\phi for proper TeX
ax.set_ylabel(r"$A_{UL}^{\sin^{n}\phi}$")

# Grid and ticks
ax.grid(True, axis="y", linestyle="--", linewidth=0.8, alpha=0.5)
ax.tick_params(direction="in", which="both", top=True, right=True)

# --------------------------
# Plot error bars (markers only, no connecting lines)
# --------------------------
ax.errorbar(
    x, y1, yerr=e1,
    fmt="o", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="red", label=r"$n=1$"
)

ax.errorbar(
    x, y2, yerr=e2,
    fmt="s", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="blue", label=r"$n=2$"
)

# Full border around plot area
for side in ax.spines.values():
    side.set_linewidth(1.2)
# endfor

# Add a thin border around the whole figure canvas
fig_border = plt.Rectangle(
    (0, 0), 1, 1, transform=fig.transFigure, fill=False, linewidth=1.2, edgecolor="black"
)
fig.add_artist(fig_border)

# Legend (top-right)
leg = ax.legend(loc="upper right", frameon=True, framealpha=1.0, fancybox=True, borderpad=0.8)
leg.get_frame().set_linewidth(1.0)

# Tight layout
plt.tight_layout()

# --------------------------
# Save to PDF
# --------------------------
out_dir = Path("output/enpi+/")
out_dir.mkdir(parents=True, exist_ok=True)
pdf_path = out_dir / "AUL_modulations_enpiPlus_markers.pdf"
fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
print(f"Saved: {pdf_path}")