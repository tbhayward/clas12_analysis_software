# -*- coding: utf-8 -*-
# Markers-only errorbar plot for A_{UL}^{sin^n φ}; removes figure-wide border,
# y-axis in [-0.4, 0.1]; adds weighted-constant fits as thin dashed lines.
# Saves PDF to output/enpi+/

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
# Weighted constant fits (weighted means)
# --------------------------
def weighted_constant_fit(y, sigma):
    w = 1.0 / np.maximum(sigma, 1e-12)**2
    mu = np.sum(w * y) / np.sum(w)
    s_mu = 1.0 / np.sqrt(np.sum(w))
    return mu, s_mu

mu1, smu1 = weighted_constant_fit(y1, e1)
mu2, smu2 = weighted_constant_fit(y2, e2)

# --------------------------
# Figure and axes
# --------------------------
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=False)
ax.set_xlim(0.5, len(x) + 0.5)
ax.set_ylim(-0.4, 0.1)  # requested range
ax.set_xticks(x)
ax.set_xlabel(r"Test $\#$")
ax.set_ylabel(r"$F_{UL}^{\sin^{n\phi}}/F_{UU}$")

# Grid and ticks
ax.grid(True, axis="y", linestyle="--", linewidth=0.8, alpha=0.5)
ax.tick_params(direction="in", which="both", top=True, right=True)

# --------------------------
# Plot error bars (markers only, no connecting lines)
# --------------------------
p1 = ax.errorbar(
    x, y1, yerr=e1,
    fmt="o", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="red", label=r"$n=1$"
)
p2 = ax.errorbar(
    x, y2, yerr=e2,
    fmt="s", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="blue", label=r"$n=2$"
)

# Full border around plot area (keep axis spines)
for side in ax.spines.values():
    side.set_linewidth(1.2)
# endfor

# Plot the fitted constants as thin dashed lines across the bin span
xmin, xmax = 0.5, len(x) + 0.5
c1 = ax.hlines(mu1, xmin, xmax, colors="red", linestyles="--", linewidth=1.0, label=r"$\mathrm{const}\ (n=1)$")
c2 = ax.hlines(mu2, xmin, xmax, colors="blue", linestyles="--", linewidth=1.0, label=r"$\mathrm{const}\ (n=2)$")

# Legend (top-right) including constants
leg = ax.legend(
    handles=[p1.lines[0], p2.lines[0], c1, c2],
    labels=[r"$n=1$", r"$n=2$"],
    loc="upper right", frameon=True, framealpha=1.0, fancybox=True, borderpad=0.8
)
leg.get_frame().set_linewidth(1.0)

# Tight layout
plt.tight_layout()

# --------------------------
# Save to PDF
# --------------------------
out_dir = Path("output/enpi+/")
out_dir.mkdir(parents=True, exist_ok=True)
pdf_path = out_dir / "AUL_modulations_enpiPlus_markers_const.pdf"
fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
print(f"Saved: {pdf_path}")


# =========================
# New canvas: F_LL^{cos^n φ} (n=0,1)
# =========================
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# Data (triples: [xB_mean (unused for x), value, err])
f0raw = np.array([
    [0.240240319,  0.489057877, 0.033797988],
    [0.240240319,  0.494490688, 0.038428438],
    [0.240240319,  0.482235425, 0.034154262],
    [0.240240319,  0.466426246, 0.030542406],
    [0.240240319,  0.489266927, 0.038678149],
], dtype=float)

f1raw = np.array([
    [0.240240319, -0.074598801, 0.058106918],
    [0.240240319,  0.051498923, 0.096899984],
    [0.240240319, -0.032504916, 0.058115016],
    [0.240240319, -0.234254151, 0.041499988],
    [0.240240319,  0.048625587, 0.095899916],
], dtype=float)

# Remap x → 1..N, extract y and ey
N0, N1 = len(f0raw), len(f1raw)
if N0 != N1:
    raise ValueError(f"Dataset size mismatch: n0={N0}, n1={N1}")  # endif
# endif

x  = np.arange(1, N0 + 1, dtype=int)
y0, e0 = f0raw[:, 1], f0raw[:, 2]
y1, e1 = f1raw[:, 1], f1raw[:, 2]

# Weighted constant fits (for dashed reference lines)
def weighted_constant_fit(y, sigma):
    w = 1.0 / np.maximum(sigma, 1e-12)**2
    mu = np.sum(w * y) / np.sum(w)
    s_mu = 1.0 / np.sqrt(np.sum(w))
    return mu, s_mu

mu0, smu0 = weighted_constant_fit(y0, e0)
mu1, smu1 = weighted_constant_fit(y1, e1)

# Figure
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=False)
ax.set_xlim(0.5, len(x) + 0.5)
ax.set_ylim(-0.30, 0.80)  # covers both series comfortably
ax.set_xticks(x)
ax.set_xlabel(r"Test $\#$")
ax.set_ylabel(r"$F_{LL}^{\cos^{n\phi}}/F_{UU}$")

# Grid + ticks; full axis border
ax.grid(True, axis="y", linestyle="--", linewidth=0.8, alpha=0.5)
ax.tick_params(direction="in", which="both", top=True, right=True)
for side in ax.spines.values():
    side.set_linewidth(1.2)
# endfor

# Markers only (no connecting lines)
p0 = ax.errorbar(
    x, y0, yerr=e0,
    fmt="o", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="red", label=r"$n=0$"
)
p1 = ax.errorbar(
    x, y1, yerr=e1,
    fmt="s", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="blue", label=r"$n=1$"
)

# Thin dashed constant fits (not included in legend)
xmin, xmax = 0.5, len(x) + 0.5
ax.hlines(mu0, xmin, xmax, colors="red",  linestyles="--", linewidth=1.0)
ax.hlines(mu1, xmin, xmax, colors="blue", linestyles="--", linewidth=1.0)

# Legend (points only)
leg = ax.legend(
    handles=[p0.lines[0], p1.lines[0]],
    labels=[r"$n=0$", r"$n=1$"],
    loc="upper right", frameon=True, framealpha=1.0, fancybox=True, borderpad=0.8
)
leg.get_frame().set_linewidth(1.0)

plt.tight_layout()

# Save
out_dir = Path("output/enpi+/")
out_dir.mkdir(parents=True, exist_ok=True)
pdf_path = out_dir / "FLL_modulations_enpiPlus_markers_const.pdf"
fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
print(f"Saved: {pdf_path}")


# =========================
# New canvas: F_LU^{sin φ}
# =========================
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# Data (triples: [xB_mean (unused for x), value, err])
flu_raw = np.array([
    [0.241878358, 0.097017905, 0.007507156],
    [0.241878358, 0.095457393, 0.013153418],
    [0.241878358, 0.110146974, 0.009050841],
    [0.242012442, 0.103718537, 0.026284315],
    [0.242012442, 0.102247621, 0.062188413],
    [0.242012442, 0.125901585, 0.030662820],
], dtype=float)

# Remap x → 1..N, extract y and ey
x  = np.arange(1, len(flu_raw) + 1, dtype=int)
y  = flu_raw[:, 1]
ey = flu_raw[:, 2]

# Weighted constant fit (reuse your earlier function if present)
def _weighted_constant_fit(y, sigma):
    w = 1.0 / np.maximum(sigma, 1e-12)**2
    mu = np.sum(w * y) / np.sum(w)
    s_mu = 1.0 / np.sqrt(np.sum(w))
    return mu, s_mu

mu, smu = _weighted_constant_fit(y, ey)

# Figure
fig, ax = plt.subplots(figsize=(7.5, 4.5), constrained_layout=False)
ax.set_xlim(0.5, len(x) + 0.5)
ax.set_ylim(-0.02, 0.18)  # adjust if you want tighter/looser bounds
ax.set_xticks(x)
ax.set_xlabel(r"Test $\#$")
ax.set_ylabel(r"$F_{LU}^{\sin\phi}/F_{UU}$")

# Grid + ticks; full axis border
ax.grid(True, axis="y", linestyle="--", linewidth=0.8, alpha=0.5)
ax.tick_params(direction="in", which="both", top=True, right=True)
for side in ax.spines.values():
    side.set_linewidth(1.2)
# endfor

# Markers only (no connecting lines)
p = ax.errorbar(
    x, y, yerr=ey,
    fmt="o", linestyle="none", markersize=6, linewidth=2, capsize=3,
    color="black"
)

# Thin dashed constant fit (not included in legend)
xmin, xmax = 0.5, len(x) + 0.5
ax.hlines(mu, xmin, xmax, colors="black", linestyles="--", linewidth=1.0)

plt.tight_layout()

# Save
out_dir = Path("output/enpi+/")
out_dir.mkdir(parents=True, exist_ok=True)
pdf_path = out_dir / "FLU_modulation_enpiPlus_markers_const.pdf"
fig.savefig(pdf_path, dpi=300, bbox_inches="tight")
print(f"Saved: {pdf_path}")