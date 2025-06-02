#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Hard-coded data for all three run periods
#
# Note: Fill in the Fa22 and Sp23 lists once you have those results.
# -----------------------------------------------------------------------------

# RGC Su22 data
# Updated Su22 data for the new iteration
enpichi2FitsALUsinphi_Su22 = [
    [0.094266731, 1.184739999, 0.126228827],
    [0.168622194, 0.082358508, 0.012946927],
    [0.254887003, 0.115297657, 0.008854289],
    [0.348247742, 0.144713802, 0.009273995],
    [0.441277333, 0.117431442, 0.012353808],
    [0.535009222, 0.094723130, 0.021626204]
]

enpichi2FitsAULsinphi_Su22 = [
    [0.094266731, -0.184631847, 0.191657162],
    [0.168622194, -0.044874584, 0.009740722],
    [0.254887003, -0.003269898, 0.004667702],
    [0.348247742,  0.026541647, 0.004654548],
    [0.441277333,  0.034840513, 0.006110932],
    [0.535009222,  0.037791240, 0.010864224]
]

enpichi2FitsAULsin2phi_Su22 = [
    [0.094266731, -0.571432390, 0.500337273],
    [0.168622194, -0.019945862, 0.019491545],
    [0.254887003, -0.056427290, 0.010940032],
    [0.348247742, -0.072941551, 0.010565592],
    [0.441277333, -0.047906347, 0.012682486],
    [0.535009222, -0.056154311, 0.022573166]
]

enpichi2FitsALL_Su22 = [
    [0.094266731, -0.306180156, 0.310512623],
    [0.168622194,  0.143004808, 0.020412174],
    [0.254887003,  0.247780646, 0.020903540],
    [0.348247742,  0.360906216, 0.028176244],
    [0.441277333,  0.455831998, 0.035994583],
    [0.535009222,  0.561823491, 0.047924768]
]

enpichi2FitsALLcosphi_Su22 = [
    [0.094266731, -0.527234710, 0.538136648],
    [0.168622194,  0.134872856, 0.032581310],
    [0.254887003,  0.046520308, 0.032689097],
    [0.348247742, -0.009177298, 0.043351783],
    [0.441277333,  0.006301889, 0.055863447],
    [0.535009222, -0.033540762, 0.074341147]
]

# RGC Fa22 data (PLACEHOLDER: Replace with your actual Fa22 results)
enpichi2FitsALUsinphi_Fa22 = []
enpichi2FitsAULsinphi_Fa22 = []
enpichi2FitsAULsin2phi_Fa22 = []
enpichi2FitsALL_Fa22 = []
enpichi2FitsALLcosphi_Fa22 = []

# RGC Sp23 data (PLACEHOLDER: Replace with your actual Sp23 results)
enpichi2FitsALUsinphi_Sp23 = []
enpichi2FitsAULsinphi_Sp23 = []
enpichi2FitsAULsin2phi_Sp23 = []
enpichi2FitsALL_Sp23 = []
enpichi2FitsALLcosphi_Sp23 = []

# -----------------------------------------------------------------------------
# Organize data by run period
# -----------------------------------------------------------------------------
periods = {
    "Su22": {
        "ALUsinphi": enpichi2FitsALUsinphi_Su22,
        "AULsinphi": enpichi2FitsAULsinphi_Su22,
        "AULsin2phi": enpichi2FitsAULsin2phi_Su22,
        "ALL_n0": enpichi2FitsALL_Su22,
        "ALLcosphi": enpichi2FitsALLcosphi_Su22
    },
    "Fa22": {
        "ALUsinphi": enpichi2FitsALUsinphi_Fa22,
        "AULsinphi": enpichi2FitsAULsinphi_Fa22,
        "AULsin2phi": enpichi2FitsAULsin2phi_Fa22,
        "ALL_n0": enpichi2FitsALL_Fa22,
        "ALLcosphi": enpichi2FitsALLcosphi_Fa22
    },
    "Sp23": {
        "ALUsinphi": enpichi2FitsALUsinphi_Sp23,
        "AULsinphi": enpichi2FitsAULsinphi_Sp23,
        "AULsin2phi": enpichi2FitsAULsin2phi_Sp23,
        "ALL_n0": enpichi2FitsALL_Sp23,
        "ALLcosphi": enpichi2FitsALLcosphi_Sp23
    }
}

# Colors for each run period
colors = {
    "Su22": "tab:blue",
    "Fa22": "tab:orange",
    "Sp23": "tab:green"
}

# -----------------------------------------------------------------------------
# Convert each list in each period to NumPy arrays
# -----------------------------------------------------------------------------
for p in periods:
    for key in periods[p]:
        data_list = periods[p][key]
        if len(data_list) > 0:
            arr = np.array(data_list)
            x = arr[:, 0]
            y = arr[:, 1]
            yerr = arr[:, 2]
            periods[p][key] = {"x": x, "y": y, "yerr": yerr}
        else:
            periods[p][key] = None
        # endif
    # end for key
# end for p

# -----------------------------------------------------------------------------
# Create output directory if it does not exist
# -----------------------------------------------------------------------------
out_dir = os.path.join("output", "enpi+")
if not os.path.isdir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
# endif

# -----------------------------------------------------------------------------
# Plotting: 1×3 figure, all three run periods on each subplot
# -----------------------------------------------------------------------------
plt.figure(figsize=(15, 5))
plt.suptitle(
    r"$ep \rightarrow en\pi^{+}$, $|t| < 1$, $M_{x}^{2} < 1.1\ (\mathrm{GeV}^{2})$",
    fontsize=16,
    y=0.96
)

# Increase base font size for axes labels
label_fontsize = 13

# -------------------------
# Subplot 1: ALU sinφ (all 3 periods)
# -------------------------
ax1 = plt.subplot(1, 3, 1)

for p in ["Su22", "Fa22", "Sp23"]:
    data = periods[p]["ALUsinphi"]
    if data is not None:
        ax1.errorbar(
            data["x"],
            data["y"],
            yerr=data["yerr"],
            fmt="o",
            color=colors[p],
            ecolor=colors[p],
            capsize=3,
            label=p
        )
    # endif
# end for

ax1.set_xlim(0, 0.7)
ax1.set_ylim(-0.2, 0.2)
ax1.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax1.set_ylabel(r"$F_{LU}^{\sin\phi}/F_{UU}$", fontsize=label_fontsize)
ax1.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax1.grid(True, linestyle="--", alpha=0.6)

# Legend for run periods (Su22, Fa22, Sp23)
legend1 = ax1.legend(
    title="Run Period",
    frameon=True,
    edgecolor="black",
    fontsize=11,
    title_fontsize=12
)
legend1.get_frame().set_alpha(0.9)

# -------------------------
# Subplot 2: AUL sinφ (n=1, open) & sin2φ (n=2, closed) for all periods
# -------------------------
ax2 = plt.subplot(1, 3, 2, sharex=ax1)

for p in ["Su22", "Fa22", "Sp23"]:
    d1 = periods[p]["AULsinphi"]
    d2 = periods[p]["AULsin2phi"]
    if d1 is not None:
        ax2.errorbar(
            d1["x"],
            d1["y"],
            yerr=d1["yerr"],
            fmt="o",
            mfc="none",        # open circle
            mec=colors[p],
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=1"
        )
    # endif
    if d2 is not None:
        ax2.errorbar(
            d2["x"],
            d2["y"],
            yerr=d2["yerr"],
            fmt="o",
            color=colors[p],   # filled circle
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=2"
        )
    # endif
# end for

ax2.set_xlim(0, 0.7)
ax2.set_ylim(-0.2, 0.2)
ax2.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax2.set_ylabel(r"$F_{UL}^{\sin\,n\phi}/F_{UU}$", fontsize=label_fontsize)
ax2.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax2.grid(True, linestyle="--", alpha=0.6)

# Legend for n-values
legend_n2 = ax2.legend(
    handles=[
        Line2D([0], [0], marker='o', mfc='none', mec='black', linestyle='', label='n=1'),
        Line2D([0], [0], marker='o', color='black', linestyle='', label='n=2')
    ],
    title="Harmonic n",
    frameon=True,
    edgecolor="black",
    loc='upper right',
    fontsize=11,
    title_fontsize=12
)
legend_n2.get_frame().set_alpha(0.9)

# A second legend for run periods (placed in lower right)
legend_runs2 = ax2.legend(
    handles=[
        Line2D([0], [0], marker='o', color=colors["Su22"], linestyle='', label='Su22'),
        Line2D([0], [0], marker='o', color=colors["Fa22"], linestyle='', label='Fa22'),
        Line2D([0], [0], marker='o', color=colors["Sp23"], linestyle='', label='Sp23')
    ],
    title="Run Period",
    frameon=True,
    edgecolor="black",
    loc='lower right',
    fontsize=11,
    title_fontsize=12
)
legend_runs2.get_frame().set_alpha(0.9)

# -------------------------
# Subplot 3: ALL n=0 (open) & cosφ (n=1, closed) for all periods
# -------------------------
ax3 = plt.subplot(1, 3, 3, sharex=ax1)

for p in ["Su22", "Fa22", "Sp23"]:
    d0 = periods[p]["ALL_n0"]
    d1 = periods[p]["ALLcosphi"]
    if d0 is not None:
        ax3.errorbar(
            d0["x"],
            d0["y"],
            yerr=d0["yerr"],
            fmt="o",
            mfc="none",        # open circle for n=0
            mec=colors[p],
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=0"
        )
    # endif
    if d1 is not None:
        ax3.errorbar(
            d1["x"],
            d1["y"],
            yerr=d1["yerr"],
            fmt="o",
            color=colors[p],   # filled circle for n=1
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=1"
        )
    # endif
# end for

ax3.set_xlim(0, 0.7)
ax3.set_ylim(-0.8, 0.8)
ax3.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax3.set_ylabel(r"$F_{LL}^{\cos\,n\phi}/F_{UU}$", fontsize=label_fontsize)
ax3.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax3.grid(True, linestyle="--", alpha=0.6)

# Legend for n-values
legend_n3 = ax3.legend(
    handles=[
        Line2D([0], [0], marker='o', mfc='none', mec='black', linestyle='', label='n=0'),
        Line2D([0], [0], marker='o', color='black', linestyle='', label='n=1')
    ],
    title="Harmonic n",
    frameon=True,
    edgecolor="black",
    loc='upper right',
    fontsize=11,
    title_fontsize=12
)
legend_n3.get_frame().set_alpha(0.9)

# Second legend for run periods (lower right)
legend_runs3 = ax3.legend(
    handles=[
        Line2D([0], [0], marker='o', color=colors["Su22"], linestyle='', label='Su22'),
        Line2D([0], [0], marker='o', color=colors["Fa22"], linestyle='', label='Fa22'),
        Line2D([0], [0], marker='o', color=colors["Sp23"], linestyle='', label='Sp23')
    ],
    title="Run Period",
    frameon=True,
    edgecolor="black",
    loc='lower right',
    fontsize=11,
    title_fontsize=12
)
legend_runs3.get_frame().set_alpha(0.9)

plt.tight_layout(rect=[0, 0, 1, 0.93])

# Save the figure under output/enpi+/
output_filename = os.path.join(out_dir, "rgc_enpi+_AllPeriods.pdf")
plt.savefig(output_filename)
print(f"Plot saved to '{output_filename}'")