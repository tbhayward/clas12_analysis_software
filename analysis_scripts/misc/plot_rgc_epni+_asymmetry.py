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
enpichi2FitsALUsinphi_Su22 = [
    [0.093969228, 1.121102910, 0.103682047],
    [0.168352864, 0.075343496, 0.011124971],
    [0.255015669, 0.110872309, 0.007668837],
    [0.348460186, 0.136797553, 0.008004000],
    [0.441445130, 0.122327908, 0.010538986],
    [0.535418362, 0.095188929, 0.018040875]
]
enpichi2FitsAULsinphi_Su22 = [
    [0.093969228, -0.338857579, 0.194599332],
    [0.168352864, -0.051397484, 0.010823311],
    [0.255015669, -0.003922942, 0.005337434],
    [0.348460186, 0.035992215, 0.005302819],
    [0.441445130, 0.039496051, 0.006732558],
    [0.535418362, 0.031684185, 0.011741457]
]
enpichi2FitsAULsin2phi_Su22 = [
    [0.093969228, -0.933145689, 0.505068832],
    [0.168352864, -0.020302810, 0.022454560],
    [0.255015669, -0.068209117, 0.012126491],
    [0.348460186, -0.086572841, 0.011778306],
    [0.441445130, -0.051343918, 0.014025767],
    [0.535418362, -0.036623704, 0.024819806]
]
enpichi2FitsALL_Su22 = [
    [0.093969228, -0.226411668, 0.322623368],
    [0.168352864, 0.191460695, 0.023617359],
    [0.255015669, 0.302080680, 0.022686092],
    [0.348460186, 0.434065099, 0.029778569],
    [0.441445130, 0.522671691, 0.036871963],
    [0.535418362, 0.634253359, 0.048909323]
]
enpichi2FitsALLcosphi_Su22 = [
    [0.093969228, -0.622700113, 0.579882138],
    [0.168352864, 0.187035564, 0.037593926],
    [0.255015669, 0.067043637, 0.035191562],
    [0.348460186, -0.028525676, 0.044653070],
    [0.441445130, -0.024721894, 0.056791830],
    [0.535418362, -0.041561041, 0.074395483]
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
    r"$ep \rightarrow en\pi^{+}$: Su22, Fa22, Sp23",
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