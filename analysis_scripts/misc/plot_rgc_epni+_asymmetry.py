#!/usr/bin/env python3

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# -----------------------------------------------------------------------------
# Hard-coded data for all three run periods
#
# Updated Su22, Fa22, and Sp23 data
# -----------------------------------------------------------------------------

# RGC Su22 data (new iteration)
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

# RGC Fa22 data (provided results)
enpichi2FitsALUsinphi_Fa22 = [
    [0.094825648, 1.221850467, 0.095501230],
    [0.168249360, 0.081525475, 0.009989312],
    [0.255099359, 0.118300134, 0.007065332],
    [0.348334763, 0.140055462, 0.007378878],
    [0.441201661, 0.140928523, 0.009773392],
    [0.534977263, 0.110737425, 0.016852914]
]

enpichi2FitsAULsinphi_Fa22 = [
    [0.094825648, -0.376535201, 0.252579804],
    [0.168249360, -0.011194402, 0.014338757],
    [0.255099359,  0.030969050, 0.009579115],
    [0.348334763,  0.095187183, 0.012142933],
    [0.441201661,  0.127867020, 0.015579627],
    [0.534977263,  0.070995111, 0.017459805]
]

enpichi2FitsAULsin2phi_Fa22 = [
    [0.094825648, -0.566783688, 0.642963670],
    [0.168249360, -0.064458370, 0.031522996],
    [0.255099359, -0.129267858, 0.020633944],
    [0.348334763, -0.119990899, 0.022631878],
    [0.441201661, -0.113559008, 0.029206404],
    [0.534977263, -0.094101903, 0.035253298]
]

enpichi2FitsALL_Fa22 = [
    [0.094966860, 0.245869621, 0.274476863],
    [0.168066023, 0.011443390, 0.019760262],
    [0.255206095, 0.013521409, 0.018272795],
    [0.348425967, 0.024734125, 0.018401430],
    [0.441160429, 0.009835744, 0.022337638],
    [0.535321960, 0.039381503, 0.030368159]
]

enpichi2FitsALLcosphi_Fa22 = [
    [0.094966860,  0.149000765, 0.501680893],
    [0.168066023, -0.006541656, 0.031406204],
    [0.255206095,  0.025607394, 0.024171661],
    [0.348425967, -0.001640641, 0.022943155],
    [0.441160429, -0.001361157, 0.029325971],
    [0.535321960,  0.025618792, 0.043025987]
]

# RGC Sp23 data (newly provided)
enpichi2FitsALUsinphi_Sp23 = [
    [0.091314751, 0.060825979, 0.081001615],
    [0.166606555, 0.070263513, 0.012697034],
    [0.251219557, 0.107175756, 0.009614248],
    [0.346388218, 0.141913308, 0.010857952],
    [0.441235552, 0.118073067, 0.014919450],
    [0.535187978, 0.115828125, 0.026600146]
]

enpichi2FitsAULsinphi_Sp23 = [
    [0.091314751, -0.162890858, 0.168725297],
    [0.166606555, -0.048730818, 0.014200007],
    [0.251219557,  0.021636751, 0.008788846],
    [0.346388218,  0.061416533, 0.011107487],
    [0.441235552,  0.083394476, 0.015168372],
    [0.535187978,  0.054499349, 0.024740870]
]

enpichi2FitsAULsin2phi_Sp23 = [
    [0.091314751, -0.457165997, 0.391519108],
    [0.166606555, -0.085380031, 0.031549856],
    [0.251219557, -0.144550462, 0.021494309],
    [0.346388218, -0.136857506, 0.023524901],
    [0.441235552, -0.111260784, 0.031197954],
    [0.535187978, -0.028997084, 0.053477696]
]

enpichi2FitsALL_Sp23 = [
    [0.091314751, -0.304360525, 0.346625601],
    [0.166606555,  0.225550203, 0.041348592],
    [0.251219557,  0.467703808, 0.049215221],
    [0.346388218,  0.682918822, 0.066882965],
    [0.441235552,  0.755524113, 0.077509172],
    [0.535187978,  1.030950644, 0.113721860]
]

enpichi2FitsALLcosphi_Sp23 = [
    [0.091314751, -0.201466363, 0.561259035],
    [0.166606555,  0.210346307, 0.063563879],
    [0.251219557,  0.113239109, 0.078644153],
    [0.346388218, -0.020487997, 0.105744642],
    [0.441235552, -0.090960470, 0.117496619],
    [0.535187978, -0.216175985, 0.175101458]
]

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
    r"$ep \rightarrow en\pi^{+}$, $|t| < 1$, $0.75 < M_{x}^{2} < 1.05\ (\mathrm{GeV}^{2})$",
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
ax1.set_ylim(-0.3, 0.3)
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
            mfc="none",        # open circle for n=1
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
            color=colors[p],   # filled circle for n=2
            ecolor=colors[p],
            capsize=3,
            label=f"{p}, n=2"
        )
    # endif
# end for

ax2.set_xlim(0, 0.7)
ax2.set_ylim(-0.3, 0.3)
ax2.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax2.set_ylabel(r"$F_{UL}^{\sin\,n\phi}/F_{UU}$", fontsize=label_fontsize)
ax2.axhline(0, color="black", linestyle="--", linewidth=1.2)
ax2.grid(True, linestyle="--", alpha=0.6)

# Legend for harmonic n (subplot 2)
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
ax2.add_artist(legend_n2)

# Legend for run periods (subplot 2)
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

# Legend for harmonic n (subplot 3)
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
ax3.add_artist(legend_n3)

# Legend for run periods (subplot 3)
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

# Save the asymmetry figure under output/enpi+/
asymmetry_filename = os.path.join(out_dir, "rgc_enpi+_AllPeriods.pdf")
plt.savefig(asymmetry_filename)
print(f"Asymmetry plot saved to '{asymmetry_filename}'")

# -----------------------------------------------------------------------------
# New functionality: Plot dilution factor vs xB for Su22 and Fa22
# -----------------------------------------------------------------------------

# Dilution factor data for Su22
x_Su22 = np.array([row[0] for row in enpichi2FitsALUsinphi_Su22])
dil_Su22 = np.array([0.711119, 0.376873, 0.390790, 0.401177, 0.410563, 0.416077])
dil_err_Su22 = np.array([0.173083, 0.0175818, 0.00942491, 0.00801711, 0.0108668, 0.0225813])

# Dilution factor data for Fa22
x_Fa22 = np.array([row[0] for row in enpichi2FitsALUsinphi_Fa22])
dil_Fa22 = np.array([0.374532, 0.406101, 0.386122, 0.397153, 0.415526, 0.428239])
dil_err_Fa22 = np.array([0.117283, 0.00608866, 0.00354573, 0.00301186, 0.00400862, 0.00824607])

plt.figure(figsize=(6, 5))
plt.errorbar(
    x_Su22,
    dil_Su22,
    yerr=dil_err_Su22,
    fmt="o",
    color=colors["Su22"],
    ecolor=colors["Su22"],
    capsize=3,
    label="Su22"
)
plt.errorbar(
    x_Fa22,
    dil_Fa22,
    yerr=dil_err_Fa22,
    fmt="o",
    color=colors["Fa22"],
    ecolor=colors["Fa22"],
    capsize=3,
    label="Fa22"
)

plt.xlabel(r"$x_{B}$", fontsize=label_fontsize)
plt.ylabel(r"$D_{f}$", fontsize=label_fontsize)
plt.ylim(0, 1)
plt.xlim(0, 0.7)
plt.grid(True, linestyle="--", alpha=0.6)
plt.legend(frameon=True, edgecolor="black", fontsize=11, title="Run Period", title_fontsize=12)

plt.tight_layout()
dilution_filename = os.path.join(out_dir, "dilution_vs_xB.pdf")
plt.savefig(dilution_filename)
print(f"Dilution factor plot saved to '{dilution_filename}'")



# -----------------------------------------------------------------------------
# New functionality: Plot mean Q², W, z, and x_F vs x_B for Su22 and Sp23
# -----------------------------------------------------------------------------

# Mean kinematic values for Su22
xB_Su22 = np.array([0.094, 0.169, 0.255, 0.348, 0.441, 0.535])
Q2_Su22 = np.array([1.325, 1.862, 2.240, 2.592, 3.374, 4.589])
W_Su22  = np.array([3.689, 3.179, 2.722, 2.381, 2.256, 2.197])
z_Su22  = np.array([0.972, 0.971, 0.955, 0.932, 0.919, 0.912])
xF_Su22 = np.array([0.883, 0.862, 0.812, 0.756, 0.745, 0.757])

# Mean kinematic values for Sp23
xB_Sp23 = np.array([0.091, 0.167, 0.251, 0.346, 0.441, 0.535])
Q2_Sp23 = np.array([1.132, 1.502, 1.864, 2.513, 3.379, 4.607])
W_Sp23  = np.array([3.487, 2.896, 2.517, 2.356, 2.258, 2.199])
z_Sp23  = np.array([0.961, 0.960, 0.943, 0.929, 0.919, 0.912])
xF_Sp23 = np.array([0.854, 0.824, 0.771, 0.750, 0.745, 0.757])

# Create a 2×2 figure
fig, axes = plt.subplots(2, 2, figsize=(10, 8))

# Common plot settings
marker_Su22 = 'o'
marker_Sp23 = 's'
color_Su22 = colors["Su22"]
color_Sp23 = colors["Sp23"]
label_fontsize = 12
tick_fontsize = 10

# Top-left: Q² vs x_B
ax = axes[0, 0]
ax.plot(xB_Su22, Q2_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='-', label='Su22')
ax.plot(xB_Sp23, Q2_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='--', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle Q^{2} \rangle\ \mathrm{(GeV^{2})}$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Top-right: W vs x_B
ax = axes[0, 1]
ax.plot(xB_Su22, W_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='-', label='Su22')
ax.plot(xB_Sp23, W_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='--', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle W \rangle\ \mathrm{(GeV)}$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Bottom-left: z vs x_B
ax = axes[1, 0]
ax.plot(xB_Su22, z_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='-', label='Su22')
ax.plot(xB_Sp23, z_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='--', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle z \rangle$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.set_ylim(0.90, 1.0)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Bottom-right: x_F vs x_B
ax = axes[1, 1]
ax.plot(xB_Su22, xF_Su22, marker_Su22, color=color_Su22, markersize=6, linestyle='-', label='Su22')
ax.plot(xB_Sp23, xF_Sp23, marker_Sp23, color=color_Sp23, markersize=6, linestyle='--', label='Sp23')
ax.set_xlabel(r"$x_{B}$", fontsize=label_fontsize)
ax.set_ylabel(r"$\langle x_{F} \rangle$", fontsize=label_fontsize)
ax.set_xlim(0.08, 0.56)
ax.grid(True, linestyle="--", alpha=0.5)
ax.tick_params(axis='both', labelsize=tick_fontsize)

# Add a single legend for the entire figure
handles = [
    Line2D([0], [0], marker=marker_Su22, color=color_Su22, linestyle='-', label='Su22', markersize=6),
    Line2D([0], [0], marker=marker_Sp23, color=color_Sp23, linestyle='--', label='Sp23', markersize=6)
]
fig.legend(
    handles=handles,
    loc='upper center',
    ncol=2,
    frameon=True,
    edgecolor="black",
    fontsize=11,
    title="Data Set",
    title_fontsize=12
)

plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save the kinematic comparison figure
kinematic_filename = os.path.join(out_dir, "kinematic_comparison_Su22_Sp23.pdf")
plt.savefig(kinematic_filename)
print(f"Kinematic comparison plot saved as '{kinematic_filename}'")