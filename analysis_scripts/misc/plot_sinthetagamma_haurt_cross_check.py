#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Quick cross-check plot:
Left  panel:  F_UL^{sin(phi)}/F_UU  vs  sin(theta_gamma)
Right panel:  A_tg^{sin(phi)}       vs  sin(theta_gamma)

Saves:
  output/enpi+/sinthetagamma_harut_cross_check.pdf
"""

import os
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────
# Input data (x = sin(theta_gamma), y, err)
# ─────────────────────────────────────────────────────────────────────
AULsinphi_harut1 = [
    (0.092612812, -0.052871958, 0.043540426),
    (0.130241094, -0.019239326, 0.006427457),
    (0.175812177, -0.005511188, 0.005219611),
    (0.224199380,  0.001744197, 0.004392396),
    (0.271114953,  0.016805510, 0.005184247),
]
AULsinphi_harut2 = [
    (0.092183435, -0.128052544, 0.030232566),
    (0.128957361, -0.076749577, 0.011913601),
    (0.176173908, -0.056635366, 0.008055780),
    (0.224994591, -0.019165618, 0.006426446),
    (0.271997590, -0.001764200, 0.007751942),
]
AULsinphi_harut3 = [
    (0.092105726,  0.007148937, 0.055761073),
    (0.128731779,  0.006518379, 0.010708090),
    (0.175886625,  0.000708323, 0.009987674),
    (0.224624894,  0.040006724, 0.007578386),
    (0.271768931,  0.066407147, 0.010135559),
]

ATGsinphi_harut1 = [
    (0.092612812,  0.049201099, 0.035113873),
    (0.130241094,  0.044869614, 0.019217713),
    (0.175812177,  0.005290574, 0.015876026),
    (0.224199380,  0.001368505, 0.009006246),
    (0.271114953,  0.010644864, 0.008452191),
]
ATGsinphi_harut2 = [
    (0.092183435,  0.114796057, 0.022594542),
    (0.128957361, -0.007801625, 0.018673315),
    (0.176173908, -0.023095460, 0.024630955),
    (0.224994591,  0.005236748, 0.013931279),
    (0.271997590,  0.002978469, 0.013140461),
]
ATGsinphi_harut3 = [
    (0.092105726, -0.043533569, 0.062025788),
    (0.128731779, -0.017256251, 0.015421493),
    (0.175886625,  0.027435918, 0.023630164),
    (0.224624894, -0.006292719, 0.003934911),
    (0.271768931,  0.001489976, 0.009015100),
]

# Labels, markers, and colors for the three t' ranges
series_meta = [
    ("0.05 < -t < 0.35", "o", "tab:blue"),   # circles
    ("0.35 < -t < 0.85", "s", "tab:orange"), # squares
    ("0.85 < -t < 1.35", "^", "tab:green"),  # triangles
]

# Group the series for easy iteration
AULsinphi_series = [AULsinphi_harut1, AULsinphi_harut2, AULsinphi_harut3]
ATGsinphi_series = [ATGsinphi_harut1, ATGsinphi_harut2, ATGsinphi_harut3]

# ─────────────────────────────────────────────────────────────────────
# Plot
# ─────────────────────────────────────────────────────────────────────
fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 4.5), sharex=True)

def plot_panel(ax, series_list, y_label):
    for (label, marker, color), series in zip(series_meta, series_list):
        x = [t[0] for t in series]
        y = [t[1] for t in series]
        e = [t[2] for t in series]
        ax.errorbar(
            x, y, yerr=e,
            fmt=marker, color=color, ecolor=color,
            markersize=5.0, capsize=3, linestyle="None",
            label=label
        )
    ax.set_xlabel(r"$\sin\theta_\gamma$")
    ax.set_ylabel(y_label)
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(frameon=True, edgecolor="black", fontsize=10)

# Left: F_UL^{sin phi}/F_UU
plot_panel(axL, AULsinphi_series, r"$F_{UL}^{\sin\phi}/F_{UU}$")

# Right: A_tg^{sin phi}
plot_panel(axR, ATGsinphi_series, r"$A_{tg}^{\sin\phi}$")

plt.tight_layout()

# Ensure output directory exists and save
out_path = "output/enpi+/sinthetagamma_harut_cross_check.pdf"
os.makedirs(os.path.dirname(out_path), exist_ok=True)
plt.savefig(out_path)
plt.close()
print(f"Saved: {out_path}")