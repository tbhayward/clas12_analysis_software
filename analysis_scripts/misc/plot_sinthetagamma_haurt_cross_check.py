#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Cross-check plot (1×3):
  Left   : F_UL^{sin(phi)}/F_UU  vs  sin(theta_gamma)
  Middle : F_UL^{sin(2phi)}/F_UU vs  sin(theta_gamma)
  Right  : A_tg^{sin(phi)}       vs  sin(theta_gamma)

Y-axes standardized to [-0.25, 0.25] with a thin dashed gray y=0 line.

Output:
  output/enpi+/sinthetagamma_harut_cross_check.pdf
"""

import os
import matplotlib.pyplot as plt

# ─────────────────────────────────────────────────────────────────────
# Input data: (x = sin(theta_gamma), y, err)
# ─────────────────────────────────────────────────────────────────────
# Left panel (AUL sinφ)
AULsinphi_h1 = [
    (0.092612812, -0.052871958, 0.043540426),
    (0.130241094, -0.019239326, 0.006427457),
    (0.175812177, -0.005511188, 0.005219611),
    (0.224199380,  0.001744197, 0.004392396),
    (0.271114953,  0.016805510, 0.005184247),
]
AULsinphi_h2 = [
    (0.092183435, -0.128052544, 0.030232566),
    (0.128957361, -0.076749577, 0.011913601),
    (0.176173908, -0.056635366, 0.008055780),
    (0.224994591, -0.019165618, 0.006426446),
    (0.271997590, -0.001764200, 0.007751942),
]
AULsinphi_h3 = [
    (0.092105726,  0.007148937, 0.055761073),
    (0.128731779,  0.006518379, 0.010708090),
    (0.175886625,  0.000708323, 0.009987674),
    (0.224624894,  0.040006724, 0.007578386),
    (0.271768931,  0.066407147, 0.010135559),
]

# Middle panel (AUL sin2φ)
AULsin2phi_h1 = [
    (0.092612812, -0.069683110, 0.071661373),
    (0.130241094, -0.057887185, 0.019352245),
    (0.175812177, -0.023098918, 0.008888072),
    (0.224199380, -0.021983563, 0.007592150),
    (0.271114953, -0.027256545, 0.010393592),
]
AULsin2phi_h2 = [
    (0.092183435,  0.292095704, 0.126174746),
    (0.128957361, -0.105595299, 0.022410160),
    (0.176173908, -0.070311525, 0.017466723),
    (0.224994591, -0.087189004, 0.012293623),
    (0.271997590, -0.081677540, 0.013565507),
]
AULsin2phi_h3 = [
    (0.092105726,  0.021877622, 0.130314600),
    (0.128731779, -0.020769192, 0.021187657),
    (0.175886625, -0.072040283, 0.033027682),
    (0.224624894, -0.084214169, 0.014214842),
    (0.271768931, -0.113232832, 0.016407114),
]

# Right panel (A_tg sinφ)
ATGsinphi_h1 = [
    (0.092612812,  0.049201099, 0.035113873),
    (0.130241094,  0.044869614, 0.019217713),
    (0.175812177,  0.005290574, 0.015876026),
    (0.224199380,  0.001368505, 0.009006246),
    (0.271114953,  0.010644864, 0.008452191),
]
ATGsinphi_h2 = [
    (0.092183435,  0.114796057, 0.022594542),
    (0.128957361, -0.007801625, 0.018673315),
    (0.176173908, -0.023095460, 0.024630955),
    (0.224994591,  0.005236748, 0.013931279),
    (0.271997590,  0.002978469, 0.013140461),
]
ATGsinphi_h3 = [
    (0.092105726, -0.043533569, 0.062025788),
    (0.128731779, -0.017256251, 0.015421493),
    (0.175886625,  0.027435918, 0.023630164),
    (0.224624894, -0.006292719, 0.003934911),
    (0.271768931,  0.001489976, 0.009015100),
]

# ─────────────────────────────────────────────────────────────────────
# Plotting
# ─────────────────────────────────────────────────────────────────────
def unzip(triples):
    x = [t[0] for t in triples]
    y = [t[1] for t in triples]
    e = [t[2] for t in triples]
    return x, y, e

def main():
    # Style knobs
    ylims = (-0.25, 0.25)
    xlabel = r"$\sin\theta_\gamma$"
    left_ylabel   = r"$F_{UL}^{\sin\phi}/F_{UU}$"
    middle_ylabel = r"$F_{UL}^{\sin2\phi}/F_{UU}$"
    right_ylabel  = r"$A_{tg}^{\sin\phi}$"

    # Legend mapping (order matters)
    series = [
        ("0.05 < -t < 0.35", "o", "tab:blue"),
        ("0.35 < -t < 0.85", "s", "tab:orange"),
        ("0.85 < -t < 1.35", "^", "tab:green"),
    ]
    legend_title = r"$0.1 < x_{B} < 0.3$"

    # Gather x-range across all panels for a pleasant shared xlim
    all_x = []
    for block in (AULsinphi_h1, AULsinphi_h2, AULsinphi_h3,
                  AULsin2phi_h1, AULsin2phi_h2, AULsin2phi_h3,
                  ATGsinphi_h1, ATGsinphi_h2, ATGsinphi_h3):
        all_x.extend([t[0] for t in block])
    xmin, xmax = 0, 0.3
    pad = 0.05 * (xmax - xmin if xmax > xmin else 0.05)
    xlims = (xmin - pad, xmax + pad)

    fig, axes = plt.subplots(1, 3, figsize=(14.5, 4.6), sharex=True)
    axL, axM, axR = axes

    # ---- Left: AUL sinφ ----
    for (lab, mkr, clr), data in zip(series, [AULsinphi_h1, AULsinphi_h2, AULsinphi_h3]):
        x, y, e = unzip(data)
        axL.errorbar(x, y, yerr=e, fmt=mkr, color=clr, ecolor=clr,
                     capsize=3, markersize=5.0, linestyle="None", label=lab)
    axL.set(xlim=xlims, ylim=ylims, xlabel=xlabel, ylabel=left_ylabel)
    axL.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    axL.grid(True, linestyle="--", alpha=0.5)
    axL.legend(frameon=True, edgecolor="black", fontsize=10, title=legend_title)

    # ---- Middle: AUL sin2φ ----
    for (lab, mkr, clr), data in zip(series, [AULsin2phi_h1, AULsin2phi_h2, AULsin2phi_h3]):
        x, y, e = unzip(data)
        axM.errorbar(x, y, yerr=e, fmt=mkr, color=clr, ecolor=clr,
                     capsize=3, markersize=5.0, linestyle="None", label=lab)
    axM.set(xlim=xlims, ylim=ylims, xlabel=xlabel, ylabel=middle_ylabel)
    axM.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    axM.grid(True, linestyle="--", alpha=0.5)
    axM.legend(frameon=True, edgecolor="black", fontsize=10, title=legend_title)

    # ---- Right: A_tg sinφ ----
    for (lab, mkr, clr), data in zip(series, [ATGsinphi_h1, ATGsinphi_h2, ATGsinphi_h3]):
        x, y, e = unzip(data)
        axR.errorbar(x, y, yerr=e, fmt=mkr, color=clr, ecolor=clr,
                     capsize=3, markersize=5.0, linestyle="None", label=lab)
    axR.set(xlim=xlims, ylim=ylims, xlabel=xlabel, ylabel=right_ylabel)
    axR.axhline(0, color="gray", linestyle="--", linewidth=0.8)
    axR.grid(True, linestyle="--", alpha=0.5)
    axR.legend(frameon=True, edgecolor="black", fontsize=10, title=legend_title)

    fig.tight_layout()

    out_dir = os.path.join("output", "enpi+")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "sinthetagamma_harut_cross_check.pdf")
    plt.savefig(out_path)
    plt.close(fig)
    print(f"Saved: {out_path}")

if __name__ == "__main__":
    main()