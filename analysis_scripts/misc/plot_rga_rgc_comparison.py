#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt

def parse_to_arrays(rows):
    """Convert [[tneg, y, e], ...] -> x=-t, y, e numpy arrays."""
    x, y, e = [], [], []
    for tneg, yy, ee in rows:
        x.append(-float(tneg))   # -t
        y.append(float(yy))
        e.append(float(ee))
    #endfor
    return np.asarray(x), np.asarray(y), np.asarray(e)
#endfor

def main():
    # --------------------------
    # LEFT PANEL: x_B in [0.10, 0.60]
    # --------------------------
    data_rga_all = [
        [-1.199099929, 0.106094505, 0.014567704],
        [-1.098606216, 0.150010511, 0.013594042],
        [-0.998715748, 0.126960020, 0.012689497],
        [-0.898456194, 0.142971560, 0.011798005],
        [-0.798613367, 0.139964185, 0.010980846],
        [-0.698573609, 0.101431440, 0.010246570],
        [-0.598691251, 0.126407286, 0.009617647],
        [-0.498604049, 0.107684259, 0.008984012],
        [-0.398688895, 0.128108519, 0.008583905],
        [-0.298635247, 0.125430412, 0.008167959],
        [-0.199615240, 0.098419617, 0.007572824],
        [-0.106069254, 0.066021980, 0.007586575],
    ]
    data_rgc_all = [
        [-1.199254858, 0.115608756, 0.026206482],
        [-1.098419750, 0.065197368, 0.024542325],
        [-0.998794627, 0.092106925, 0.023140667],
        [-0.898830973, 0.156998810, 0.021486700],
        [-0.798462126, 0.137654081, 0.020014104],
        [-0.698657933, 0.134173438, 0.018416641],
        [-0.598504529, 0.098279741, 0.017362193],
        [-0.498563527, 0.146441720, 0.016207771],
        [-0.398570196, 0.147899288, 0.015109343],
        [-0.298828009, 0.135025430, 0.014143565],
        [-0.200246087, 0.097417894, 0.013434122],
        [-0.109031497, 0.058112801, 0.014815651],
    ]

    # --------------------------
    # MIDDLE PANEL: x_B in [0.10, 0.25]
    # --------------------------
    data_rga_low = [
        [-1.145084992, 0.054622672, 0.023743843],
        [-0.945084171, 0.051959110, 0.020716086],
        [-0.745308617, 0.071202111, 0.018329393],
        [-0.542407324, 0.068584436, 0.015417029],
        [-0.338122978, 0.103860346, 0.012533450],
        [-0.127154681, 0.083377226, 0.007795620],
    ]
    data_rgc_low = [
        [-1.144558703, 0.034687712, 0.043386411],
        [-0.944756340, -0.015494051, 0.038461909],
        [-0.742782191, 0.034221919, 0.032816930],
        [-0.542254429, 0.083657502, 0.027904226],
        [-0.338530753, 0.150299892, 0.021888257],
        [-0.133457081, 0.059036585, 0.014071653],
    ]

    # --------------------------
    # RIGHT PANEL: x_B in [0.45, 0.60]
    # --------------------------
    data_rga_high = [
        [-1.145403261, 0.140435559, 0.018536696],
        [-0.944694169, 0.156121547, 0.016271461],
        [-0.745924501, 0.139093691, 0.014378928],
        [-0.550344018, 0.125560504, 0.014023515],
        [-0.386017760, 0.058567029, 0.021135413],
        [-0.245797394, -0.098866164, 0.709697519],
    ]
    data_rgc_high = [
        [-1.145182109, 0.100559171, 0.032409370],
        [-0.945414028, 0.153935709, 0.028465377],
        [-0.744747106, 0.171656072, 0.025044578],
        [-0.549817284, 0.120362491, 0.024436794],
        [-0.388625205, 0.071432401, 0.038096994],
        [-0.248053003, -8.700734592, 0.058707248],
    ]

    # Convert to arrays
    xA, yA, eA = parse_to_arrays(data_rga_all)
    xB, yB, eB = parse_to_arrays(data_rgc_all)

    xL_A, yL_A, eL_A = parse_to_arrays(data_rga_low)
    xL_B, yL_B, eL_B = parse_to_arrays(data_rgc_low)

    xH_A, yH_A, eH_A = parse_to_arrays(data_rga_high)
    xH_B, yH_B, eH_B = parse_to_arrays(data_rgc_high)

    # Figure and axes
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.5), sharey=True)
    (ax1, ax2, ax3) = axes

    # Common styling
    xlim = (0.0, 1.3)
    ylim = (-0.1, 0.2)

    # Left panel
    l1 = ax1.errorbar(
        xA, yA, yerr=eA, fmt='o', ms=4, lw=1.0, capsize=2,
        color='blue', label=r'RGA Fa18 Inb'
    )
    l2 = ax1.errorbar(
        xB, yB, yerr=eB, fmt='s', ms=4, lw=1.0, capsize=2,
        color='red', label=r'RGC Fa22 NH$_{3}$ Inb'
    )
    ax1.set_title(r'$x_{B}\ \in\ [0.10,\,0.60]$')
    ax1.set_xlabel(r'$-t$ (GeV$^{2}$)')
    ax1.set_ylabel(r'$F_{LU}^{\sin\phi}$')
    ax1.set_xlim(*xlim)
    ax1.set_ylim(*ylim)
    ax1.grid(True, alpha=0.3, linestyle='--')

    # Middle panel
    ax2.errorbar(
        xL_A, yL_A, yerr=eL_A, fmt='o', ms=4, lw=1.0, capsize=2, color='blue'
    )
    ax2.errorbar(
        xL_B, yL_B, yerr=eL_B, fmt='s', ms=4, lw=1.0, capsize=2, color='red'
    )
    ax2.set_title(r'$x_{B}\ \in\ [0.10,\,0.25]$')
    ax2.set_xlabel(r'$-t$ (GeV$^{2}$)')
    ax2.set_xlim(*xlim)
    ax2.set_ylim(*ylim)
    ax2.grid(True, alpha=0.3, linestyle='--')

    # Right panel
    ax3.errorbar(
        xH_A, yH_A, yerr=eH_A, fmt='o', ms=4, lw=1.0, capsize=2, color='blue'
    )
    ax3.errorbar(
        xH_B, yH_B, yerr=eH_B, fmt='s', ms=4, lw=1.0, capsize=2, color='red'
    )
    ax3.set_title(r'$x_{B}\ \in\ [0.45,\,0.60]$')
    ax3.set_xlabel(r'$-t$ (GeV$^{2}$)')
    ax3.set_xlim(*xlim)
    ax3.set_ylim(*ylim)
    ax3.grid(True, alpha=0.3, linestyle='--')

    # One shared legend for the whole figure
    fig.legend(
        handles=[l1, l2],
        labels=[r'RGA Fa18 Inb', r'RGC Fa22 NH$_{3}$ Inb'],
        loc='upper center', ncol=2, frameon=False, bbox_to_anchor=(0.5, 1.03)
    )

    plt.tight_layout()

    # Save
    outdir = os.path.join('output', 'enpi+')
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, 'rga_rgc_comparison_1x3.pdf')
    plt.savefig(outpath, bbox_inches='tight')
    print('Saved:', outpath)
#endif

if __name__ == '__main__':
    main()
#endif