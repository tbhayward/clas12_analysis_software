#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import matplotlib.pyplot as plt

def main():
    # --------------------------
    # Raw data: each row = [t (negative), F_LU^{sin phi}, err]
    # We'll plot versus -t, so we negate the first column.
    # --------------------------
    data_rga = [
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

    data_rgc = [
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
    # Parse into arrays
    # --------------------------
    x_rga, y_rga, e_rga = [], [], []
    for tneg, y, e in data_rga:
        x_rga.append(-tneg)  # -t
        y_rga.append(y)
        e_rga.append(e)
    #endfor

    x_rgc, y_rgc, e_rgc = [], [], []
    for tneg, y, e in data_rgc:
        x_rgc.append(-tneg)  # -t
        y_rgc.append(y)
        e_rgc.append(e)
    #endfor

    x_rga = np.asarray(x_rga, dtype=float)
    y_rga = np.asarray(y_rga, dtype=float)
    e_rga = np.asarray(e_rga, dtype=float)

    x_rgc = np.asarray(x_rgc, dtype=float)
    y_rgc = np.asarray(y_rgc, dtype=float)
    e_rgc = np.asarray(e_rgc, dtype=float)

    # --------------------------
    # Plot
    # --------------------------
    fig, ax = plt.subplots(figsize=(6.5, 4.5))

    ax.errorbar(
        x_rga, y_rga, yerr=e_rga,
        fmt='o', ms=4, lw=1.0, capsize=2,
        color='blue', label=r'RGA Fa18 Inb'
    )
    ax.errorbar(
        x_rgc, y_rgc, yerr=e_rgc,
        fmt='s', ms=4, lw=1.0, capsize=2,
        color='red', label=r'RGC Fa22 NH$_{3}$ Inb'
    )

    ax.set_xlabel(r'$-t$ (GeV$^{2}$)')
    ax.set_ylabel(r'$F_{LU}^{\sin\phi}$')
    ax.set_xlim(0.0, 1.3)
    ax.set_ylim(-0.1, 0.2)
    ax.legend(frameon=False)
    ax.grid(True, alpha=0.3, linestyle='--')

    # --------------------------
    # Save
    # --------------------------
    outdir = os.path.join('output', 'enpi+')
    if not os.path.isdir(outdir):
        os.makedirs(outdir, exist_ok=True)
    #endif

    outpath = os.path.join(outdir, 'rga_rgc_comparison.pdf')
    plt.tight_layout()
    fig.savefig(outpath, bbox_inches='tight')
    print('Saved:', outpath)
#endfor

if __name__ == '__main__':
    main()
#endif