#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import gepard as g
from gepard.fits import th_KM15

def km15_model(xB, Q2, t_pos, phi_deg, beam_E=10.604):
    """
    Compute the BKM-KM15 prediction for ALU at a single data point.
    """
    t_km15     = -abs(t_pos)
    phi_rad    = np.radians(phi_deg)
    phi_trento = np.pi - phi_rad

    pt = g.DataPoint(
        xB             = xB,
        t              = t_km15,
        Q2             = Q2,
        phi            = phi_trento,
        observable     = 'ALU',
        frame          = 'trento',
        process        = 'ep2epgamma',
        exptype        = 'fixed target',
        in1energy      = beam_E,
        in1charge      = -1,
        in1polarization= +1,
        in2polarization= -1,
    )
    pt.prepare()
    return th_KM15.predict(pt)

def main():
    data_file = '/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/import/rga_dvcs_prl/ALL18.txt'
    cols = ['phi_rad', 'Q2', 'xB', 't', 'Eb', 'A', 'sigA']
    df = pd.read_csv(data_file, sep=r'\s+', comment='#', header=None, names=cols)

    # detect bin boundaries where phi resets
    phi = df['phi_rad'].to_numpy()
    reset_idxs = np.where(phi[1:] < phi[:-1])[0] + 1
    idx_groups = np.split(np.arange(len(df)), reset_idxs)

    # collect info per (Q2, xB, t) segment
    seg_info = []
    for idxs in idx_groups:
        seg = df.iloc[idxs]
        mean_xB = seg['xB'].mean()
        mean_Q2 = seg['Q2'].mean()
        mean_t  = seg['t'].mean()
        mean_Eb = seg['Eb'].mean()

        phi_deg = np.degrees(seg['phi_rad'].to_numpy())
        A_vals  = seg['A'].to_numpy()
        sigA    = seg['sigA'].to_numpy()

        # model on a fine phi grid
        phi_grid = np.linspace(0, 360, 100)
        model_vals = np.array([
            km15_model(mean_xB, mean_Q2, mean_t, phi, beam_E=mean_Eb)
            for phi in phi_grid
        ])

        seg_info.append({
            'mean_xB':  mean_xB,
            'mean_Q2':  mean_Q2,
            'mean_t':   mean_t,
            'phi_deg':  phi_deg,
            'A_vals':   A_vals,
            'sigA':     sigA,
            'phi_grid': phi_grid,
            'model':    model_vals
        })

    out_dir = 'output/rga_dvcs_BSA_model_comparison'
    os.makedirs(out_dir, exist_ok=True)

    # group by t in 0.1 steps
    by_tbin = {}
    for si in seg_info:
        tbin = round(si['mean_t'], 1)  # 0.1 GeV^2 bins
        by_tbin.setdefault(tbin, []).append(si)

    for tbin, group in sorted(by_tbin.items()):
        Q2_vals = sorted({round(si['mean_Q2'],6) for si in group})
        xB_vals = sorted({round(si['mean_xB'],6) for si in group})

        nrows, ncols = len(Q2_vals), len(xB_vals)
        fig, axes = plt.subplots(
            nrows, ncols,
            figsize=(3*ncols, 3*nrows),
            sharex=True, sharey=True
        )
        axes = np.atleast_2d(axes)

        first_legend = True
        for si in group:
            i = Q2_vals.index(round(si['mean_Q2'],6))
            j = xB_vals.index(round(si['mean_xB'],6))
            ax = axes[nrows-1-i, j]

            h_data = ax.errorbar(si['phi_deg'], si['A_vals'], yerr=si['sigA'], fmt='o')
            h_model, = ax.plot(si['phi_grid'], si['model'], '-')

            if first_legend:
                ax.legend([h_data, h_model], ['RGA Fa18 Data', 'KM15 Model'], loc='upper right')
                first_legend = False

            if j == 0:
                ax.set_ylabel(rf"$Q^2={si['mean_Q2']:.2f}\,\mathrm{{GeV}}^2$")
            if nrows-1-i == nrows-1:
                ax.set_xlabel(r'$\phi$ (deg)')

            ax.set_xlim(0, 360)
            ax.set_ylim(-0.6, 0.6)

            if nrows-1-i == nrows-1:
                ax.set_title(rf"$x_B={si['mean_xB']:.3f}$")

        fig.suptitle(rf"$t\approx{tbin:.1f}\,\mathrm{{GeV}}^2$", fontsize=10)
        plt.tight_layout(rect=[0,0,1,0.95])
        fig.savefig(os.path.join(out_dir, f"tbin_{tbin:.1f}.pdf"))
        plt.close(fig)

if __name__ == '__main__':
    main()