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
    # Read first 29 points from ALL18.txt
    data_file = '/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/import/rga_dvcs_prl/ALL18.txt'
    cols = ['phi_rad', 'Q2', 'xB', 't', 'Eb', 'A', 'sigA']
    df = pd.read_csv(
        data_file,
        sep=r'\s+',
        comment='#',
        header=None,
        names=cols
    ).iloc[:29]  # 29 points

    # Compute bin-averaged kinematics
    mean_xB = df['xB'].mean()
    mean_Q2 = df['Q2'].mean()
    mean_t  = df['t'].mean()
    mean_Eb = df['Eb'].mean()

    # Data for plotting
    phi_deg = np.degrees(df['phi_rad'].values)
    A_vals  = df['A'].values
    sigA    = df['sigA'].values

    # Create a phi grid (100 points) from 0–360°
    phi_grid = np.linspace(0, 360, 100)

    # Compute KM15 model on that grid using mean kinematics
    model_vals = np.array([
        km15_model(mean_xB, mean_Q2, mean_t, phi, beam_E=mean_Eb)
        for phi in phi_grid
    ])

    # Plot data points with error bars
    plt.errorbar(phi_deg, A_vals, yerr=sigA, fmt='o', label='RGA Fa18 Data')

    # Plot KM15 model curve
    plt.plot(phi_grid, model_vals, '-', label='KM15 Model')

    # Labels, limits, legend, title
    plt.xlabel(r'$\phi$ (deg)')
    plt.ylabel(r'$A_{LU}$')
    plt.xlim(0, 360)
    plt.ylim(-0.6, 0.6)
    plt.legend(loc='upper right')
    plt.title(
        rf"$\langle x_B\rangle={mean_xB:.3f},\ \langle Q^2\rangle={mean_Q2:.3f}\,\mathrm{{GeV}}^2,\ "
        rf"\langle t\rangle={mean_t:.3f}\,\mathrm{{GeV}}^2,\ E_b={mean_Eb:.3f}\,\mathrm{{GeV}}$",
        fontsize=10
    )

    # Save figure
    out_dir = 'output/rga_dvcs_BSA_model_comparison'
    os.makedirs(out_dir, exist_ok=True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'first_test.pdf'))
    plt.close()

if __name__ == '__main__':
    main()