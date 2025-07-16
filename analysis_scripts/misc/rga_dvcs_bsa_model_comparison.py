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
    # path to the ALL18.txt file
    data_file = '/u/home/thayward/clas12_analysis_software/analysis_scripts/misc/import/rga_dvcs_prl/ALL18.txt'

    # columns: phi(rad), Q2, xB, t, Eb, A, sigA
    cols = ['phi_rad', 'Q2', 'xB', 't', 'Eb', 'A', 'sigA']
    df = pd.read_csv(
        data_file,
        delim_whitespace=True,
        comment='#',
        header=None,
        names=cols
    )

    # select the first 19 points
    df = df.iloc[:19]

    # extract arrays
    phi_rad = df['phi_rad'].values
    phi_deg = np.degrees(phi_rad)
    Q2_vals = df['Q2'].values
    xB_vals = df['xB'].values
    t_vals  = df['t'].values
    Eb_vals = df['Eb'].values
    A_vals  = df['A'].values
    sigA    = df['sigA'].values

    # compute KM15 predictions
    model_vals = [
        km15_model(xb, q2, t0, phi, beam_E=Eb)
        for xb, q2, t0, phi, Eb in zip(xB_vals, Q2_vals, t_vals, phi_deg, Eb_vals)
    ]
    model_vals = np.array(model_vals)

    # sort by phi for a smooth line
    sort_idx = np.argsort(phi_deg)

    # plotting
    plt.errorbar(phi_deg, A_vals, yerr=sigA, fmt='o', label='RGA Fa18 Data')
    plt.plot(phi_deg[sort_idx], model_vals[sort_idx], '-', label='KM15 Model')

    plt.xlabel(r'$\phi$ (deg)')
    plt.ylabel(r'$A_{LU}$')
    plt.legend(loc='upper right')

    # save output
    out_dir = 'output/rga_dvcs_BSA_model_comparison'
    os.makedirs(out_dir, exist_ok=True)
    plt.tight_layout()
    plt.savefig(os.path.join(out_dir, 'first_test.pdf'))
    plt.close()

if __name__ == '__main__':
    main()