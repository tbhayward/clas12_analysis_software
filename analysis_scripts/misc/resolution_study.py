#!/usr/bin/env python3

import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
from concurrent.futures import ProcessPoolExecutor

# Physical constants (GeV units)
M_PROTON = 0.938272
M_ELECTRON = 0.000511
M_PHOTON = 0.0

# Beam energies by run
BEAM_ENERGIES = {
    'rga_fa18_inb': 10.6041,
    'rga_fa18_out': 10.6041,
    'rga_sp19_inb': 10.1998
}

# List of runs with full paths
RUNS = [
    {
        'name': 'rga_fa18_inb',
        'title': 'RGA Fa18 Inb',
        'mc_file':  '/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/'
                    'dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root',
        'data_file':'/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/'
                    'dvcs/rga_fa18_inb_epgamma.root'
    },
    {
        'name': 'rga_fa18_out',
        'title': 'RGA Fa18 Out',
        'mc_file':  '/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/'
                    'dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root',
        'data_file':'/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/'
                    'dvcs/rga_fa18_out_epgamma.root'
    },
    {
        'name': 'rga_sp19_inb',
        'title': 'RGA Sp19 Inb',
        'mc_file':  '/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/'
                    'dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root',
        'data_file':'/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/'
                    'dvcs/rga_sp19_inb_epgamma.root'
    }
]

# Branches to compute: missing mass types
BRANCH_SETTINGS = [
    ('Mx2_1', (-0.5, 0.5), r'$M_{x(p)}^2$ (GeV$^2$)', 'proton'),
    ('Mx2_2', (0.4,  1.6), r'$M_{x(\gamma)}^2$ (GeV$^2$)', 'photon')
]

# Detector topology definitions
TOPOLOGIES = [
    {'det1':1,'det2':0,'label':'FD–FT'},
    {'det1':1,'det2':1,'label':'FD–FD'},
    {'det1':2,'det2':0,'label':'CD–FT'},
    {'det1':2,'det2':1,'label':'CD–FD'}
]

# Resolution functions (corrected per Eqns 1.1-1.2)
def p_resolution(p, th):
    pS1 = 0.018429 - 0.011008*th + 0.0022766*th**2 - 0.00014015*th**3 + 3.07424e-6*th**4
    return 0.02 * np.sqrt((pS1 * p)**2 + (0.02 * th)**2)

def theta_resolution(p, th):
    # Polar-angle resolution θ_R
    num = (0.004*th + 0.1)**2 * p**2 + 0.139572**2
    return 2.5 * np.sqrt(num) / p

def phi_resolution(p, th):
    # Azimuthal-angle resolution φ_R
    phiS1 = 0.85 - 0.015*th
    phiS2 = 0.17 - 0.003*th
    return 3.5 * np.sqrt((phiS1 * np.sqrt(p**2 + 0.139572**2) / p)**2 + phiS2**2)

# Four-vector builder for missing mass calculation
def four_vector(p, theta, phi, mass):
    px = p * np.sin(theta) * np.cos(phi)
    py = p * np.sin(theta) * np.sin(phi)
    pz = p * np.cos(theta)
    E  = np.sqrt(p**2 + mass**2)
    return np.vstack([E, px, py, pz])

# Compute missing mass squared after smearing
def compute_missing_mass2(e_p, e_th, e_ph,
                          pr_p, pr_th, pr_ph,
                          ph_p, ph_th, ph_ph,
                          beam_e):
    p4_beam = np.array([beam_e, 0.0, 0.0, beam_e])[:,None]
    p4_tgt  = np.array([M_PROTON, 0.0, 0.0, 0.0])[:,None]
    p4_init = p4_beam + p4_tgt

    p4_e  = four_vector(e_p,  e_th,  e_ph,  M_ELECTRON)
    p4_pr = four_vector(pr_p, pr_th, pr_ph, M_PROTON)
    p4_ph = four_vector(ph_p, ph_th, ph_ph, M_PHOTON)

    p4_m1  = p4_init - p4_e - p4_pr
    mx2_1  = p4_m1[0]**2 - np.sum(p4_m1[1:]**2, axis=0)
    p4_m2  = p4_init - p4_e - p4_ph
    mx2_2  = p4_m2[0]**2 - np.sum(p4_m2[1:]**2, axis=0)
    return mx2_1.ravel(), mx2_2.ravel()

# Process one run: apply cuts, smear, compute histograms
def process_run(run):
    tree_mc = uproot.open(run['mc_file'])['PhysicsEvents']
    tree_dt = uproot.open(run['data_file'])['PhysicsEvents']
    beam_e  = BEAM_ENERGIES[run['name']]
    # Read arrays
    keys = ['e_p','e_theta','e_phi', 'p1_p','p1_theta','p1_phi', 'p2_p','p2_theta','p2_phi',
            'detector1','detector2','t1','theta_gamma_gamma','pTmiss','Emiss2','Mx2_1','Mx2_2']
    arr_mc = {k: tree_mc[k].array(library='np') for k in keys}
    arr_dt = {k: tree_dt[k].array(library='np') for k in keys}
    # Kinematic cuts
    mask_dt = ((np.abs(arr_dt['t1'])<1) & (arr_dt['theta_gamma_gamma']<0.4) &
               (arr_dt['pTmiss']<0.05) & (arr_dt['Emiss2']<1))
    mask_mc = ((np.abs(arr_mc['t1'])<1) & (arr_mc['theta_gamma_gamma']<0.4) &
               (arr_mc['pTmiss']<0.05) & (arr_mc['Emiss2']<1))
    # Smear MC (AS=1)
    rng     = np.random.randn
    e_p_sm  = arr_mc['e_p']     + p_resolution(arr_mc['e_p'],arr_mc['e_theta']) * rng(len(arr_mc['e_p']))
    e_th_sm = arr_mc['e_theta'] + theta_resolution(arr_mc['e_p'],arr_mc['e_theta']) * rng(len(arr_mc['e_theta']))
    e_ph_sm = arr_mc['e_phi']   + phi_resolution(arr_mc['e_p'],arr_mc['e_theta'])   * rng(len(arr_mc['e_phi']))
    pr_p_sm  = arr_mc['p1_p']    + p_resolution(arr_mc['p1_p'],arr_mc['p1_theta'])* rng(len(arr_mc['p1_p']))
    pr_th_sm = arr_mc['p1_theta']+ theta_resolution(arr_mc['p1_p'],arr_mc['p1_theta'])* rng(len(arr_mc['p1_theta']))
    pr_ph_sm = arr_mc['p1_phi']  + phi_resolution(arr_mc['p1_p'],arr_mc['p1_theta'])* rng(len(arr_mc['p1_phi']))
    ph_p_sm  = arr_mc['p2_p']    + p_resolution(arr_mc['p2_p'],arr_mc['p2_theta'])* rng(len(arr_mc['p2_p']))
    ph_th_sm = arr_mc['p2_theta']+ theta_resolution(arr_mc['p2_p'],arr_mc['p2_theta'])* rng(len(arr_mc['p2_theta']))
    ph_ph_sm = arr_mc['p2_phi']  + phi_resolution(arr_mc['p2_p'],arr_mc['p2_theta'])* rng(len(arr_mc['p2_phi']))
    mx2_1_sm, mx2_2_sm = compute_missing_mass2(e_p_sm,e_th_sm,e_ph_sm,
                                               pr_p_sm,pr_th_sm,pr_ph_sm,
                                               ph_p_sm,ph_th_sm,ph_ph_sm,
                                               beam_e)
    # Build histograms and stats
    results = {b:{} for b,*_ in BRANCH_SETTINGS}
    for branch, xlim, xlabel, _ in BRANCH_SETTINGS:
        bins = np.linspace(xlim[0],xlim[1],101)
        for topo in TOPOLOGIES:
            mask_td = (mask_dt & (arr_dt['detector1']==topo['det1']) & (arr_dt['detector2']==topo['det2']))
            mask_tm = (mask_mc & (arr_mc['detector1']==topo['det1']) & (arr_mc['detector2']==topo['det2']))
            dt_vals   = arr_dt[branch][mask_td]
            mc_before = arr_mc[branch][mask_tm]
            mc_after  = mx2_1_sm[mask_tm] if branch=='Mx2_1' else mx2_2_sm[mask_tm]
            dt_c,_ = np.histogram(dt_vals,   bins=bins, density=True)
            mb_c,_ = np.histogram(mc_before, bins=bins, density=True)
            ma_c,_ = np.histogram(mc_after,  bins=bins, density=True)
            dv      = dt_vals[(dt_vals>=xlim[0])&(dt_vals<=xlim[1])]
            bv      = mc_before[(mc_before>=xlim[0])&(mc_before<=xlim[1])]
            av      = mc_after[(mc_after>=xlim[0])&(mc_after<=xlim[1])]
            results[branch][topo['label']] = {
                'bins':bins, 'dt':(dt_c, dv.mean(), dv.std()),
                'mc_before':(mb_c, bv.mean(), bv.std()),
                'mc_after':(ma_c, av.mean(), av.std())
            }
    return run['name'], results

# Plotting
def plot_results(all_results):
    for branch, xlim, xlabel, _ in BRANCH_SETTINGS:
        fig,axes = plt.subplots(len(RUNS),len(TOPOLOGIES),
                               figsize=(4*len(TOPOLOGIES),3*len(RUNS)),
                               sharex=True,sharey=True)
        for i,run in enumerate(RUNS):
            res = all_results[run['name']][branch]
            for j,topo in enumerate(TOPOLOGIES):
                ax = axes[i,j]
                r  = res[topo['label']]
                centers = 0.5*(r['bins'][:-1]+r['bins'][1:])
                dt_c,mu_dt,sig_dt = r['dt']
                mb_c,mu_mb,sig_mb = r['mc_before']
                ma_c,mu_ma,sig_ma = r['mc_after']
                ax.step(centers,dt_c,where='mid', label=f"Data (μ={mu_dt:.3f},σ={sig_dt:.3f})")
                ax.step(centers,mb_c,where='mid',linestyle='--',label=f"MC before (μ={mu_mb:.3f},σ={sig_mb:.3f})")
                ax.step(centers,ma_c,where='mid',linestyle='-.',label=f"MC after  (μ={mu_ma:.3f},σ={sig_ma:.3f})")
                ax.set_xlim(xlim); ax.set_title(f"{run['title']}\n{topo['label']}",fontsize=10)
                ax.legend(loc='upper right',fontsize=7)
                if i==len(RUNS)-1: ax.set_xlabel(xlabel)
                if j==0:          ax.set_ylabel('normalized counts')
        fig.tight_layout()
        out = f"output/resolution_study/{branch}_smearing_comparison.pdf"
        os.makedirs(os.path.dirname(out),exist_ok=True)
        fig.savefig(out); plt.close(fig)

# Main
def main():
    with ProcessPoolExecutor() as execu:
        futures = [execu.submit(process_run,run) for run in RUNS]
        all_res = {name:res for name,res in (f.result() for f in futures)}
    plot_results(all_res)

if __name__=='__main__':
    main()
#endif