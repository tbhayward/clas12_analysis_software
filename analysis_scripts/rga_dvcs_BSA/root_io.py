# root_io.py

import os
import ROOT

# Configure ROOT to run in batch mode so we don't pop up windows
ROOT.gROOT.SetBatch(True)

def load_root_files(period):
    """
    Load ROOT files with only required branches. We handle DVCS or eppi0
    by checking if period starts with 'eppi0' or 'DVCS'.
    Returns (period, {'data': TChain, 'mc': TChain}).
    """
    # Mapping from period -> file paths
    # We'll prefix them with "DVCS_" or "eppi0_" to differentiate
    file_map = {
        "DVCS_Fa18_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_10604MeV.root"
        },
        "DVCS_Fa18_out": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_10604MeV.root"
        },
        "DVCS_Sp19_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_10200MeV.root"
        },
        "DVCS_Sp18_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_inb_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_inb_10594MeV.root"
        },
        "DVCS_Sp18_out": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp18_out_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp18_out_10594MeV.root"
        },
        "eppi0_Fa18_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_inb_eppi0.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_inb_10604MeV.root"
        },
        "eppi0_Fa18_out": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_fa18_out_eppi0.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_fa18_out_10604MeV.root"
        },
        "eppi0_Sp19_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp19_inb_eppi0.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp19_inb_10200MeV.root"
        },
        "eppi0_Sp18_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_inb_eppi0.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_inb_10594MeV.root"
        },
        "eppi0_Sp18_out": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/data/rga_sp18_out_eppi0.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_sp18_out_10594MeV.root"
        },
        "eppi0_bkg_Fa18_inb": {
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_epgamma.root"
        },
        "eppi0_bkg_Fa18_out": {
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_epgamma.root"
        },
        "eppi0_bkg_Sp19_inb": {
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_epgamma.root"
        },
        "eppi0_bkg_Sp18_inb": {
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_inb_epgamma.root"
        },
        "eppi0_bkg_Sp19_out": {
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp18_out_epgamma.root"
        }
    }

    if period not in file_map:
        raise ValueError(f"Unrecognized period: {period}. Check your file_map keys.")

    trees = {}
    for category in file_map[period]:  # Only iterate over available categories
        chain = ROOT.TChain("PhysicsEvents")
        chain.Add(file_map[period][category])

        # If period starts with 'eppi0', we assume we need 'theta_pi0_pi0' branch
        # Otherwise we assume DVCS with 'theta_gamma_gamma'
        if period.startswith("eppi0_bkg") or period.startswith("DVCS"):
            # Background files, even though they start with eppi0, are DVCS-like.
            branches = [
                "beam_pol", "helicity", "x", "Q2", "detector1", "detector2", "t1", 
                "open_angle_ep2", "theta_gamma_gamma",
                "Emiss2", "Mx2", "Mx2_1", "Mx2_2", "pTmiss", "xF", "phi2"
            ]
        elif period.startswith("eppi0"):
            # Regular eppi0 files.
            branches = [
                "beam_pol", "helicity", "x", "Q2", "detector1", "detector2", "t1", 
                "open_angle_ep2", "theta_pi0_pi0",
                "Emiss2", "Mx2", "Mx2_1", "Mx2_2", "pTmiss", "xF", "phi2"
            ]

        # Activate only these required branches:
        for br in branches:
            chain.SetBranchStatus(br, 1)
        #endfor

        trees[category] = chain
    #endfor

    return period, trees
#enddef