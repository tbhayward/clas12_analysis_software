# root_io.py

import os
import ROOT

# Configure ROOT to run in batch mode
ROOT.gROOT.SetBatch(True)

def load_root_files(period):
    """
    Load ROOT files with only required branches.
    Returns (period, {'data': TChain, 'mc': TChain}).
    """
    file_map = {
        "DVCS_Fa18_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root"
        },
        "DVCS_Fa18_out": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root"
        },
        "DVCS_Sp19_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
            "mc":   "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"
        }
    }

    trees = {}
    for category in ["data", "mc"]:
        chain = ROOT.TChain("PhysicsEvents")
        chain.Add(file_map[period][category])

        # Activate only the required branches:
        branches = [
            "detector1", "detector2", "t1", "open_angle_ep2", "theta_gamma_gamma",
            "Emiss2", "Mx2", "Mx2_1", "Mx2_2", "pTmiss", "xF"
        ]
        for br in branches:
            chain.SetBranchStatus(br, 1)
        #endfor

        trees[category] = chain
    #endfor

    return period, trees
#enddef