#!/usr/bin/env python3

import os
import json
import math
import ROOT

# Ensure ROOT runs in batch mode
ROOT.gROOT.SetBatch(True)

# --------------------------------------------------------------------------------------
# 1) apply_kinematic_cuts
# --------------------------------------------------------------------------------------
def apply_kinematic_cuts(t_value, open_angle_ep2_value, theta_neutral_neutral_value,
                         Emiss2_value, Mx2_value, Mx2_1_value, Mx2_2_value, pTmiss_value,
                         xF_value, analysis_type, data_type, run_period, topology):
    return True  # Replace with actual cuts

# --------------------------------------------------------------------------------------
# 2) format_label_name
# --------------------------------------------------------------------------------------
def format_label_name(variable, analysis_type):
    labels = {
        "open_angle_ep2": "#theta_{e'#gamma} (deg)",
        "Mx2_2": "Miss M_{e'#gamma}^{2} (GeV^{2})",
        "theta_gamma_gamma": "#theta_{#gamma#gamma}",
        "theta_pi0_pi0": "#theta_{#pi^{0}#pi^{0}}",
        "xF": "x_{F}",
        "Emiss2": "Miss E^{2} (GeV^{2})",
        "Mx2": "Miss M^{2} (GeV^{2})",
        "Mx2_1": "Miss M_{e'p'}^{2} (GeV^{2})",
        "pTmiss": "p_{T}^{miss} (GeV)"
    }
    label = labels.get(variable, variable)
    if analysis_type == "eppi0":
        label = label.replace("#gamma", "#pi^{0}")
    return label

# --------------------------------------------------------------------------------------
# 3) load_root_files
# --------------------------------------------------------------------------------------
def load_root_files():
    """
    Loads your ROOT files into TChains keyed by (category, run period).
    """
    file_dict = {
        # dvcs data:
        "dvcs_data": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"
        }, 
        # dvcs mc (reconstructed):
        "dvcs_mc": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"
        }
        # eppi0-related files omitted, since we're not using them yet
    }

    tree_dict = {}
    for category, subdict in file_dict.items():
        tree_dict[category] = {}
        for period, filename in subdict.items():
            chain = ROOT.TChain("PhysicsEvents")
            chain.Add(filename)
            tree_dict[category][period] = chain
        #endfor
    #endfor

    return tree_dict

# --------------------------------------------------------------------------------------
# 4) determine_exclusivity_cuts (UPDATED)
# --------------------------------------------------------------------------------------
def determine_exclusivity_cuts(analysis_type, topology, data_tree, mc_tree,
                               output_dir, run_period, plot_title):
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.05, "XY")
    ROOT.gStyle.SetLabelSize(0.04, "XY")

    # Create output directory
    analysis_dir = "dvcs" if analysis_type == "dvcs" else "eppi0"
    final_output_dir = os.path.join(output_dir, analysis_dir)
    os.makedirs(final_output_dir, exist_ok=True)

    # Initialize variables and histograms
    variables = ["open_angle_ep2", "theta_gamma_gamma", "pTmiss", "xF",
                 "Emiss2", "Mx2", "Mx2_1", "Mx2_2"] if analysis_type == "dvcs" else []

    hist_configs = {
        "open_angle_ep2": (100, 0, 40),
        "theta_gamma_gamma": (100, 0, 40),
        "pTmiss": (100, 0, 1.5),
        "xF": (100, -1, 1),
        "Emiss2": (100, -0.5, 1),
        "Mx2": (100, -0.5, 1.5),
        "Mx2_1": (100, -0.5, 1.5),
        "Mx2_2": (100, -0.5, 1.5)
    }

    # Initialize all histograms
    hists = {
        'data': {var: ROOT.TH1D(f"data_{var}", "", *hist_configs[var]) for var in variables},
        'mc': {var: ROOT.TH1D(f"mc_{var}", "", *hist_configs[var]) for var in variables},
        'data_loose': {var: ROOT.TH1D(f"data_loose_{var}", "", *hist_configs[var]) for var in variables},
        'mc_loose': {var: ROOT.TH1D(f"mc_loose_{var}", "", *hist_configs[var]) for var in variables}
    }

    # Style histograms
    for key in ['data', 'data_loose']:
        for hist in hists[key].values():
            hist.SetMarkerStyle(20)
            hist.SetMarkerColor(ROOT.kBlue)
            hist.SetLineColor(ROOT.kBlue)
    for key in ['mc', 'mc_loose']:
        for hist in hists[key].values():
            hist.SetMarkerStyle(21)
            hist.SetMarkerColor(ROOT.kRed)
            hist.SetLineColor(ROOT.kRed)

    # Process data tree
    for event in data_tree:
        det1 = event.detector1
        det2 = event.detector2
        if not ((topology == "(FD,FD)" and det1 == 1 and det2 == 1) or
                (topology == "(CD,FD)" and det1 == 2 and det2 == 1) or
                (topology == "(CD,FT)" and det1 == 2 and det2 == 0)):
            continue

        # Get cut variables
        t_val = event.t1
        open_ang = event.open_angle_ep2
        theta_nn = event.theta_gamma_gamma if analysis_type == "dvcs" else event.theta_pi0_pi0
        cuts_passed = apply_kinematic_cuts(
            t_val, open_ang, theta_nn, event.Emiss2, event.Mx2,
            event.Mx2_1, event.Mx2_2, event.pTmiss, event.xF,
            analysis_type, "data", run_period, topology
        )

        # Fill histograms
        for var in variables:
            val = getattr(event, var)
            hists['data'][var].Fill(val)
            if cuts_passed:
                hists['data_loose'][var].Fill(val)

    # Process MC tree (similar structure)
    # ... [Repeat same structure for MC tree] ...

    # Normalize and draw
    canvas = ROOT.TCanvas("canvas", "Exclusivity Plots", 1600, 800)
    canvas.Divide(4, 2)
    
    for i, var in enumerate(variables):
        canvas.cd(i+1)
        hists['data'][var].Scale(1/hists['data'][var].Integral() if hists['data'][var].Integral() > 0 else 1)
        hists['mc'][var].Scale(1/hists['mc'][var].Integral() if hists['mc'][var].Integral() > 0 else 1)
        hists['data'][var].Draw("E1")
        hists['mc'][var].Draw("E1 SAME")
        # Add legends/labels

    # Save and clean up
    canvas.SaveAs(os.path.join(final_output_dir, "plots.png"))
    return cuts_dictionary

# --------------------------------------------------------------------------------------
# 5) main (unchanged)
# --------------------------------------------------------------------------------------
def main():
    # Same as original

    if __name__ == "__main__":
    main()