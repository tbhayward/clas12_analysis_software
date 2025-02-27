#!/usr/bin/env python3

import os
import json
import ROOT

# Configure ROOT to run in batch mode
ROOT.gROOT.SetBatch(True)

# --------------------------------------------------------------------------------------
# 1) apply_kinematic_cuts
# --------------------------------------------------------------------------------------
def apply_kinematic_cuts(t_value, open_angle_ep2_value, theta_neutral_neutral_value,
                         Emiss2_value, Mx2_value, Mx2_1_value, Mx2_2_value, pTmiss_value,
                         xF_value, analysis_type, data_type, run_period, topology):
    """Placeholder for actual kinematic cuts logic"""
    return True

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
    """Load ROOT files into TChains"""
    file_dict = {
        "dvcs_data": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root"
        },
        "dvcs_mc": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"
        }
    }

    tree_dict = {}
    for category, subdict in file_dict.items():
        tree_dict[category] = {}
        for period, filename in subdict.items():
            chain = ROOT.TChain("PhysicsEvents")
            chain.Add(filename)
            tree_dict[category][period] = chain
    return tree_dict

# --------------------------------------------------------------------------------------
# 4) determine_exclusivity_cuts (UPDATED)
# --------------------------------------------------------------------------------------
def determine_exclusivity_cuts(analysis_type, topology, data_tree, mc_tree,
                               output_dir, run_period, plot_title):
    """Main analysis function with corrected event processing"""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.05, "XY")
    ROOT.gStyle.SetLabelSize(0.04, "XY")

    # Create output directory
    analysis_dir = "dvcs" if analysis_type == "dvcs" else "eppi0"
    final_output_dir = os.path.join(output_dir, analysis_dir)
    os.makedirs(final_output_dir, exist_ok=True)

    # Variables and histogram configurations
    variables = ["open_angle_ep2", "theta_gamma_gamma", "pTmiss", "xF",
                 "Emiss2", "Mx2", "Mx2_1", "Mx2_2"]
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

    # Initialize histograms
    hists = {
        'data': {var: ROOT.TH1D(f"data_{var}", "", *hist_configs[var]) for var in variables},
        'mc': {var: ROOT.TH1D(f"mc_{var}", "", *hist_configs[var]) for var in variables},
        'data_loose': {var: ROOT.TH1D(f"data_loose_{var}", "", *hist_configs[var]) for var in variables},
        'mc_loose': {var: ROOT.TH1D(f"mc_loose_{var}", "", *hist_configs[var]) for var in variables}
    }

    # Style configuration
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
        cuts_passed = apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, analysis_type, "data", run_period, topology
        )

        # Fill histograms
        for var in variables:
            val = getattr(event, var)
            hists['data'][var].Fill(val)
            if cuts_passed:
                hists['data_loose'][var].Fill(val)

    # Process MC tree
    for event in mc_tree:
        det1 = event.detector1
        det2 = event.detector2
        if not ((topology == "(FD,FD)" and det1 == 1 and det2 == 1) or
                (topology == "(CD,FD)" and det1 == 2 and det2 == 1) or
                (topology == "(CD,FT)" and det1 == 2 and det2 == 0)):
            continue

        # Get cut variables
        cuts_passed = apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, analysis_type, "mc", run_period, topology
        )

        # Fill histograms
        for var in variables:
            val = getattr(event, var)
            hists['mc'][var].Fill(val)
            if cuts_passed:
                hists['mc_loose'][var].Fill(val)

    # Create canvases
    canvas = ROOT.TCanvas("canvas", "Exclusivity Plots", 1600, 800)
    canvas.Divide(4, 2)
    canvas_loose = ROOT.TCanvas("canvas_loose", "With Kinematic Cuts", 1600, 800)
    canvas_loose.Divide(4, 2)

    cuts_dictionary = {"data": {}, "mc": {}}

    # Draw and process results
    for i, var in enumerate(variables):
        # Normalize
        for key in ['data', 'mc', 'data_loose', 'mc_loose']:
            if hists[key][var].Integral() > 0:
                hists[key][var].Scale(1.0 / hists[key][var].Integral())

        # Original plots
        canvas.cd(i+1)
        hists['data'][var].SetTitle(f"{plot_title} - {var}")
        hists['data'][var].SetXTitle(format_label_name(var, analysis_type))
        hists['data'][var].SetYTitle("Normalized Counts")
        hists['data'][var].Draw("E1")
        hists['mc'][var].Draw("E1 SAME")

        # Loose cuts plots
        canvas_loose.cd(i+1)
        hists['data_loose'][var].SetTitle(f"{plot_title} - {var} (Cuts)")
        hists['data_loose'][var].SetXTitle(format_label_name(var, analysis_type))
        hists['data_loose'][var].SetYTitle("Normalized Counts")
        hists['data_loose'][var].Draw("E1")
        hists['mc_loose'][var].Draw("E1 SAME")

        # Store cuts
        cuts_dictionary["data"][var] = {
            "mean": hists['data_loose'][var].GetMean(),
            "std": hists['data_loose'][var].GetStdDev()
        }
        cuts_dictionary["mc"][var] = {
            "mean": hists['mc_loose'][var].GetMean(),
            "std": hists['mc_loose'][var].GetStdDev()
        }

    # Save outputs
    clean_title = plot_title.replace(" ", "_")
    canvas.SaveAs(os.path.join(final_output_dir, f"{clean_title}_{topology}_nocuts.png"))
    canvas_loose.SaveAs(os.path.join(final_output_dir, f"{clean_title}_{topology}_cuts.png"))

    return cuts_dictionary

# --------------------------------------------------------------------------------------
# 5) Main function
# --------------------------------------------------------------------------------------
def main():
    """Main driver function"""
    output_dir = "exclusivity_plots"
    os.makedirs(output_dir, exist_ok=True)
    tree_dict = load_root_files()

    analysis_type = "dvcs"
    topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    periods = {
        "Fa18_inb": ("RGA Fa18 Inb", "Fa18 Inb"),
        "Fa18_out": ("RGA Fa18 Out", "Fa18 Out"),
        "Sp19_inb": ("RGA Sp19 Inb", "Sp19 Inb")
    }

    combined_cuts = {}

    for period_code, (run_period, period_text) in periods.items():
        data_tree = tree_dict["dvcs_data"][period_code]
        mc_tree = tree_dict["dvcs_mc"][period_code]

        for topology in topologies:
            plot_title = f"{period_text} DVCS {topology}"
            cuts = determine_exclusivity_cuts(
                analysis_type, topology, data_tree, mc_tree,
                output_dir, run_period, plot_title
            )
            combined_cuts[f"{period_code}_{topology}"] = cuts

    # Save combined cuts
    with open(os.path.join(output_dir, "exclusivity_cuts.json"), "w") as f:
        json.dump(combined_cuts, f, indent=2)

if __name__ == "__main__":
    main()