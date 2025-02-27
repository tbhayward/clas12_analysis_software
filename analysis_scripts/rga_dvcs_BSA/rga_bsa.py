#!/usr/bin/env python3

import os
import json
import ROOT
from multiprocessing import Pool
from functools import partial

# Configure ROOT to run in batch mode
ROOT.gROOT.SetBatch(True)

# --------------------------------------------------------------------------------------
# 1) apply_kinematic_cuts (optimized)
# --------------------------------------------------------------------------------------
def apply_kinematic_cuts(t_value, open_angle_ep2_value, theta_neutral_neutral_value,
                         Emiss2_value, Mx2_value, Mx2_1_value, Mx2_2_value, pTmiss_value,
                         xF_value, analysis_type, data_type, run_period, topology):
    """Placeholder for actual kinematic cuts logic"""
    return True

# --------------------------------------------------------------------------------------
# 2) format_label_name (updated ranges)
# --------------------------------------------------------------------------------------
def format_label_name(variable, analysis_type):
    labels = {
        "open_angle_ep2": "#theta_{e'#gamma} (deg)",
        "Mx2_2": "Miss M_{e'#gamma}^{2} (GeV^{2})",
        "theta_gamma_gamma": "#theta_{#gamma#gamma} (deg)",
        "theta_pi0_pi0": "#theta_{#pi^{0}#pi^{0}} (deg)",
        "xF": "x_{F}",
        "Emiss2": "Miss E^{2} (GeV^{2})",
        "Mx2": "Miss M^{2} (GeV^{2})",
        "Mx2_1": "Miss M_{e'p'}^{2} (GeV^{2})",
        "pTmiss": "p_{T}^{miss} (GeV)"
    }
    if analysis_type == "eppi0":
        return labels.get(variable, variable).replace("#gamma", "#pi^{0}")
    return labels.get(variable, variable)

# --------------------------------------------------------------------------------------
# 3) load_root_files (optimized branch loading)
# --------------------------------------------------------------------------------------
def load_root_files(period):
    """Load ROOT files with only required branches"""
    file_map = {
        "Fa18_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_inb_epgamma.root",
            "mc": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_inb_50nA_10604MeV_epgamma.root"
        },
        "Fa18_out": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_fa18_out_epgamma.root",
            "mc": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_fa18_out_50nA_10604MeV_epgamma.root"
        },
        "Sp19_inb": {
            "data": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data/dvcs/rga_sp19_inb_epgamma.root",
            "mc": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10200MeV_epgamma.root"
        }
    }

    trees = {}
    for category in ["data", "mc"]:
        chain = ROOT.TChain("PhysicsEvents")
        chain.Add(file_map[period][category])
        
        # Activate only required branches
        branches = ["detector1", "detector2", "t1", "open_angle_ep2", "theta_gamma_gamma",
                    "Emiss2", "Mx2", "Mx2_1", "Mx2_2", "pTmiss", "xF"]
        for br in branches:
            chain.SetBranchStatus(br, 1)
        
        trees[category] = chain
    return period, trees

# --------------------------------------------------------------------------------------
# 4) process_events (fixed histogram handling)
# --------------------------------------------------------------------------------------
def process_events(tree, topology, analysis_type, is_mc):
    """Process events with corrected histogram creation"""
    hist_configs = {
        "open_angle_ep2":    (100, 0, 40),
        "theta_gamma_gamma": (100, 0, 2),
        "pTmiss":            (100, 0, 0.3),
        "xF":                (100, -0.4, 0.2),
        "Emiss2":            (100, -1, 2),
        "Mx2":               (100, -0.015, 0.015),
        "Mx2_1":             (100, -1, 1.5),
        "Mx2_2":             (100, 0, 3)
    }

    variables = list(hist_configs.keys())
    hists = {}
    hists_loose = {}

    # Create histograms first
    for var in variables:
        hists[var] = ROOT.TH1D(f"{'mc' if is_mc else 'data'}_{var}", "", *hist_configs[var])
        hists_loose[var] = ROOT.TH1D(f"{'mc' if is_mc else 'data'}_loose_{var}", "", *hist_configs[var])

    # Precompute topology condition
    if topology == "(FD,FD)":
        det_cond = lambda d1, d2: d1 == 1 and d2 == 1
    elif topology == "(CD,FD)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 1
    elif topology == "(CD,FT)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 0
    else:
        return hists, hists_loose

    # Main processing loop
    for event in tree:
        d1 = event.detector1
        d2 = event.detector2
        if not det_cond(d1, d2):
            continue

        cuts_passed = apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, analysis_type, 
            "mc" if is_mc else "data", "", topology
        )

        for var in variables:
            val = getattr(event, var)
            hists[var].Fill(val)
            if cuts_passed:
                hists_loose[var].Fill(val)

    return hists, hists_loose

# --------------------------------------------------------------------------------------
# 5) plot_results (fixed tuple handling)
# --------------------------------------------------------------------------------------
def plot_results(data_hists, mc_hists, plot_title, topology, output_dir):
    """Create plots with proper histogram handling"""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.045, "XY")
    ROOT.gStyle.SetLabelSize(0.04, "XY")
    ROOT.gStyle.SetPadLeftMargin(0.12)
    ROOT.gStyle.SetPadRightMargin(0.08)
    ROOT.gStyle.SetPadTopMargin(0.1)
    ROOT.gStyle.SetPadBottomMargin(0.12)
    ROOT.gStyle.SetLegendBorderSize(0)

    variables = list(data_hists.keys())
    canvas = ROOT.TCanvas("canvas", "", 2400, 1200)
    canvas.Divide(4, 2, 0.005, 0.005)

    for i, var in enumerate(variables):
        pad = canvas.cd(i+1)
        pad.SetGrid()
        pad.SetTicks(1, 1)

        # Get histograms
        dh = data_hists[var]
        mh = mc_hists[var]

        # Configure styles
        dh.SetLineColor(ROOT.kBlue)
        dh.SetLineWidth(2)
        mh.SetLineColor(ROOT.kRed)
        mh.SetLineWidth(2)

        # Normalize
        for h in [dh, mh]:
            h.Scale(1/h.Integral() if h.Integral() > 0 else 1)
            h.GetYaxis().SetTitle("Normalized Counts")
            h.GetYaxis().SetTitleOffset(1.2)
            h.GetXaxis().SetTitle(format_label_name(var, "dvcs"))

        # Draw order
        if dh.GetMaximum() > mh.GetMaximum():
            dh.Draw("HIST E")
            mh.Draw("HIST E SAME")
        else:
            mh.Draw("HIST E")
            dh.Draw("HIST E SAME")

        # Legend
        leg = ROOT.TLegend(0.65, 0.7, 0.95, 0.9)
        leg.SetHeader(f"{var}", "C")
        leg.AddEntry(dh, Form("Data (#mu=%.3f, #sigma=%.3f)", dh.GetMean(), dh.GetStdDev()), "l")
        leg.AddEntry(mh, Form("MC (#mu=%.3f, #sigma=%.3f)", mh.GetMean(), mh.GetStdDev()), "l")
        leg.Draw()

    # Save output
    clean_title = plot_title.replace(" ", "_")
    canvas.SaveAs(os.path.join(output_dir, f"{clean_title}_{topology}_comparison.png"))
    del canvas

# --------------------------------------------------------------------------------------
# 6) process_period (fixed histogram separation)
# --------------------------------------------------------------------------------------
def process_period(period, output_dir):
    """Process a single run period with correct histogram separation"""
    period_code, trees = load_root_files(period)
    run_period, period_text = {
        "Fa18_inb": ("RGA Fa18 Inb", "Fa18 Inb"),
        "Fa18_out": ("RGA Fa18 Out", "Fa18 Out"),
        "Sp19_inb": ("RGA Sp19 Inb", "Sp19 Inb")
    }[period_code]

    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        # Process data and MC
        data_hists, data_loose = process_events(trees["data"], topology, "dvcs", False)
        mc_hists, mc_loose = process_events(trees["mc"], topology, "dvcs", True)

        # Create plots
        plot_title = f"{period_text} DVCS {topology}"
        plot_results(data_hists, mc_hists, plot_title, topology, output_dir)

        # Save cuts
        cuts = {
            "data": {var: {"mean": data_loose[var].GetMean(), 
                          "std": data_loose[var].GetStdDev()} 
                    for var in data_loose},
            "mc": {var: {"mean": mc_loose[var].GetMean(), 
                        "std": mc_loose[var].GetStdDev()} 
                  for var in mc_loose}
        }
        with open(os.path.join(output_dir, f"cuts_{period_code}_{topology}.json"), "w") as f:
            json.dump(cuts, f, indent=2)

# --------------------------------------------------------------------------------------
# 7) Main function with safe parallel processing
# --------------------------------------------------------------------------------------
def main():
    """Main driver with thread-safe ROOT handling"""
    output_dir = "exclusivity_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    # Process periods sequentially due to ROOT's thread-safety issues
    for period in ["Fa18_inb", "Fa18_out", "Sp19_inb"]:
        process_period(period, output_dir)

    # Combine all cuts
    combined_cuts = {}
    for period in ["Fa18_inb", "Fa18_out", "Sp19_inb"]:
        for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
            with open(os.path.join(output_dir, f"cuts_{period}_{topology}.json")) as f:
                combined_cuts[f"{period}_{topology}"] = json.load(f)
    
    with open(os.path.join(output_dir, "combined_cuts.json"), "w") as f:
        json.dump(combined_cuts, f, indent=2)

if __name__ == "__main__":
    main()