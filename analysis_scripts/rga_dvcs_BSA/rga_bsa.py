#!/usr/bin/env python3

import os
import json
import ROOT

# Configure ROOT to run in batch mode
ROOT.gROOT.SetBatch(True)

# ------------------------------ Global Style Settings ------------------------------
def configure_global_style():
    style = ROOT.gStyle
    style.SetOptStat(0)
    style.SetTitleSize(0.045, "XY")
    style.SetLabelSize(0.04, "XY")
    style.SetPadLeftMargin(0.12)
    style.SetPadRightMargin(0.08)
    style.SetPadTopMargin(0.08)
    style.SetPadBottomMargin(0.12)
    style.SetLegendBorderSize(1)
    style.SetLegendTextSize(0.035)
    style.SetLegendFillColor(0)
    style.SetGridColor(ROOT.kGray+1)
    style.SetGridStyle(2)
    style.SetGridWidth(1)
#enddef

configure_global_style()

# --------------------------------------------------------------------------------------
# 1) apply_kinematic_cuts (placeholder + any additional logic you desire)
# --------------------------------------------------------------------------------------
def apply_kinematic_cuts(t_value, open_angle_ep2_value, theta_neutral_neutral_value,
                         Emiss2_value, Mx2_value, Mx2_1_value, Mx2_2_value, pTmiss_value,
                         xF_value, analysis_type, data_type, run_period, topology):
    """
    Placeholder for your actual physics kinematic cuts logic.
    Currently always returns True.
    """
    return True
#enddef

# --------------------------------------------------------------------------------------
# 1b) passes_3sigma_cuts (new function to enforce mu +/- 3*sigma from JSON)
# --------------------------------------------------------------------------------------
def passes_3sigma_cuts(event, is_mc, cuts_dict):
    """
    Enforce a mu +/- 3*sigma cut for *each* variable in the JSON dictionary.
    We must pick the 'data' or 'mc' sub-dict depending on is_mc.
    If event.var is outside mu +/- 3*sigma for any variable, return False.
    """
    category = "mc" if is_mc else "data"
    for var, stats in cuts_dict[category].items():
        mean = stats["mean"]
        stdv = stats["std"]
        val = getattr(event, var, None)
        # If this variable isn't in the event tree for some reason, skip safely.
        if val is None:
            continue
        #endfor
        if val < (mean - 3.0*stdv) or val > (mean + 3.0*stdv):
            return False
        #endif
    #endfor
    return True
#enddef

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
#enddef

# --------------------------------------------------------------------------------------
# 3) load_root_files (optimized branch loading)
# --------------------------------------------------------------------------------------
def load_root_files(period):
    """Load ROOT files with only required branches."""
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

# --------------------------------------------------------------------------------------
# 4) process_events (first pass, no JSON cuts)
# --------------------------------------------------------------------------------------
def process_events(tree, topology, analysis_type, is_mc):
    """Process events (first pass) with only the original kinematic_cuts."""
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

    # Create histograms
    for var in variables:
        hname = f"{'mc' if is_mc else 'data'}_{var}"
        hname_loose = f"{'mc' if is_mc else 'data'}_loose_{var}"
        hists[var] = ROOT.TH1D(hname, "", *hist_configs[var])
        hists_loose[var] = ROOT.TH1D(hname_loose, "", *hist_configs[var])
    #endfor

    # Precompute topology condition
    if topology == "(FD,FD)":
        det_cond = lambda d1, d2: d1 == 1 and d2 == 1
    elif topology == "(CD,FD)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 1
    elif topology == "(CD,FT)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 0
    else:
        return hists, hists_loose
    #endif

    # Main loop over events
    for event in tree:
        d1 = event.detector1
        d2 = event.detector2
        if not det_cond(d1, d2):
            continue
        #endif

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
            #endif
        #endfor
    #endfor

    return hists, hists_loose
#enddef

# --------------------------------------------------------------------------------------
# 4b) process_events_with_json_cuts (second pass, uses mu +/- 3 sigma)
# --------------------------------------------------------------------------------------
def process_events_with_json_cuts(tree, topology, analysis_type, is_mc, cuts_dict):
    """
    Process events (second pass) using:
      1) The original kinematic_cuts (placeholder).
      2) Additional mu +/- 3*sigma from the loaded JSON dictionary.
    """
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

    # Create histograms
    for var in variables:
        hname = f"{'mc' if is_mc else 'data'}_3sigma_{var}"
        hists[var] = ROOT.TH1D(hname, "", *hist_configs[var])
    #endfor

    # Precompute topology condition
    if topology == "(FD,FD)":
        det_cond = lambda d1, d2: d1 == 1 and d2 == 1
    elif topology == "(CD,FD)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 1
    elif topology == "(CD,FT)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 0
    else:
        return hists
    #endif

    # Main loop over events
    for event in tree:
        d1 = event.detector1
        d2 = event.detector2
        if not det_cond(d1, d2):
            continue
        #endif

        # First apply the original placeholder kinematic cuts
        basic_cuts_passed = apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, analysis_type,
            "mc" if is_mc else "data", "", topology
        )
        if not basic_cuts_passed:
            continue
        #endif

        # Now enforce the mu +/- 3*sigma cuts from JSON
        if not passes_3sigma_cuts(event, is_mc, cuts_dict):
            continue
        #endif

        # If we get here, the event passes all criteria
        for var in variables:
            val = getattr(event, var)
            hists[var].Fill(val)
        #endfor
    #endfor

    return hists
#enddef

# --------------------------------------------------------------------------------------
# 5) plot_results (shared between first pass and second pass)
# --------------------------------------------------------------------------------------
def plot_results(data_hists, mc_hists, plot_title, topology, output_dir):
    """
    Create publication-quality plots with proper styling.
      - Normalizes each histogram before computing the maximum.
      - Sets the y-axis range to 1.2x the largest normalized bin.
      - Plots data in blue, MC in red.
      - Appends each histogram’s (mu, sigma) to the legend entry.
    """
    variables = list(data_hists.keys())
    canvas = ROOT.TCanvas("canvas", "", 2400, 1200)
    canvas.Divide(4, 2, 0.002, 0.002)
    
    # Shifted the legend to the left to avoid clipping mu/sigma text
    base_legend = ROOT.TLegend(0.45, 0.75, 0.92, 0.89)
    base_legend.SetFillStyle(0)
    base_legend.SetTextFont(42)
    base_legend.SetBorderSize(1)
    base_legend.SetMargin(0.12)

    for i, var in enumerate(variables):
        pad = canvas.cd(i + 1)
        pad.SetTicks(1, 1)  # Ticks only (no grid).

        # Retrieve histograms
        dh = data_hists[var]
        mh = mc_hists[var]

        # Set colors/styles: Data in blue, MC in red
        dh.SetLineColor(ROOT.kBlue)
        dh.SetMarkerColor(ROOT.kBlue)
        dh.SetLineWidth(3)
        dh.SetMarkerStyle(20)
        dh.SetMarkerSize(1.2)

        mh.SetLineColor(ROOT.kRed)
        mh.SetMarkerColor(ROOT.kRed)
        mh.SetLineWidth(3)
        mh.SetMarkerStyle(21)
        mh.SetMarkerSize(1.2)

        # Scale (normalize) the histograms before setting maximum
        for h in [dh, mh]:
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())
            #endif
            h.GetYaxis().SetTitle("Normalized Counts")
            h.GetYaxis().SetTitleOffset(1.4)
            h.GetXaxis().SetTitle(format_label_name(var, "dvcs"))
        #endfor

        # Now find the largest bin content after scaling
        max_val = max(dh.GetMaximum(), mh.GetMaximum()) * 1.2
        dh.SetMaximum(max_val)
        mh.SetMaximum(max_val)

        # Compute mu and sigma for each
        mu_data = dh.GetMean()
        sigma_data = dh.GetStdDev()
        mu_mc   = mh.GetMean()
        sigma_mc = mh.GetStdDev()

        # Draw the histograms with error bars
        dh.Draw("E1")
        mh.Draw("E1 SAME")

        # Clone and position a legend for each pad; add mu and sigma text
        pad_leg = base_legend.Clone()
        pad_leg.AddEntry(dh, f"Data (#mu={mu_data:.3f}, #sigma={sigma_data:.3f})", "lep")
        pad_leg.AddEntry(mh, f"MC (#mu={mu_mc:.3f}, #sigma={sigma_mc:.3f})", "lep")
        pad_leg.Draw()
    #endfor

    # Save output
    clean_title = plot_title.replace(" ", "_").replace("(", "").replace(")", "")
    canvas.SaveAs(os.path.join(output_dir, f"{clean_title}_{topology}_comparison.png"))
    del canvas
#enddef

# --------------------------------------------------------------------------------------
# 6) process_period (with second pass after saving JSON)
# --------------------------------------------------------------------------------------
def process_period(period, output_dir):
    """Process a single run period with resource awareness."""
    print(f"⏳ Starting processing for {period}...")
    
    period_code, trees = load_root_files(period)
    run_info = {
        "Fa18_inb": ("RGA Fa18 Inb", "Fa18 Inb"),
        "Fa18_out": ("RGA Fa18 Out", "Fa18 Out"),
        "Sp19_inb": ("RGA Sp19 Inb", "Sp19 Inb")
    }[period_code]

    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        print(f"  🔄 Processing topology {topology}")
        
        # 1) First pass: fill histograms, produce comparison plots, save JSON
        data_hists, data_loose = process_events(trees["data"], topology, "dvcs", False)
        mc_hists, mc_loose = process_events(trees["mc"], topology, "dvcs", True)

        # Create the first comparison plot
        plot_title = f"{run_info[1]} DVCS"
        plot_results(data_hists, mc_hists, plot_title, topology, output_dir)

        # Save the cuts (mu,std) from the "loose" histograms
        save_cuts(period_code, topology, output_dir, data_loose, mc_loose)

        # 2) Second pass: read the JSON, enforce mu +/- 3*sigma, produce new plots
        json_path = os.path.join(output_dir, f"cuts_{period_code}_{topology}.json")
        if not os.path.exists(json_path):
            print(f"   ⚠️ JSON cuts file not found, skipping second pass: {json_path}")
            continue
        #endif

        with open(json_path, "r") as infile:
            cuts_dict = json.load(infile)
        #endwith

        # Apply 3σ cuts
        data_hists_3s = process_events_with_json_cuts(trees["data"], topology, "dvcs", False, cuts_dict)
        mc_hists_3s   = process_events_with_json_cuts(trees["mc"],   topology, "dvcs", True,  cuts_dict)

        # Plot results with 3σ cuts
        plot_title_3s = f"{run_info[1]} DVCS (3sigma Cuts)"
        plot_results(data_hists_3s, mc_hists_3s, plot_title_3s, topology, output_dir)
    #endfor

    print(f"✅ Completed {period}\n")
#enddef

# --------------------------------------------------------------------------------------
# 6b) save_cuts: writes out means and stddev for each variable
# --------------------------------------------------------------------------------------
def save_cuts(period, topology, output_dir, data, mc):
    """Save cuts (mu,sigma) to JSON with error handling."""
    cuts = {
        "data": {var: {"mean": data[var].GetMean(), "std": data[var].GetStdDev()}
                 for var in data},
        "mc":   {var: {"mean": mc[var].GetMean(),   "std": mc[var].GetStdDev()}
                 for var in mc}
    }
    try:
        out_path = os.path.join(output_dir, f"cuts_{period}_{topology}.json")
        with open(out_path, "w") as f:
            json.dump(cuts, f, indent=2)
        #endwith
    except Exception as e:
        print(f"⚠️ Error saving cuts for {period}_{topology}: {str(e)}")
    #endtry
#enddef

# --------------------------------------------------------------------------------------
# 7) Main function with improved progress tracking
# --------------------------------------------------------------------------------------
def main():
    """Main driver with enhanced resource management."""
    output_dir = "exclusivity_plots"
    os.makedirs(output_dir, exist_ok=True)
    
    print("🚀 Starting analysis with sequential processing")
    print("❗ Note: ROOT's thread-safety limitations require sequential processing")
    print("      We achieve speed through optimized single-pass event handling\n")
    
    for period in ["Fa18_inb", "Fa18_out", "Sp19_inb"]:
        process_period(period, output_dir)
    #endfor

    print("🧩 Combining results...")
    combine_results(output_dir)
    print("\n🎉 Analysis complete!")
#enddef

def combine_results(output_dir):
    """Combine all JSON cuts into a single file."""
    combined = {}
    for period in ["Fa18_inb", "Fa18_out", "Sp19_inb"]:
        for topology in ["FD_FD", "CD_FD", "CD_FT"]:  # Adjusted for filename consistency
            try:
                fpath = os.path.join(output_dir, f"cuts_{period}_{topology}.json")
                with open(fpath, "r") as f:
                    combined[f"{period}_{topology}"] = json.load(f)
                #endwith
            except FileNotFoundError:
                print(f"⚠️ Missing cuts file for {period}_{topology}")
            #endtry
        #endfor
    #endfor

    combined_path = os.path.join(output_dir, "combined_cuts.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    #endwith
#enddef

if __name__ == "__main__":
    main()
#endif