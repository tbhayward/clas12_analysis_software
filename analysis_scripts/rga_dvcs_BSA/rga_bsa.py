#!/usr/bin/env python3

import os
import json
import math
import ROOT

# Make sure ROOT doesn't pop up any GUI windows:
ROOT.gROOT.SetBatch(True)

# --------------------------------------------------------------------------------------
# 1) apply_kinematic_cuts
# --------------------------------------------------------------------------------------
def apply_kinematic_cuts(t_value,
                         open_angle_ep2_value,
                         theta_neutral_neutral_value,
                         Emiss2_value,
                         Mx2_value,
                         Mx2_1_value,
                         Mx2_2_value,
                         pTmiss_value,
                         xF_value,
                         analysis_type,
                         data_type,
                         run_period,
                         topology):
    """
    Placeholder function that implements your real logic for 'apply_kinematic_cuts'.
    Currently returns True for all events.
    Replace with your actual cut thresholds/conditions as needed.
    """
    return True
#enddef


# --------------------------------------------------------------------------------------
# 2) format_label_name
# --------------------------------------------------------------------------------------
def format_label_name(variable, analysis_type):
    """
    Mimics your C++ 'formatLabelName' function with LaTeX-style replacements.
    We define a dictionary of common labels for both DVCS and eppi0. If analysisType
    is "eppi0", we do a string replacement of '#gamma' -> '#pi^{0}'.
    """
    common_labels = {
        "open_angle_ep2"    : "#theta_{e'#gamma} (deg)",
        "Mx2_2"             : "Miss M_{e'#gamma}^{2} (GeV^{2})",
        "theta_gamma_gamma" : "#theta_{#gamma(det)#gamma(rec)}",
        "theta_pi0_pi0"     : "#theta_{#pi^{0}(det)#pi^{0}(rec)}",
        "xF"                : "x_{F}",
        "Emiss2"            : "Miss E_{e'p'#gamma} (GeV)",
        "Mx2"               : "Miss M_{e'p'#gamma}^{2} (GeV^{2})",
        "Mx2_1"             : "Miss M_{e'p'}^{2} (GeV^{2})",
        "pTmiss"            : "Miss P_{T(e'p'#gamma)} (GeV)"
    }

    label = common_labels.get(variable, variable)  # fallback if not in dict

    if analysis_type == "eppi0":
        # Replace any occurrences of "#gamma" with "#pi^{0}"
        label = label.replace("#gamma", "#pi^{0}")
    #endif

    return label
#enddef


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
#enddef


# --------------------------------------------------------------------------------------
# 4) determine_exclusivity_cuts
# --------------------------------------------------------------------------------------
def determine_exclusivity_cuts(analysis_type,
                               topology,
                               data_tree,
                               mc_tree,
                               output_dir,
                               run_period,
                               plot_title):
    """
    Creates TCanvas, defines variables, loops over events to fill histograms,
    applies "loose" kinematic cuts, normalizes, draws, and saves histograms.
    Returns a dictionary of mean ± 3σ cuts for each variable in data vs. MC.
    """

    # 1) Global style options
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.05, "XY")
    ROOT.gStyle.SetLabelSize(0.04, "XY")
    ROOT.gStyle.SetLegendTextSize(0.04)

    # 2) Subdirectory under output_dir for dvcs or eppi0
    analysis_dir = "dvcs" if analysis_type == "dvcs" else "eppi0"
    final_output_dir = os.path.join(output_dir, analysis_dir)
    os.makedirs(final_output_dir, exist_ok=True)

    # Create canvases (2x4 = 8 pads)
    canvas = ROOT.TCanvas("exclusivity_canvas", "Exclusivity Plots", 1600, 800)
    canvas.Divide(4, 2)

    canvas_loose_cuts = ROOT.TCanvas("exclusivity_canvas_loose_cuts", 
                                     "Exclusivity Plots with Kinematic Cuts", 
                                     1600, 800)
    canvas_loose_cuts.Divide(4, 2)

    # Variables
    if analysis_type == "dvcs":
        variables = [
            "open_angle_ep2", 
            "theta_gamma_gamma", 
            "pTmiss", 
            "xF", 
            "Emiss2", 
            "Mx2", 
            "Mx2_1", 
            "Mx2_2"
        ]
    else:
        # eppi0 (not in use)
        variables = [
            "open_angle_ep2",
            "theta_pi0_pi0",
            "pTmiss",
            "xF",
            "Emiss2",
            "Mx2",
            "Mx2_1",
            "Mx2_2"
        ]
    #endif

    # Dictionary to store final cuts
    cuts_dictionary = {"data": {}, "mc": {}}

    # Hist binning
    hist_configs = {
        "open_angle_ep2"   : (100, 0.0, 40.0),
        "theta_gamma_gamma": (100, 0.0, 40.0),
        "theta_pi0_pi0"    : (100, 0.0, 40.0),
        "pTmiss"           : (100, 0.0, 1.5),
        "xF"               : (100, -1.0, 1.0),
        "Emiss2"           : (100, -0.5, 1.0),
        "Mx2"              : (100, -0.5, 1.5),
        "Mx2_1"            : (100, -0.5, 1.5),
        "Mx2_2"            : (100, -0.5, 1.5)
    }

    # Loop over variables
    for i, var in enumerate(variables):
        bins, hmin, hmax = hist_configs[var]

        hist_data       = ROOT.TH1D(f"data_{var}",       "", bins, hmin, hmax)
        hist_mc         = ROOT.TH1D(f"mc_{var}",         "", bins, hmin, hmax)
        hist_data_loose = ROOT.TH1D(f"data_loose_{var}", "", bins, hmin, hmax)
        hist_mc_loose   = ROOT.TH1D(f"mc_loose_{var}",   "", bins, hmin, hmax)

        for hist_ in [hist_data, hist_mc, hist_data_loose, hist_mc_loose]:
            hist_.SetDirectory(0)
        #endfor

        # Style
        hist_data.SetMarkerStyle(20)
        hist_data.SetMarkerColor(ROOT.kBlue)
        hist_data.SetLineColor(ROOT.kBlue)
        hist_data.SetMarkerSize(0.5)

        hist_mc.SetMarkerStyle(21)
        hist_mc.SetMarkerColor(ROOT.kRed)
        hist_mc.SetLineColor(ROOT.kRed)
        hist_mc.SetMarkerSize(0.5)

        hist_data_loose.SetMarkerStyle(20)
        hist_data_loose.SetMarkerColor(ROOT.kBlue)
        hist_data_loose.SetLineColor(ROOT.kBlue)
        hist_data_loose.SetMarkerSize(0.5)

        hist_mc_loose.SetMarkerStyle(21)
        hist_mc_loose.SetMarkerColor(ROOT.kRed)
        hist_mc_loose.SetLineColor(ROOT.kRed)
        hist_mc_loose.SetMarkerSize(0.5)

        # ---------------------
        # Fill data histograms
        # ---------------------
        for event in data_tree:
            det1 = getattr(event, "detector1", -1)
            det2 = getattr(event, "detector2", -1)

            # Match chosen topology
            if (topology == "(FD,FD)" and det1 == 1 and det2 == 1) or \
               (topology == "(CD,FD)" and det1 == 2 and det2 == 1) or \
               (topology == "(CD,FT)" and det1 == 2 and det2 == 0):

                val_data = getattr(event, var, None)
                if val_data is None:
                    continue
                hist_data.Fill(val_data)

                # "loose" cuts
                t_val_data    = getattr(event, "t1", 0)
                open_ang_data = getattr(event, "open_angle_ep2", 0)
                if analysis_type == "dvcs":
                    theta_nn_data = getattr(event, "theta_gamma_gamma", 0)
                else:
                    theta_nn_data = getattr(event, "theta_pi0_pi0", 0)

                Emiss2_d  = getattr(event, "Emiss2", 0)
                Mx2_d     = getattr(event, "Mx2", 0)
                Mx2_1_d   = getattr(event, "Mx2_1", 0)
                Mx2_2_d   = getattr(event, "Mx2_2", 0)
                pTmiss_d  = getattr(event, "pTmiss", 0)
                xF_d      = getattr(event, "xF", 0)

                # If your logic says the event passes cuts:
                if apply_kinematic_cuts(t_val_data,
                                        open_ang_data,
                                        theta_nn_data,
                                        Emiss2_d,
                                        Mx2_d,
                                        Mx2_1_d,
                                        Mx2_2_d,
                                        pTmiss_d,
                                        xF_d,
                                        analysis_type,
                                        "data",
                                        run_period,
                                        topology):
                    hist_data_loose.Fill(val_data)
                #endif
            #endif
        #endfor

        # ---------------
        # Fill MC histos
        # ---------------
        for event in mc_tree:
            det1_mc = getattr(event, "detector1", -1)
            det2_mc = getattr(event, "detector2", -1)

            if (topology == "(FD,FD)" and det1_mc == 1 and det2_mc == 1) or \
               (topology == "(CD,FD)" and det1_mc == 2 and det2_mc == 1) or \
               (topology == "(CD,FT)" and det1_mc == 2 and det2_mc == 0):

                val_mc = getattr(event, var, None)
                if val_mc is None:
                    continue
                hist_mc.Fill(val_mc)

                t_val_mc    = getattr(event, "t1", 0)
                open_ang_mc = getattr(event, "open_angle_ep2", 0)
                if analysis_type == "dvcs":
                    theta_nn_mc = getattr(event, "theta_gamma_gamma", 0)
                else:
                    theta_nn_mc = getattr(event, "theta_pi0_pi0", 0)

                Emiss2_m = getattr(event, "Emiss2", 0)
                Mx2_m    = getattr(event, "Mx2", 0)
                Mx2_1_m  = getattr(event, "Mx2_1", 0)
                Mx2_2_m  = getattr(event, "Mx2_2", 0)
                pTmiss_m = getattr(event, "pTmiss", 0)
                xF_m     = getattr(event, "xF", 0)

                if apply_kinematic_cuts(t_val_mc,
                                        open_ang_mc,
                                        theta_nn_mc,
                                        Emiss2_m,
                                        Mx2_m,
                                        Mx2_1_m,
                                        Mx2_2_m,
                                        pTmiss_m,
                                        xF_m,
                                        analysis_type,
                                        "mc",
                                        run_period,
                                        topology):
                    hist_mc_loose.Fill(val_mc)
                #endif
            #endif
        #endfor

        # ----------------
        # Normalize hists
        # ----------------
        if hist_data.Integral() > 0:
            hist_data.Scale(1.0 / hist_data.Integral())
        if hist_mc.Integral() > 0:
            hist_mc.Scale(1.0 / hist_mc.Integral())
        if hist_data_loose.Integral() > 0:
            hist_data_loose.Scale(1.0 / hist_data_loose.Integral())
        if hist_mc_loose.Integral() > 0:
            hist_mc_loose.Scale(1.0 / hist_mc_loose.Integral())

        # y-range
        y_max_no_cuts = 1.75 * max(hist_data.GetMaximum(), hist_mc.GetMaximum())
        y_max_loose   = 1.75 * max(hist_data_loose.GetMaximum(), hist_mc_loose.GetMaximum())

        # -------------------------------------------------
        # Draw "no cuts" histograms in pad i+1 of first canvas
        # -------------------------------------------------
        canvas.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)

        hist_data.SetTitle(plot_title)
        hist_data.SetXTitle(format_label_name(var, analysis_type))
        hist_data.SetYTitle("Normalized counts")
        hist_data.GetYaxis().SetRangeUser(0, y_max_no_cuts)
        hist_data.Draw("E1")
        hist_mc.Draw("E1 SAME")

        legend_nocuts = ROOT.TLegend(0.45, 0.70, 0.85, 0.90)
        legend_nocuts.AddEntry(hist_data, f"Data ({int(hist_data.GetEntries())} evts)", "lep")
        legend_nocuts.AddEntry(hist_mc,   f"MC ({int(hist_mc.GetEntries())} evts)",   "lep")
        legend_nocuts.SetTextSize(0.03)
        legend_nocuts.Draw()

        # -----------------------------------------------------
        # Draw "loose cuts" histograms in pad i+1 of 2nd canvas
        # -----------------------------------------------------
        canvas_loose_cuts.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)

        hist_data_loose.SetTitle(plot_title)
        hist_data_loose.SetXTitle(format_label_name(var, analysis_type))
        hist_data_loose.SetYTitle("Normalized counts")
        hist_data_loose.GetYaxis().SetRangeUser(0, y_max_loose)
        hist_data_loose.Draw("E1")
        hist_mc_loose.Draw("E1 SAME")

        mean_data = hist_data_loose.GetMean()
        std_data  = hist_data_loose.GetStdDev()
        mean_mc   = hist_mc_loose.GetMean()
        std_mc    = hist_mc_loose.GetStdDev()

        data_min_cut = mean_data - 3.0 * std_data
        data_max_cut = mean_data + 3.0 * std_data
        mc_min_cut   = mean_mc   - 3.0 * std_mc
        mc_max_cut   = mean_mc   + 3.0 * std_mc

        cuts_dictionary["data"][var] = {
            "mean":    mean_data,
            "std":     std_data,
            "min_cut": data_min_cut,
            "max_cut": data_max_cut
        }
        cuts_dictionary["mc"][var] = {
            "mean":    mean_mc,
            "std":     std_mc,
            "min_cut": mc_min_cut,
            "max_cut": mc_max_cut
        }

        legend_loose = ROOT.TLegend(0.35, 0.65, 0.85, 0.90)
        # Note: Use #mu and #sigma for Greek letters in ROOT
        legend_loose.AddEntry(hist_data_loose, 
                              f"Data (#mu={mean_data:.3f}, #sigma={std_data:.3f})", "lep")
        legend_loose.AddEntry(hist_mc_loose, 
                              f"MC   (#mu={mean_mc:.3f}, #sigma={std_mc:.3f})", "lep")
        legend_loose.SetTextSize(0.03)
        legend_loose.Draw()
    #endfor

    # Save canvases
    clean_title = plot_title.replace(" ", "_")
    outname_nocuts = os.path.join(
        final_output_dir, 
        f"exclusivity_plots_{clean_title}_{topology}_no_cuts.png"
    )
    outname_cuts   = os.path.join(
        final_output_dir, 
        f"exclusivity_plots_{clean_title}_{topology}_kinematic_cuts.png"
    )

    canvas.SaveAs(outname_nocuts)
    canvas_loose_cuts.SaveAs(outname_cuts)

    # Cleanup
    del canvas
    del canvas_loose_cuts

    return cuts_dictionary
#enddef


# --------------------------------------------------------------------------------------
# 5) main
# --------------------------------------------------------------------------------------
def main():
    """
    Main driver: no parallel processing, just a loop over run periods and topologies.
    """

    # Use a single top-level output directory
    output_dir = "exclusivity_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Load the ROOT files (DVCS data + MC)
    tree_dict = load_root_files()

    # We'll just do DVCS:
    analysis_type = "dvcs"
    topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    dvcs_period_keys = ["Fa18_inb", "Fa18_out", "Sp19_inb"]

    combined_cuts = {}

    # Loop over run periods
    for rp in dvcs_period_keys:
        # Convert rp => user-friendly run_period and plot title snippet
        if rp == "Fa18_inb":
            run_period  = "RGA Fa18 Inb"
            period_text = "Fa18 Inb"
        elif rp == "Fa18_out":
            run_period  = "RGA Fa18 Out"
            period_text = "Fa18 Out"
        elif rp == "Sp19_inb":
            run_period  = "RGA Sp19 Inb"
            period_text = "Sp19 Inb"
        else:
            run_period  = "Unknown"
            period_text = rp
        #endif

        # Retrieve the DVCS data + MC TChains for this run period
        data_tree = tree_dict["dvcs_data"][rp]
        mc_tree   = tree_dict["dvcs_mc"][rp]

        for topo in topologies:
            # Build a descriptive plot title
            plot_title = f"{period_text} DVCS"

            # Perform the exclusivity analysis for (analysis_type, topo, data_tree, mc_tree)
            cuts_dict = determine_exclusivity_cuts(
                analysis_type=analysis_type,
                topology=topo,
                data_tree=data_tree,
                mc_tree=mc_tree,
                output_dir=output_dir,
                run_period=run_period,
                plot_title=plot_title
            )

            # Store the returned dictionary of cuts
            combined_cuts[(analysis_type, rp, topo)] = cuts_dict
        #endfor
    #endfor

    # Write out combined cuts to JSON
    out_json = os.path.join(output_dir, "combined_exclusivity_cuts.json")
    jdict = {}
    for (analysis, rp, topo), val in combined_cuts.items():
        key_str = f"{analysis}__{rp}__{topo}"
        jdict[key_str] = val
    #endfor

    with open(out_json, "w") as f:
        json.dump(jdict, f, indent=4)
    #endwith

    print(f"All DVCS cuts saved to {out_json}")
#enddef


if __name__ == "__main__":
    main()
#endif