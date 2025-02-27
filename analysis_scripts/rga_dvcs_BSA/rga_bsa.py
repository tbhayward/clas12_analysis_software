#!/usr/bin/env python3

import os
import json
import math
import ROOT
import concurrent.futures  # for parallel processing with futures

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
    This function implements your real logic for 'apply_kinematic_cuts'.
    Right now it's a placeholder that always returns True. 
    Replace with your actual cut thresholds/conditions.
    """
    return True
#enddef


# --------------------------------------------------------------------------------------
# 2) format_label_name
# --------------------------------------------------------------------------------------
def format_label_name(variable, analysis_type):
    """
    Mirrors 'formatLabelName' from your C++ code. 
    Customize if certain variable names should get special LaTeX formatting, etc.
    """
    if variable == "theta_gamma_gamma" and analysis_type == "dvcs":
        return "#theta_{#gamma#gamma} (deg)"
    elif variable == "theta_pi0_pi0" and analysis_type == "eppi0":
        return "#theta_{#pi^{0}#pi^{0}} (deg)"
    elif variable == "pTmiss":
        return "p_{T}^{miss} (GeV/c)"
    else:
        return variable
#enddef


# --------------------------------------------------------------------------------------
# 3) load_root_files
# --------------------------------------------------------------------------------------
def load_root_files():
    """
    Loads all your ROOT files into TChains keyed by (category, run period).
    No lines omitted -- every path you gave is included explicitly.
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
        }, 
        # eppi0 reconstructed mc (for background subtraction):
        "eppi0_mc_recon": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_rga_sp19_inb_50nA_10200MeV_eppi0.root"
        },
        # eppi0 generated mc (for background subtraction):
        "eppi0_mc_gen": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_rga_fa18_inb_50nA_10604MeV_eppi0.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_rga_fa18_out_50nA_10604MeV_eppi0.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/gen_aaogen_norad_rga_sp19_inb_50nA_10200MeV_eppi0.root"
        },
        # eppi0 background files:
        "eppi0_bkg": {
            "Fa18_inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_inb_epgamma.root",
            "Fa18_out": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_fa18_out_epgamma.root",
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/eppi0/data/pass2/mc/hipo_files/eppi0_bkg_aaogen_norad_rga_sp19_inb_epgamma.root"
        }
    }
    #enddict

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
    Replicates your C++ determine_exclusivity(...) logic:
     - Creates canvases
     - Defines variables (dvcs/eppi0)
     - Loops over events to fill histograms
     - Applies "loose" cuts (apply_kinematic_cuts)
     - Normalizes, draws, saves histograms
     - Returns mean±3σ cuts for each variable in data vs. MC
    """

    # 1) Global style options
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.05, "XY")
    ROOT.gStyle.SetLabelSize(0.04, "XY")
    ROOT.gStyle.SetLegendTextSize(0.04)

    # -------------------------------------------------------
    #  FIX #1: Remove the extra "exclusivity_plots" component
    # -------------------------------------------------------
    # Instead of os.path.join(output_dir, "exclusivity_plots", analysis_dir),
    # we just do: final_output_dir = os.path.join(output_dir, analysis_dir)
    # Because output_dir is already "exclusivity_plots" in your main() function
    analysis_dir = "dvcs" if analysis_type == "dvcs" else "eppi0"
    final_output_dir = os.path.join(output_dir, analysis_dir)

    # If multiple processes try to create the same directory, we can use exist_ok=True
    os.makedirs(final_output_dir, exist_ok=True)

    # 3) Create canvases
    canvas = ROOT.TCanvas("exclusivity_canvas", "Exclusivity Plots", 1600, 800)
    canvas.Divide(4, 2)  # 2 rows × 4 columns

    canvas_loose_cuts = ROOT.TCanvas("exclusivity_canvas_loose_cuts", 
                                     "Exclusivity Plots with Kinematic Cuts", 
                                     1600, 800)
    canvas_loose_cuts.Divide(4, 2)

    # 4) Define variables
    if analysis_type == "dvcs":
        variables = ["open_angle_ep2", 
                     "theta_gamma_gamma", 
                     "pTmiss", 
                     "xF", 
                     "Emiss2", 
                     "Mx2", 
                     "Mx2_1", 
                     "Mx2_2"]
    else:
        variables = ["open_angle_ep2",
                     "theta_pi0_pi0",
                     "pTmiss",
                     "xF",
                     "Emiss2",
                     "Mx2",
                     "Mx2_1",
                     "Mx2_2"]
    #endif

    # 6) Prepare dictionary to store final cuts
    cuts_dictionary = {
        "data": {},
        "mc":   {}
    }

    # 7) Hist binning
    hist_configs = {
        "open_angle_ep2":   (100,  0.0,  40.0),
        "theta_gamma_gamma":(100,  0.0,  40.0),
        "theta_pi0_pi0":    (100,  0.0,  40.0),
        "pTmiss":           (100,  0.0,  1.5),
        "xF":               (100, -1.0,  1.0),
        "Emiss2":           (100, -0.5,  1.0),
        "Mx2":              (100, -0.5,  1.5),
        "Mx2_1":            (100, -0.5,  1.5),
        "Mx2_2":            (100, -0.5,  1.5),
    }

    # 8) Loop over variables
    for i, var in enumerate(variables):
        bins, hmin, hmax = hist_configs[var]

        hist_data       = ROOT.TH1D(f"data_{var}",       "", bins, hmin, hmax)
        hist_mc         = ROOT.TH1D(f"mc_{var}",         "", bins, hmin, hmax)
        hist_data_loose = ROOT.TH1D(f"data_loose_{var}", "", bins, hmin, hmax)
        hist_mc_loose   = ROOT.TH1D(f"mc_loose_{var}",   "", bins, hmin, hmax)

        hist_data.SetDirectory(0)
        hist_mc.SetDirectory(0)
        hist_data_loose.SetDirectory(0)
        hist_mc_loose.SetDirectory(0)

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

        # Fill data histograms
        for event in data_tree:
            det1 = getattr(event, "detector1")
            det2 = getattr(event, "detector2")

            if (topology == "(FD,FD)" and det1 == 1 and det2 == 1) or \
               (topology == "(CD,FD)" and det1 == 2 and det2 == 1) or \
               (topology == "(CD,FT)" and det1 == 2 and det2 == 0):
                val_data = getattr(event, var)
                hist_data.Fill(val_data)

                t_val_data    = getattr(event, "t1")
                open_ang_data = getattr(event, "open_angle_ep2")
                if analysis_type == "dvcs":
                    theta_nn_data = getattr(event, "theta_gamma_gamma")
                else:
                    theta_nn_data = getattr(event, "theta_pi0_pi0")
                #endif

                Emiss2_d  = getattr(event, "Emiss2")
                Mx2_d     = getattr(event, "Mx2")
                Mx2_1_d   = getattr(event, "Mx2_1")
                Mx2_2_d   = getattr(event, "Mx2_2")
                pTmiss_d  = getattr(event, "pTmiss")
                xF_d      = getattr(event, "xF")

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

        # Fill MC histograms
        for event in mc_tree:
            det1_mc = getattr(event, "detector1")
            det2_mc = getattr(event, "detector2")
            if (topology == "(FD,FD)" and det1_mc == 1 and det2_mc == 1) or \
               (topology == "(CD,FD)" and det1_mc == 2 and det2_mc == 1) or \
               (topology == "(CD,FT)" and det1_mc == 2 and det2_mc == 0):
                val_mc = getattr(event, var)
                hist_mc.Fill(val_mc)

                t_val_mc    = getattr(event, "t1")
                open_ang_mc = getattr(event, "open_angle_ep2")
                if analysis_type == "dvcs":
                    theta_nn_mc = getattr(event, "theta_gamma_gamma")
                else:
                    theta_nn_mc = getattr(event, "theta_pi0_pi0")
                #endif

                Emiss2_m = getattr(event, "Emiss2")
                Mx2_m    = getattr(event, "Mx2")
                Mx2_1_m  = getattr(event, "Mx2_1")
                Mx2_2_m  = getattr(event, "Mx2_2")
                pTmiss_m = getattr(event, "pTmiss")
                xF_m     = getattr(event, "xF")

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

        # Normalize
        if hist_data.Integral() > 0:
            hist_data.Scale(1.0 / hist_data.Integral())
        #endif
        if hist_mc.Integral() > 0:
            hist_mc.Scale(1.0 / hist_mc.Integral())
        #endif
        if hist_data_loose.Integral() > 0:
            hist_data_loose.Scale(1.0 / hist_data_loose.Integral())
        #endif
        if hist_mc_loose.Integral() > 0:
            hist_mc_loose.Scale(1.0 / hist_mc_loose.Integral())
        #endif

        # Determine y-range
        y_max_no_cuts = 1.75 * max(hist_data.GetMaximum(), hist_mc.GetMaximum())
        y_max_loose   = 1.75 * max(hist_data_loose.GetMaximum(), hist_mc_loose.GetMaximum())

        # Draw "no cuts" histograms
        canvas.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)

        hist_data.SetTitle(plot_title)
        hist_data.SetXTitle(format_label_name(var, analysis_type))
        hist_data.SetYTitle("Normalized counts")
        hist_data.GetYaxis().SetRangeUser(0, y_max_no_cuts)
        hist_data.Draw("E1")
        hist_mc.Draw("E1 SAME")

        legend_nocuts = ROOT.TLegend(0.4, 0.7, 0.9, 0.9)
        legend_nocuts.AddEntry(hist_data, f"Data ({int(hist_data.GetEntries())} evts)", "lep")
        legend_nocuts.AddEntry(hist_mc,   f"MC ({int(hist_mc.GetEntries())} evts)",   "lep")
        legend_nocuts.SetTextSize(0.03)
        legend_nocuts.Draw()

        # Draw "loose cuts" histograms
        canvas_loose_cuts.cd(i + 1)
        ROOT.gPad.SetLeftMargin(0.15)
        ROOT.gPad.SetBottomMargin(0.15)

        hist_data_loose.SetTitle(plot_title)
        hist_data_loose.SetXTitle(format_label_name(var, analysis_type))
        hist_data_loose.SetYTitle("Normalized counts")
        hist_data_loose.GetYaxis().SetRangeUser(0, y_max_loose)
        hist_data_loose.Draw("E1")
        hist_mc_loose.Draw("E1 SAME")

        # Compute mean + stddev => define ±3σ
        mean_data = hist_data_loose.GetMean()
        std_data  = hist_data_loose.GetStdDev()
        mean_mc   = hist_mc_loose.GetMean()
        std_mc    = hist_mc_loose.GetStdDev()

        data_min_cut = mean_data - 3.0*std_data
        data_max_cut = mean_data + 3.0*std_data
        mc_min_cut   = mean_mc   - 3.0*std_mc
        mc_max_cut   = mean_mc   + 3.0*std_mc

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

        legend_loose = ROOT.TLegend(0.275, 0.6, 0.9, 0.9)
        legend_loose.AddEntry(hist_data_loose, 
                              f"Data (μ={mean_data:.3f}, σ={std_data:.3f})", "lep")
        legend_loose.AddEntry(hist_mc_loose, 
                              f"MC   (μ={mean_mc:.3f}, σ={std_mc:.3f})", "lep")
        legend_loose.SetTextSize(0.03)
        legend_loose.Draw()
    #endfor

    # Save canvases
    clean_title = plot_title.replace(" ", "_")
    outname_nocuts = os.path.join(final_output_dir, 
                                  f"exclusivity_plots_{clean_title}_{topology}_no_cuts.png")
    outname_cuts   = os.path.join(final_output_dir, 
                                  f"exclusivity_plots_{clean_title}_{topology}_kinematic_cuts.png")

    canvas.SaveAs(outname_nocuts)
    canvas_loose_cuts.SaveAs(outname_cuts)

    # Cleanup
    del canvas
    del canvas_loose_cuts

    return cuts_dictionary
#enddef


# --------------------------------------------------------------------------------------
# 5) process_combination
# --------------------------------------------------------------------------------------
def process_combination(analysis_type, run_period_key, topology, tree_dict, output_dir):
    """
    Runs determine_exclusivity_cuts for a given (analysis, run_period, topology).
    Returns (analysis_type, run_period_key, topology, cuts_dictionary).
    """

    # ---------------------------------------------------------
    #  FIX #2: Directly define run_period & plot_title from key
    # ---------------------------------------------------------
    # Instead of relying on searching "Fa18 Inb" in the string, 
    # we know the dictionary key is "Fa18_inb", "Fa18_out", or "Sp19_inb".
    # We'll build the "user-friendly" run_period string accordingly:
    if run_period_key == "Fa18_inb":
        run_period  = "RGA Fa18 Inb"
        period_text = "Fa18 Inb"
    elif run_period_key == "Fa18_out":
        run_period  = "RGA Fa18 Out"
        period_text = "Fa18 Out"
    elif run_period_key == "Sp19_inb":
        run_period  = "RGA Sp19 Inb"
        period_text = "Sp19 Inb"
    else:
        # fallback:
        run_period  = "Unknown"
        period_text = run_period_key
    #endif

    # Build a better plot title:
    if analysis_type == "dvcs":
        data_tree = tree_dict["dvcs_data"][run_period_key]
        mc_tree   = tree_dict["dvcs_mc"][run_period_key]
        plot_title = f"{period_text} DVCS"
    else:
        # eppi0
        data_tree = tree_dict["eppi0_mc_recon"][run_period_key]
        mc_tree   = tree_dict["eppi0_bkg"][run_period_key]
        plot_title = f"{period_text} eppi0"
    #endif

    cuts = determine_exclusivity_cuts(analysis_type=analysis_type,
                                      topology=topology,
                                      data_tree=data_tree,
                                      mc_tree=mc_tree,
                                      output_dir=output_dir,
                                      run_period=run_period,
                                      plot_title=plot_title)
    return (analysis_type, run_period_key, topology, cuts)
#enddef


# --------------------------------------------------------------------------------------
# 6) main
# --------------------------------------------------------------------------------------
def main():
    # We'll create a single "exclusivity_plots" folder
    output_dir = "exclusivity_plots"
    os.makedirs(output_dir, exist_ok=True)

    # Load all ROOT files
    tree_dict = load_root_files()

    topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    dvcs_period_keys   = ["Fa18_inb", "Fa18_out", "Sp19_inb"]
    eppi0_period_keys  = ["Fa18_inb", "Fa18_out", "Sp19_inb"]

    tasks = []
    for rp in dvcs_period_keys:
        for topo in topologies:
            tasks.append(("dvcs", rp, topo))
        #endfor
    #endfor
    for rp in eppi0_period_keys:
        for topo in topologies:
            tasks.append(("eppi0", rp, topo))
        #endfor
    #endfor

    combined_cuts = {}

    # Use a ProcessPoolExecutor for parallel runs, ignoring "already exists" errors
    # We do exist_ok=True on directory creation, so no crash even if multiple processes
    # attempt to create the same subdir. 
    with concurrent.futures.ProcessPoolExecutor() as executor:
        future_to_task = {
            executor.submit(process_combination, 
                            analysis_type=a, 
                            run_period_key=rp, 
                            topology=topo, 
                            tree_dict=tree_dict, 
                            output_dir=output_dir): (a, rp, topo)
            for (a, rp, topo) in tasks
        }
        #endfor

        for future in concurrent.futures.as_completed(future_to_task):
            (analysis, rp, topo) = future_to_task[future]
            try:
                result = future.result()
                # result is (analysis_type, run_period_key, topology, cuts)
                (_, _, _, cuts) = result
                combined_cuts[(analysis, rp, topo)] = cuts
            except Exception as exc:
                print(f"Task for {analysis}, {rp}, {topo} generated an exception: {exc}")
            #endif
        #endfor
    #endwith

    # Save combined results to JSON
    out_json = os.path.join(output_dir, "combined_exclusivity_cuts.json")
    jdict = {}
    for (analysis, rp, topo), val in combined_cuts.items():
        key_str = f"{analysis}__{rp}__{topo}"
        jdict[key_str] = val
    #endfor

    with open(out_json, "w") as f:
        json.dump(jdict, f, indent=4)
    #endwith

    print(f"All cuts saved to {out_json}")
#enddef


if __name__ == "__main__":
    main()
#endif