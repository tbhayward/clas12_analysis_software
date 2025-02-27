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
    This function implements your real logic for 'apply_kinematic_cuts', 
    checking if an event passes the "loose" cut conditions. 
    Replace this placeholder with your actual cut thresholds/conditions 
    from your C++ code.
    """

    # ---------------------------------------------------------------------------
    # PLACEHOLDER LOGIC:
    # In your real code, you would do something like:
    #
    # if Emiss2_value < 0.3 and Mx2_value > -0.1 and Mx2_value < 1.4 and pTmiss_value < 0.4:
    #     return True
    # else:
    #     return False
    #
    # Or whatever your actual logic is. Currently we just return True.
    # ---------------------------------------------------------------------------
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
    # Simple placeholder logic:
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
    This function loads all your ROOT files into TChains keyed by 
    (category, run period). 
    No lines omitted -- every path you gave is included explicitly.
    """

    # A dictionary of all file paths
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
            "Sp19_inb": "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_sp19_inb_50nA_10604MeV_epgamma.root"
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

    # Build a dictionary of TChains
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
                               plot_title):
    """
    This function replicates your C++ determine_exclusivity(...) logic:
      - Creates canvases
      - Defines variables (depends on dvcs/eppi0)
      - Loops over events to fill histograms
      - Applies "loose" kinematic cuts (apply_kinematic_cuts)
      - Normalizes, draws, and saves histograms
      - Extracts mean±3σ to define cuts for each variable in data vs. MC
    Returns a dictionary of these cuts for subsequent usage.
    """

    # 1) Global style options
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetTitleSize(0.05, "XY")
    ROOT.gStyle.SetLabelSize(0.04, "XY")
    ROOT.gStyle.SetLegendTextSize(0.04)

    # 2) Create subdirectory: "dvcs" or "eppi0"
    analysis_dir = "dvcs" if analysis_type == "dvcs" else "eppi0"
    final_output_dir = os.path.join(output_dir, "exclusivity_plots", analysis_dir)
    if not os.path.exists(final_output_dir):
        os.makedirs(final_output_dir)
    #endif

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
        # eppi0
        variables = ["open_angle_ep2",
                     "theta_pi0_pi0",
                     "pTmiss",
                     "xF",
                     "Emiss2",
                     "Mx2",
                     "Mx2_1",
                     "Mx2_2"]
    #endif

    # 5) Parse run period from plot_title (like your C++ code)
    if "Fa18 Inb" in plot_title:
        run_period = "RGA Fa18 Inb"
    elif "Fa18 Out" in plot_title:
        run_period = "RGA Fa18 Out"
    elif "Sp19 Inb" in plot_title:
        run_period = "RGA Sp19 Inb"
    else:
        print(f"Unknown run period in plotTitle: {plot_title}")
        return {}
    #endif

    # 6) Prepare dictionary to store final cuts
    #    We'll store them like: 
    #    cuts_dictionary["data"][var] = { "mean":..., "std":..., "min_cut":..., "max_cut":... }
    #    cuts_dictionary["mc"]  [var] = { "mean":..., "std":..., "min_cut":..., "max_cut":... }
    cuts_dictionary = {
        "data": {},
        "mc":   {}
    }

    # 7) Hist binning (like histConfigs in your C++ code). 
    #    Customize these as needed.
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

        # Style adjustments
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

            # Match the chosen topology
            if (topology == "(FD,FD)" and det1 == 1 and det2 == 1) or \
               (topology == "(CD,FD)" and det1 == 2 and det2 == 1) or \
               (topology == "(CD,FT)" and det1 == 2 and det2 == 0):
                val_data = getattr(event, var)
                hist_data.Fill(val_data)

                # Now check if event passes "loose" cuts
                t_val_data      = getattr(event, "t1")
                open_ang_data   = getattr(event, "open_angle_ep2")
                if analysis_type == "dvcs":
                    theta_nn_data = getattr(event, "theta_gamma_gamma")
                else:
                    theta_nn_data = getattr(event, "theta_pi0_pi0")
                #endif

                Emiss2_d    = getattr(event, "Emiss2")
                Mx2_d       = getattr(event, "Mx2")
                Mx2_1_d     = getattr(event, "Mx2_1")
                Mx2_2_d     = getattr(event, "Mx2_2")
                pTmiss_d    = getattr(event, "pTmiss")
                xF_d        = getattr(event, "xF")

                # apply cuts
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

                t_val_mc      = getattr(event, "t1")
                open_ang_mc   = getattr(event, "open_angle_ep2")
                if analysis_type == "dvcs":
                    theta_nn_mc = getattr(event, "theta_gamma_gamma")
                else:
                    theta_nn_mc = getattr(event, "theta_pi0_pi0")
                #endif

                Emiss2_m   = getattr(event, "Emiss2")
                Mx2_m      = getattr(event, "Mx2")
                Mx2_1_m    = getattr(event, "Mx2_1")
                Mx2_2_m    = getattr(event, "Mx2_2")
                pTmiss_m   = getattr(event, "pTmiss")
                xF_m       = getattr(event, "xF")

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

    # 9) Save canvases
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
# 5) Example function to run one combination (analysis, run_period, topology) 
#    and store the cuts
# --------------------------------------------------------------------------------------
def process_combination(analysis_type, 
                        run_period_key,
                        topology,
                        tree_dict,
                        output_dir):
    """
    Load the data/MC TChain from tree_dict for the given analysis/run_period,
    call determine_exclusivity_cuts, and return the resulting dictionary of cut boundaries.
    """
    # figure out which categories in tree_dict to use:
    # For DVCS, we might assume "dvcs_data" / "dvcs_mc"
    # For eppi0, we might use "eppi0_mc_recon" or "eppi0_bkg"? 
    # Adjust to your logic or pass in which categories you want. 
    # As an example, let's assume we always do dvcs_data/dvcs_mc for dvcs 
    # and eppi0_mc_recon/eppi0_bkg for eppi0 (customize as needed).
    if analysis_type == "dvcs":
        data_tree = tree_dict["dvcs_data"][run_period_key]
        mc_tree   = tree_dict["dvcs_mc"][run_period_key]
        plot_title = f"{run_period_key} DVCS"
    else:
        # eppi0
        data_tree = tree_dict["eppi0_mc_recon"][run_period_key]  # example choice
        mc_tree   = tree_dict["eppi0_bkg"][run_period_key]       # example choice
        plot_title = f"{run_period_key} eppi0"
    #endif

    # call function
    cuts = determine_exclusivity_cuts(analysis_type,
                                      topology,
                                      data_tree,
                                      mc_tree,
                                      output_dir,
                                      plot_title)
    # Return the dictionary so we can combine with others
    return (analysis_type, run_period_key, topology, cuts)
#enddef


# --------------------------------------------------------------------------------------
# 6) main
# --------------------------------------------------------------------------------------
def main():
    # We assume you run this in the "rga_dvcs_BSA" directory 
    # so we can create "exclusivity_plots" subdirectory here.
    output_dir = "exclusivity_plots"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    #endif

    # Load all ROOT files
    tree_dict = load_root_files()

    # Example topologies you want to check
    topologies = ["(FD,FD)", "(CD,FD)", "(CD,FT)"]
    # Example run period keys for dvcs or eppi0
    dvcs_period_keys   = ["Fa18_inb", "Fa18_out", "Sp19_inb"]
    eppi0_period_keys  = ["Fa18_inb", "Fa18_out", "Sp19_inb"]

    # We'll set up a list of tasks (analysis_type, run_period_key, topology)
    tasks = []
    for rp in dvcs_period_keys:
        for topo in topologies:
            tasks.append(("dvcs", rp, topo))
        #endfor
    #endfor

    # Similarly, if you want to do eppi0:
    for rp in eppi0_period_keys:
        for topo in topologies:
            tasks.append(("eppi0", rp, topo))
        #endfor
    #endfor

    # We'll run them in parallel using concurrent.futures. 
    # Each task will produce a dictionary of cuts.
    combined_cuts = {}  # e.g., combined_cuts[(analysis, run_period, topology)] = dict_of_cuts

    # We can limit the number of workers or let it default to CPU count
    with concurrent.futures.ProcessPoolExecutor() as executor:
        future_to_task = {
            executor.submit(process_combination, 
                            analysis_type=analysis, 
                            run_period_key=rp, 
                            topology=topo, 
                            tree_dict=tree_dict, 
                            output_dir=output_dir): (analysis, rp, topo)
            for (analysis, rp, topo) in tasks
        }
        #endfor

        # Collect results
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

    # Now we have a dictionary 'combined_cuts' containing all the cut info 
    # for each (analysis, run_period, topology). 
    # We can save it to a JSON if you wish:
    out_json = os.path.join(output_dir, "combined_exclusivity_cuts.json")
    with open(out_json, "w") as f:
        # We might want to convert the keys to strings for JSON
        # e.g. "dvcs___Fa18_inb___(FD,FD)" : { "data": {...}, "mc": {...} }
        jdict = {}
        for (analysis, rp, topo), val in combined_cuts.items():
            key_str = f"{analysis}__{rp}__{topo}"
            jdict[key_str] = val
        #endfor
        json.dump(jdict, f, indent=4)
    #endwith

    print(f"All cuts saved to {out_json}")
#enddef


# --------------------------------------------------------------------------------------
# 7) If run directly, call main()
# --------------------------------------------------------------------------------------
if __name__ == "__main__":
    main()
#endif