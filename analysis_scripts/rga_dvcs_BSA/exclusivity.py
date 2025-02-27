# exclusivity.py

import os
import json
import ROOT

from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from label_formatting import format_label_name
from root_io import load_root_files

# -------------------------------------------------------------------------
# process_period_multi_stage
# -------------------------------------------------------------------------
def process_period_multi_stage(period, output_dir, analysis_type):
    """
    A multi-stage approach, but we only save one JSON at the end of all stages
    for each topology. We keep a "cumulative_cuts_dict" in memory, then at the
    final stage we have all discovered variables and their mu±3σ. We do:
        save_final_cuts(...)
    to write out exactly one file per topology.
    """
    print(f"⏳ Multi-stage processing for period='{period}' (analysis_type='{analysis_type}')...")

    period_code, trees = load_root_files(period)

    # Friendly naming for your output
    run_info_map = {
        "DVCS_Fa18_inb":  ("RGA Fa18 Inb DVCS",  "Fa18 Inb DVCS"),
        "DVCS_Fa18_out":  ("RGA Fa18 Out DVCS",  "Fa18 Out DVCS"),
        "DVCS_Sp19_inb":  ("RGA Sp19 Inb DVCS",  "Sp19 Inb DVCS"),
        "eppi0_Fa18_inb": ("RGA Fa18 Inb eppi0", "Fa18 Inb eppi0"),
        "eppi0_Fa18_out": ("RGA Fa18 Out eppi0", "Fa18 Out eppi0"),
        "eppi0_Sp19_inb": ("RGA Sp19 Inb eppi0", "Sp19 Inb eppi0")
    }
    run_info = run_info_map.get(period_code, (period_code, period_code))

    # Example sets of variables to measure/cut in multiple stages:
    # Stage 0 => Mx2, Mx2_1
    # Stage 1 => Emiss2, Mx2_2
    # Stage 2 => pTmiss, xF, plus "theta_gamma_gamma" or "theta_pi0_pi0"
    stage_vars = [
        ["Mx2", "Mx2_1"],
        ["Emiss2", "Mx2_2"],
        ["pTmiss", "xF"]
    ]
    # If eppi0, we replace gamma variable with "theta_pi0_pi0"
    # If dvcs, we use "theta_gamma_gamma"
    if analysis_type == "dvcs":
        stage_vars[-1].append("theta_gamma_gamma")
    else:
        stage_vars[-1].append("theta_pi0_pi0")
    #endfor

    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        print(f"  🔄 Processing topology {topology}")

        # This dictionary accumulates the final mu±3σ for all variables discovered
        # across all stages. We do NOT write partial JSON now.
        cumulative_cuts_dict = {"data": {}, "mc": {}}

        # We'll do len(stage_vars)+1 passes:
        #   pass 0 => measure stage_vars[0]
        #   pass 1 => measure stage_vars[1]
        #   pass 2 => measure stage_vars[2]
        #   pass 3 => final pass with all cuts in place, but no new measurement
        num_stages = len(stage_vars) + 1

        for stage_index in range(num_stages):
            # Fill histograms applying all previously known cuts
            data_hists, mc_hists = fill_stage_histograms(
                trees["data"], trees["mc"], topology, analysis_type,
                cumulative_cuts_dict, stage_index
            )

            # Plot them
            cut_label = f"cut_{stage_index}"
            plot_title = f"{run_info[1]}"
            plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix=cut_label)

            # If we are not at final stage, measure new mu/sigma for stage_vars[stage_index]
            if stage_index < len(stage_vars):
                active_vars = stage_vars[stage_index]
                update_cuts_dict(data_hists, mc_hists, cumulative_cuts_dict, active_vars)
            #endif
        #endfor

        # ⭐ At this point we've completed all stages for this topology,
        # and "cumulative_cuts_dict" has the final mu/sigma for everything.
        # We only save JSON *once* here:
        save_final_cuts(period_code, topology, output_dir, cumulative_cuts_dict)
    #endfor

    print(f"✅ Completed multi-stage for {period}\n")
#enddef

# -------------------------------------------------------------------------
# fill_stage_histograms
# -------------------------------------------------------------------------
def fill_stage_histograms(data_tree, mc_tree, topology, analysis_type,
                          cuts_dict, stage_index):
    """
    Fill histograms for data & MC, applying:
      1) The universal kinematic_cuts(...) for open_angle_ep2 > 5, etc.
      2) The existing 3sigma cuts from 'cuts_dict' for all previously discovered variables.
    'stage_index' is used only for naming histograms. We do NOT save partial JSON here.
    """
    hist_configs = get_hist_configs(analysis_type)

    data_hists = {}
    mc_hists   = {}

    for var, (nbins, xlow, xhigh) in hist_configs.items():
        data_name = f"data_{var}_stage{stage_index}"
        mc_name   = f"mc_{var}_stage{stage_index}"
        data_hists[var] = ROOT.TH1D(data_name, "", nbins, xlow, xhigh)
        mc_hists[var]   = ROOT.TH1D(mc_name,   "", nbins, xlow, xhigh)
    #endfor

    def passes_topology(e):
        if topology == "(FD,FD)":
            return (e.detector1 == 1 and e.detector2 == 1)
        elif topology == "(CD,FD)":
            return (e.detector1 == 2 and e.detector2 == 1)
        elif topology == "(CD,FT)":
            return (e.detector1 == 2 and e.detector2 == 0)
        else:
            return False

    # Fill data
    for event in data_tree:
        if not passes_topology(event):
            continue
        # universal cuts
        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, 0.0,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "data", "", topology
        ):
            continue
        # apply previous 3sigma cuts
        if not passes_3sigma_cuts(event, False, cuts_dict):
            continue
        # fill
        for var in hist_configs.keys():
            val = getattr(event, var, None)
            if val is not None:
                data_hists[var].Fill(val)
            #endif
        #endfor
    #endfor

    # Fill MC
    for event in mc_tree:
        if not passes_topology(event):
            continue
        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, 0.0,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "mc", "", topology
        ):
            continue
        if not passes_3sigma_cuts(event, True, cuts_dict):
            continue
        for var in hist_configs.keys():
            val = getattr(event, var, None)
            if val is not None:
                mc_hists[var].Fill(val)
            #endif
        #endfor
    #endfor

    return data_hists, mc_hists
#enddef

# -------------------------------------------------------------------------
# get_hist_configs
# -------------------------------------------------------------------------
def get_hist_configs(analysis_type):
    """
    Returns a dictionary of exactly the variables we want to plot for either dvcs or eppi0.
    That way, dvcs => 'theta_gamma_gamma', eppi0 => 'theta_pi0_pi0'.

    This ensures we only produce 8 subplots (4x2).
    """
    if analysis_type == "dvcs":
        # 8 variables total (e.g. open_angle_ep2, theta_gamma_gamma, pTmiss, xF, Emiss2, Mx2, Mx2_1, Mx2_2)
        return {
            "open_angle_ep2":    (100, 0, 60),
            "theta_gamma_gamma": (100, 0, 2),
            "pTmiss":            (100, 0, 0.3),
            "xF":                (100, -0.4, 0.2),
            "Emiss2":            (100, -1, 2),
            "Mx2":               (100, -0.025, 0.025),
            "Mx2_1":             (100, -1.5, 1.5),
            "Mx2_2":             (100, 0, 3)
        }
    elif analysis_type == "eppi0":
        # 8 variables total for eppi0 (replace theta_gamma_gamma w/ theta_pi0_pi0)
        return {
            "open_angle_ep2":    (100, 0, 60),
            "theta_pi0_pi0":     (100, 0, 2),
            "pTmiss":            (100, 0, 0.3),
            "xF":                (100, -0.4, 0.2),
            "Emiss2":            (100, -1, 2),
            "Mx2":               (100, -0.025, 0.025),
            "Mx2_1":             (100, -1.5 , 1.5),
            "Mx2_2":             (100, 0, 3)
        }
    else:
        raise ValueError(f"Unrecognized analysis_type: {analysis_type}")
    #enddef
#enddef

# -------------------------------------------------------------------------
# update_cuts_dict
# -------------------------------------------------------------------------
def update_cuts_dict(data_hists, mc_hists, cumulative_dict, active_vars):
    """
    For each var in active_vars, measure data_hists[var].mean/std, MC likewise,
    update the dictionary. We do not write any file here.
    """
    for var in active_vars:
        d_mean = data_hists[var].GetMean()
        d_std  = data_hists[var].GetStdDev()
        m_mean = mc_hists[var].GetMean()
        m_std  = mc_hists[var].GetStdDev()

        cumulative_dict["data"][var] = {"mean": d_mean, "std": d_std}
        cumulative_dict["mc"][var]   = {"mean": m_mean, "std": m_std}
    #endfor
#enddef

# -------------------------------------------------------------------------
# save_final_cuts
# -------------------------------------------------------------------------
def save_final_cuts(period_code, topology, output_dir, cuts_dict):
    """
    Called once after all multi-stage loops. Writes the final dictionary of mu/sigma
    for each variable we've discovered. Only a SINGLE JSON file is created for each (period_code, topology).
    """
    safe_topo = topology.replace("(", "").replace(")", "")
    filename = f"cuts_{period_code}_{safe_topo}_final.json"
    out_path = os.path.join(output_dir, filename)
    with open(out_path, "w") as f:
        json.dump(cuts_dict, f, indent=2)
    #endwith
    print(f"   ✅ Wrote final JSON => {out_path}")
#enddef

# -------------------------------------------------------------------------
# plot_results
# -------------------------------------------------------------------------
def plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix="cut_0"):
    """
    Compare data vs. MC histograms: normalize them, set max, display (mu, sigma).
    We append 'suffix' in the final PNG file name so we see e.g. _cut_0, _cut_1, etc.
    """
    variables = list(data_hists.keys())
    canvas = ROOT.TCanvas("canvas", "", 2400, 1200)
    canvas.Divide(4, 2, 0.002, 0.002)

    base_legend = ROOT.TLegend(0.45, 0.75, 0.92, 0.89)
    base_legend.SetFillStyle(0)
    base_legend.SetTextFont(42)
    base_legend.SetBorderSize(1)
    base_legend.SetMargin(0.12)

    for i, var in enumerate(variables):
        pad = canvas.cd(i + 1)
        pad.SetTicks(1, 1)

        dh = data_hists[var]
        mh = mc_hists[var]

        # Optional: Set a nice histogram title: "Fa18 Inb DVCS (FD,FD) cut_0"
        the_title = f"{plot_title} {topology} {suffix}"
        dh.SetTitle(the_title)
        mh.SetTitle(the_title)

        # Style
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

        # Normalize
        if dh.Integral() > 0:
            dh.Scale(1.0 / dh.Integral())
        #endif
        if mh.Integral() > 0:
            mh.Scale(1.0 / mh.Integral())
        #endif

        dh.GetYaxis().SetTitle("Normalized Counts")
        dh.GetYaxis().SetTitleOffset(1.4)
        dh.GetXaxis().SetTitle(format_label_name(var, "dvcs"))  # or analysis_type if you prefer

        mh.GetYaxis().SetTitle("Normalized Counts")
        mh.GetYaxis().SetTitleOffset(1.4)
        mh.GetXaxis().SetTitle(format_label_name(var, "dvcs"))

        max_val = max(dh.GetMaximum(), mh.GetMaximum()) * 1.2
        dh.SetMaximum(max_val)
        mh.SetMaximum(max_val)

        # Means & sigmas
        mu_data   = dh.GetMean()
        sigma_data = dh.GetStdDev()
        mu_mc     = mh.GetMean()
        sigma_mc  = mh.GetStdDev()

        dh.Draw("E1")
        mh.Draw("E1 SAME")

        pad_leg = base_legend.Clone()
        pad_leg.AddEntry(dh, f"Data (#mu={mu_data:.3f}, #sigma={sigma_data:.3f})", "lep")
        pad_leg.AddEntry(mh, f"MC (#mu={mu_mc:.3f}, #sigma={sigma_mc:.3f})", "lep")
        pad_leg.Draw()
    #endfor

    safe_topo = topology.replace("(", "").replace(")", "")
    clean_title = plot_title.replace(" ", "_")
    out_png = f"{clean_title}_{safe_topo}_{suffix}_comparison.png"
    out_png = out_png.replace(" ", "_")
    canvas.SaveAs(os.path.join(output_dir, out_png))
    del canvas
#enddef

def combine_results(output_dir):
    """
    Combine all final JSON files of the form:
      cuts_{period_code}_{topology}_final.json
    into a single combined_cuts.json for convenience.
    """
    import json, os

    combined = {}

    # Suppose your periods are:
    all_periods = [
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "eppi0_Fa18_inb", "eppi0_Fa18_out", "eppi0_Sp19_inb"
    ]
    topologies = ["FD_FD", "CD_FD", "CD_FT"]

    for period_code in all_periods:
        for topo in topologies:
            file_name = f"cuts_{period_code}_{topo}_final.json"
            fpath = os.path.join(output_dir, file_name)
            if not os.path.exists(fpath):
                continue
            #endif
            key_name = f"{period_code}_{topo}"
            try:
                with open(fpath, "r") as f:
                    combined[key_name] = json.load(f)
                #endwith
            except FileNotFoundError:
                print(f"⚠️ Missing file {fpath}")
            #endtry
        #endfor
    #endfor

    # Write out one big combined_cuts.json
    combined_path = os.path.join(output_dir, "combined_cuts.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    #endwith
    print(f"✅ Wrote combined JSON => {combined_path}")
#enddef