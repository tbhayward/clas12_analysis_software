# exclusivity.py

import os
import json
import ROOT

# 1) Import cut logic
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts

# 2) Label function
from label_formatting import format_label_name

# 3) Root-file loader
from root_io import load_root_files

# -------------------------------------------------------------------------
# process_period_multi_stage
# -------------------------------------------------------------------------
def process_period_multi_stage(period, output_dir, analysis_type):
    """
    A multi-stage analysis approach for either 'dvcs' or 'eppi0'.
    We do repeated passes (stage 0..N), each time:
      - filling histograms,
      - plotting them,
      - measuring mu±3sigma for 'stage variables',
      - saving those cuts to a cumulative dictionary,
      - applying them in the next stage.

    The universal cut on open_angle_ep2 > 5 is enforced *inside*
    apply_kinematic_cuts(...), from the very first pass (stage 0).
    """
    print(f"⏳ Multi-stage processing for period='{period}' (analysis_type='{analysis_type}')...")

    # 1) Load data & MC trees
    period_code, trees = load_root_files(period)

    # Some user-friendly naming for output plots:
    run_info_map = {
        "DVCS_Fa18_inb":  ("RGA Fa18 Inb DVCS",  "Fa18 Inb DVCS"),
        "DVCS_Fa18_out":  ("RGA Fa18 Out DVCS",  "Fa18 Out DVCS"),
        "DVCS_Sp19_inb":  ("RGA Sp19 Inb DVCS",  "Sp19 Inb DVCS"),
        "eppi0_Fa18_inb": ("RGA Fa18 Inb eppi0", "Fa18 Inb eppi0"),
        "eppi0_Fa18_out": ("RGA Fa18 Out eppi0", "Fa18 Out eppi0"),
        "eppi0_Sp19_inb": ("RGA Sp19 Inb eppi0", "Sp19 Inb eppi0")
    }
    run_info = run_info_map.get(period_code, (period_code, period_code))

    # 2) Define stage variables we will measure & cut on incrementally.
    #    e.g. first stage: Mx2, Mx2_1; second stage: Emiss2, Mx2_2; final stage: pTmiss, xF, etc.
    #    Also note if analysis_type == "eppi0", we might want to use "theta_pi0_pi0" instead of "theta_gamma_gamma".
    stage_vars = [
        ["Mx2", "Mx2_1"],  # Stage 0 variables to measure => then cut
        ["Emiss2", "Mx2_2"],  # Stage 1
        ["theta_gamma_gamma", "pTmiss", "xF"]  # Stage 2
    ]
    if analysis_type == "eppi0":
        # Replace 'theta_gamma_gamma' with 'theta_pi0_pi0' in the final stage
        stage_vars[2][0] = "theta_pi0_pi0"
    #endif

    # 3) We'll handle multiple topologies:
    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        print(f"  🔄 Processing topology {topology}")

        # We'll keep a dictionary of all mu/std so far. Each stage updates it:
        cumulative_cuts_dict = {"data": {}, "mc": {}}

        # We'll do len(stage_vars)+1 passes => e.g. if stage_vars has 3 entries, we do passes 0..3
        # Stage 0 => measure stage_vars[0]
        # Stage 1 => measure stage_vars[1]
        # Stage 2 => measure stage_vars[2]
        # Stage 3 => final pass (no new measurement, just final distribution with all cuts applied)
        num_stages = len(stage_vars) + 1

        for stage_index in range(num_stages):
            # 3.1) Fill histograms applying all known cuts from previous stages
            data_hists, mc_hists = fill_stage_histograms(
                trees["data"], trees["mc"], topology, analysis_type,
                cumulative_cuts_dict, stage_index
            )

            # 3.2) Plot them with suffix => e.g. _cut_0, _cut_1, ...
            cut_label = f"cut_{stage_index}"
            plot_title = f"{run_info[1]}"
            plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix=cut_label)

            # 3.3) If we're not at the final stage, measure new mu/sigma for stage_vars[stage_index]
            if stage_index < len(stage_vars):
                active_vars = stage_vars[stage_index]
                # Update the cumulative dictionary with newly measured mu/sigmas
                update_cuts_dict(data_hists, mc_hists, cumulative_cuts_dict, active_vars)

                # Write them out to JSON for debugging or for post-check
                save_stage_cuts(period_code, topology, output_dir, cumulative_cuts_dict, stage_index)
            #endif
        #endfor
    #endfor

    print(f"✅ Completed multi-stage for {period}\n")
#enddef


def get_hist_configs(analysis_type):
    """
    Returns a dictionary mapping variable_name -> (nbins, xlow, xhigh).
    For DVCS: we include theta_gamma_gamma, but NOT theta_pi0_pi0.
    For eppi0: we include theta_pi0_pi0, but NOT theta_gamma_gamma.

    We also include Mx2_2 so it appears in the final plots.
    """
    if analysis_type == "dvcs":
        # 8 variables total (e.g. open_angle_ep2, theta_gamma_gamma, pTmiss, xF, Emiss2, Mx2, Mx2_1, Mx2_2)
        return {
            "open_angle_ep2":    (100, 0, 40),
            "theta_gamma_gamma": (100, 0, 2),
            "pTmiss":            (100, 0, 0.3),
            "xF":                (100, -0.4, 0.2),
            "Emiss2":            (100, -1, 2),
            "Mx2":               (100, -0.015, 0.015),
            "Mx2_1":             (100, -1, 1.5),
            "Mx2_2":             (100, 0, 3)
        }
    elif analysis_type == "eppi0":
        # 8 variables total for eppi0 (replace theta_gamma_gamma w/ theta_pi0_pi0)
        return {
            "open_angle_ep2":    (100, 0, 40),
            "theta_pi0_pi0":     (100, 0, 40),
            "pTmiss":            (100, 0, 0.3),
            "xF":                (100, -0.4, 0.2),
            "Emiss2":            (100, -1, 2),
            "Mx2":               (100, -0.015, 0.015),
            "Mx2_1":             (100, -1, 1.5),
            "Mx2_2":             (100, 0, 3)
        }
    else:
        raise ValueError(f"Unrecognized analysis_type: {analysis_type}")
    #enddef
#enddef

# -------------------------------------------------------------------------
# fill_stage_histograms
# -------------------------------------------------------------------------
def fill_stage_histograms(data_tree, mc_tree, topology, analysis_type,
                          cuts_dict, stage_index):
    """
    Fill histograms for data & MC, applying:
      1) The original apply_kinematic_cuts(...) for open_angle_ep2>5
      2) The 3sigma cuts from 'cuts_dict' for all variables discovered so far.

    'stage_index' is just used for naming the histograms (like _stage0).
    """
    # Get the correct dictionary for dvcs or eppi0
    hist_configs = get_hist_configs(analysis_type)

    # Decide which "theta" variable to fill if analysis_type == "dvcs" or "eppi0"
    # But we already handle that in 'root_io.py' by only loading the correct branch name.
    # So we can fill them all if they exist or not.

    vars_to_fill = list(hist_configs.keys())

    # Create the histograms
    data_hists = {}
    mc_hists   = {}

    for var in vars_to_fill:
        data_hists[var] = ROOT.TH1D(f"data_{var}_stage{stage_index}", "", *hist_configs[var])
        mc_hists[var]   = ROOT.TH1D(f"mc_{var}_stage{stage_index}",   "", *hist_configs[var])
    #endfor

    # Make a topology condition function
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
        # Apply our "placeholder" kinematic cuts, including open_angle_ep2>5
        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, getattr(event, "theta_gamma_gamma", 0), 
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, analysis_type, "data", "", topology
        ):
            continue
        # Now apply any discovered 3sigma cuts
        if not passes_3sigma_cuts(event, False, cuts_dict):
            continue
        # If we get here, the event passes everything
        for var in vars_to_fill:
            val = getattr(event, var, None)
            if val is not None:
                data_hists[var].Fill(val)
    #endfor

    # Fill MC
    for event in mc_tree:
        if not passes_topology(event):
            continue
        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, getattr(event, "theta_gamma_gamma", 0), 
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, analysis_type, "mc", "", topology
        ):
            continue
        if not passes_3sigma_cuts(event, True, cuts_dict):
            continue
        # If we get here, the event is accepted
        for var in vars_to_fill:
            val = getattr(event, var, None)
            if val is not None:
                mc_hists[var].Fill(val)
    #endfor

    return data_hists, mc_hists
#enddef


# -------------------------------------------------------------------------
# update_cuts_dict
# -------------------------------------------------------------------------
def update_cuts_dict(data_hists, mc_hists, cumulative_dict, active_vars):
    """
    For the given 'active_vars', we compute data_hists[var].GetMean/StdDev
    and MC hists. Then we store them in 'cumulative_dict', which is used
    by passes_3sigma_cuts(...) to enforce mu+/-3sigma in the subsequent stage.
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
# save_stage_cuts
# -------------------------------------------------------------------------
def save_stage_cuts(period_code, topology, output_dir, cuts_dict, stage_index):
    """
    Write out the current dictionary of mu±3sigma for the variables
    discovered so far. e.g. cuts_DVCS_Fa18_inb_(FD,FD)_stage0.json
    """
    safe_topology = topology.replace("(", "").replace(")", "")
    filename = f"cuts_{period_code}_{safe_topology}_stage{stage_index}.json"
    out_path = os.path.join(output_dir, filename)
    with open(out_path, "w") as f:
        json.dump(cuts_dict, f, indent=2)
    print(f"   📝 Saved stage-{stage_index} JSON => {out_path}")
#enddef


# -------------------------------------------------------------------------
# plot_results
# -------------------------------------------------------------------------
def plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix="cut_0"):
    """
    Compare data vs. MC histograms: normalize them, set max, display (mu, sigma).
    We also append 'suffix' to the final PNG file name, e.g. _cut_0, _cut_1, etc.

    'plot_title' typically is something like "Fa18 Inb DVCS" from run_info[1].
    We'll incorporate 'topology' and 'suffix' in the final filename.
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

        # Optionally set the histogram title to show period, topology, suffix:
        # e.g. "Fa18 Inb DVCS (FD,FD) cut_1"
        hist_title = f"{plot_title} {topology} {suffix}"
        dh.SetTitle(hist_title)
        mh.SetTitle(hist_title)

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
        for h in [dh, mh]:
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())
            h.GetYaxis().SetTitle("Normalized Counts")
            h.GetYaxis().SetTitleOffset(1.4)
            # We pass analysis_type='dvcs' or 'eppi0' if we want different labels
            # But let's just pass "dvcs" for demonstration or adapt as needed
            h.GetXaxis().SetTitle(format_label_name(var, "dvcs"))
        #endfor

        # Set max
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

        # Legend
        pad_leg = base_legend.Clone()
        pad_leg.AddEntry(dh, f"Data (#mu={mu_data:.3f}, #sigma={sigma_data:.3f})", "lep")
        pad_leg.AddEntry(mh, f"MC (#mu={mu_mc:.3f}, #sigma={sigma_mc:.3f})", "lep")
        pad_leg.Draw()
    #endfor

    # Save
    safe_topology = topology.replace("(", "").replace(")", "")
    clean_title = plot_title.replace(" ", "_")
    out_png = f"{clean_title}_{safe_topology}_{suffix}_comparison.png"
    out_png = out_png.replace(" ", "_")
    canvas.SaveAs(os.path.join(output_dir, out_png))
    del canvas
#enddef


# -------------------------------------------------------------------------
# combine_results
# -------------------------------------------------------------------------
def combine_results(output_dir):
    """
    Combine all stage cut JSONs (cuts_period_topology_stageX.json) 
    into one big 'combined_cuts.json' if desired.
    This is called once at the end from main.py
    """
    combined = {}

    # We'll guess your periods might be these:
    all_periods = [
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "eppi0_Fa18_inb", "eppi0_Fa18_out", "eppi0_Sp19_inb"
    ]
    topologies = ["FD_FD", "CD_FD", "CD_FT"]

    # We'll look for stage0..stageN
    for period_code in all_periods:
        for topo in topologies:
            stage_index = 0
            while True:
                fname = f"cuts_{period_code}_{topo}_stage{stage_index}.json"
                fpath = os.path.join(output_dir, fname)
                if not os.path.exists(fpath):
                    break
                #endif
                key_name = f"{period_code}_{topo}_stage{stage_index}"
                try:
                    with open(fpath, "r") as f:
                        combined[key_name] = json.load(f)
                    #endwith
                except FileNotFoundError:
                    print(f"⚠️ Missing file {fpath}")
                #endif
                stage_index += 1
            #endwhile
        #endfor
    #endfor

    combined_path = os.path.join(output_dir, "combined_cuts.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"✅ Wrote combined JSON => {combined_path}")
#enddef