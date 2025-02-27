#!/usr/bin/env python3

import os
import json
import ROOT

from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from label_formatting import format_label_name
from root_io import load_root_files

def process_period_multi_stage(period, output_dir, analysis_type):
    """
    We have 3 stages, plus a final pass:
      stage_vars = [
         ["Mx2", "Mx2_1"],        # stage 0
         ["Emiss2", "Mx2_2"],     # stage 1
         ["pTmiss", "xF", <theta> ]  # stage 2
      ]
    The final pass is stage 3 => no new measurement, just final distribution.
    """
    print(f"⏳ Multi-stage processing for period='{period}' (analysis_type='{analysis_type}')...")

    period_code, trees = load_root_files(period)

    run_info_map = {
        "DVCS_Fa18_inb":  ("RGA Fa18 Inb DVCS",  "Fa18 Inb DVCS"),
        "DVCS_Fa18_out":  ("RGA Fa18 Out DVCS",  "Fa18 Out DVCS"),
        "DVCS_Sp19_inb":  ("RGA Sp19 Inb DVCS",  "Sp19 Inb DVCS"),
        "eppi0_Fa18_inb": ("RGA Fa18 Inb eppi0", "Fa18 Inb eppi0"),
        "eppi0_Fa18_out": ("RGA Fa18 Out eppi0", "Fa18 Out eppi0"),
        "eppi0_Sp19_inb": ("RGA Sp19 Inb eppi0", "Sp19 Inb eppi0")
    }
    run_info = run_info_map.get(period_code, (period_code, period_code))

    # Define the stages for each pass
    stage_vars = [
        ["Mx2", "Mx2_1"],      # Stage 0
        ["Emiss2", "Mx2_2"],   # Stage 1
        ["pTmiss", "xF"]       # Stage 2 => we append the correct theta below
    ]
    if analysis_type == "dvcs":
        stage_vars[-1].append("theta_gamma_gamma")
    else:
        stage_vars[-1].append("theta_pi0_pi0")
    #endif

    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        print(f"  🔄 Processing topology {topology}")

        # A single dictionary for final mu±3σ across all stages
        cumulative_cuts_dict = {"data": {}, "mc": {}}

        num_stages = len(stage_vars) + 1  # e.g. 3 + 1 = 4 total passes
        for stage_index in range(num_stages):
            data_hists, mc_hists = fill_stage_histograms(
                trees["data"], trees["mc"], topology, analysis_type,
                cumulative_cuts_dict, stage_index
            )
            #endfor

            cut_label = f"cut_{stage_index}"
            plot_title = f"{run_info[1]}"
            plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix=cut_label)

            # Update the 3σ cut dictionary for the next stage
            if stage_index < len(stage_vars):
                active_vars = stage_vars[stage_index]
                update_cuts_dict(data_hists, mc_hists, cumulative_cuts_dict, active_vars)
            #endif
        #endfor

        save_final_cuts(period_code, topology, output_dir, cumulative_cuts_dict)
    #endfor

    print(f"✅ Completed multi-stage for {period}\n")
#enddef

def fill_stage_histograms(data_tree, mc_tree, topology, analysis_type, cuts_dict, stage_index):
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
        #endif
        return False
    #enddef

    # Fill data
    for event in data_tree:
        if not passes_topology(event):
            continue
        #endif

        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, 0.0,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "data", "", topology
        ):
            continue
        #endif

        if not passes_3sigma_cuts(event, False, cuts_dict):
            continue
        #endif

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
        #endif

        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, 0.0,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "mc", "", topology
        ):
            continue
        #endif

        if not passes_3sigma_cuts(event, True, cuts_dict):
            continue
        #endif

        for var in hist_configs.keys():
            val = getattr(event, var, None)
            if val is not None:
                mc_hists[var].Fill(val)
            #endif
        #endfor
    #endfor

    return data_hists, mc_hists
#enddef

def get_hist_configs(analysis_type):
    if analysis_type == "dvcs":
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
        return {
            "open_angle_ep2":    (100, 0, 60),
            "theta_pi0_pi0":     (100, 0, 2),
            "pTmiss":            (100, 0, 0.3),
            "xF":                (100, -0.4, 0.2),
            "Emiss2":            (100, -1, 2),
            "Mx2":               (100, -0.025, 0.025),
            "Mx2_1":             (100, -1.5, 1.5),
            "Mx2_2":             (100, 0, 3)
        }
    else:
        raise ValueError(f"Unrecognized analysis_type: {analysis_type}")
    #endif
#enddef

def update_cuts_dict(data_hists, mc_hists, cumulative_dict, active_vars):
    """
    For each var in active_vars, we do a left-side Gaussian fit for pTmiss & theta_* to get (mu, sigma).
    That way, the next stage's 3σ cut uses these fitted values. For other variables,
    we just store the raw histogram mean/std.
    """
    for var in active_vars:
        if var in ["theta_gamma_gamma", "theta_pi0_pi0", "pTmiss"]:
            # Use left-side Gaussian to get (mu, sigma)
            fit_mu_data, fit_sigma_data = fit_gaussian_left_side(data_hists[var], var, True)
            fit_mu_mc,   fit_sigma_mc   = fit_gaussian_left_side(mc_hists[var],   var, False)

            cumulative_dict["data"][var] = {"mean": fit_mu_data, "std": fit_sigma_data}
            cumulative_dict["mc"][var]   = {"mean": fit_mu_mc,   "std": fit_sigma_mc}
        else:
            d_mean = data_hists[var].GetMean()
            d_std  = data_hists[var].GetStdDev()
            m_mean = mc_hists[var].GetMean()
            m_std  = mc_hists[var].GetStdDev()

            cumulative_dict["data"][var] = {"mean": d_mean, "std": d_std}
            cumulative_dict["mc"][var]   = {"mean": m_mean, "std": m_std}
        #endif
    #endfor
#enddef

def save_final_cuts(period_code, topology, output_dir, cuts_dict):
    safe_topo = topology.replace("(", "").replace(")", "")
    filename = f"cuts_{period_code}_{safe_topo}_final.json"
    out_path = os.path.join(output_dir, filename)
    with open(out_path, "w") as f:
        json.dump(cuts_dict, f, indent=2)
    #endif
    print(f"   ✅ Wrote final JSON => {out_path}")
#enddef

def fit_gaussian_left_side(hist, var_name, is_data=True, return_tf1=False):
    """
    Fits a simple Gaussian to the LEFT side of the distribution:
      - from x_min (0) up to ~some fraction of the peak position.
    Returns (fit_mu, fit_sigma) or (fit_mu, fit_sigma, TF1) if return_tf1=True.
    """
    # 1) Decide the left boundary. Often 0 is a natural left edge for these variables.
    x_left = 0
    # 2) Find the bin with the maximum content:
    peak_bin = hist.GetMaximumBin()
    x_peak   = hist.GetBinCenter(peak_bin)

    # 3) Build a standard Gaussian TF1 from x_left to slightly passed the x_peak (avoid big tails)
    func_name = f"gausLeft_{var_name}_{'data' if is_data else 'mc'}"
    fgaus = ROOT.TF1(func_name, "gaus(0)", x_left, 1.2*x_peak)

    # 4) Initial parameter guesses:
    amp_guess  = hist.GetMaximum()
    mean_guess = hist.GetMean()
    sigma_guess= hist.GetRMS() / 2.0

    fgaus.SetParameter(0, amp_guess)
    fgaus.SetParameter(1, mean_guess)
    fgaus.SetParameter(2, sigma_guess)

    # Keep sigma > 0
    fgaus.SetParLimits(2, 1e-6, 5.0)

    # 5) Fit only in [x_left, 0.85*x_peak], quietly
    fit_opts = "R0Q"
    hist.Fit(fgaus, fit_opts)

    # 6) Extract results
    fit_mu    = fgaus.GetParameter(1)
    fit_sigma = fgaus.GetParameter(2)

    if return_tf1:
        return (fit_mu, fit_sigma, fgaus)
    else:
        return (fit_mu, fit_sigma)
    #endif
#enddef

def plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix="cut_0"):
    variables = list(data_hists.keys())
    canvas = ROOT.TCanvas("canvas", "", 2400, 1200)
    canvas.Divide(4, 2, 0.002, 0.002)

    base_legend = ROOT.TLegend(0.45, 0.75, 0.92, 0.89)
    base_legend.SetFillStyle(0)
    base_legend.SetTextFont(42)
    base_legend.SetBorderSize(1)
    base_legend.SetMargin(0.12)

    # We'll do left-side Gaussian fits for these variables:
    gaussian_vars = ["theta_gamma_gamma", "theta_pi0_pi0", "pTmiss"]

    for i, var in enumerate(variables):
        pad = canvas.cd(i + 1)
        pad.SetTicks(1, 1)

        dh = data_hists[var]
        mh = mc_hists[var]

        the_title = f"{plot_title} {topology} {suffix}"
        dh.SetTitle(the_title)
        mh.SetTitle(the_title)

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
        dh.GetXaxis().SetTitle(format_label_name(var, "dvcs"))
        mh.GetYaxis().SetTitle("Normalized Counts")
        mh.GetXaxis().SetTitle(format_label_name(var, "dvcs"))

        max_val = max(dh.GetMaximum(), mh.GetMaximum()) * 1.2
        dh.SetMaximum(max_val)
        mh.SetMaximum(max_val)

        if var in gaussian_vars:
            # Use the left-side Gaussian fit for both data and MC
            mu_data, sigma_data, fgaus_data = fit_gaussian_left_side(dh, var, True, return_tf1=True)
            mu_mc, sigma_mc, fgaus_mc       = fit_gaussian_left_side(mh, var, False, return_tf1=True)

            dh.Draw("E1")
            mh.Draw("E1 SAME")

            # Draw the fitted Gaussian over the histogram
            fgaus_data.SetNpx(1000)
            fgaus_data.SetLineColor(ROOT.kBlue + 1)
            fgaus_data.SetLineStyle(2)   # dashed
            fgaus_data.SetLineWidth(2)
            fgaus_data.Draw("SAME")

            fgaus_mc.SetNpx(1000)
            fgaus_mc.SetLineColor(ROOT.kRed + 1)
            fgaus_mc.SetLineStyle(2)    # dashed
            fgaus_mc.SetLineWidth(2)
            fgaus_mc.Draw("SAME")
        else:
            # If not in gaussian_vars, just draw the hist
            mu_data    = dh.GetMean()
            sigma_data = dh.GetStdDev()
            mu_mc      = mh.GetMean()
            sigma_mc   = mh.GetStdDev()

            dh.Draw("E1")
            mh.Draw("E1 SAME")
        #endif

        # Legend
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
    Optional final JSON combiner
    """
    import json

    combined = {}
    all_periods = [
        "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb",
        "eppi0_Fa18_inb", "eppi0_Fa18_out", "eppi0_Sp19_inb"
    ]
    topologies = ["FD_FD", "CD_FD", "CD_FT"]

    for period_code in all_periods:
        for topo in topologies:
            fname = f"cuts_{period_code}_{topo}_final.json"
            fpath = os.path.join(output_dir, fname)
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

    combined_path = os.path.join(output_dir, "combined_cuts.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    #endif
    print(f"✅ Wrote combined JSON => {combined_path}")
#enddef