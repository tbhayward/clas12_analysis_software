# exclusivity.py

import os
import json
import ROOT

from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from label_formatting import format_label_name
from root_io import load_root_files

# 1) Multi-stage approach, single JSON at the end
def process_period_multi_stage(period, output_dir, analysis_type):
    """
    We have 3 stages, plus a final pass:
      stage_vars = [
         ["Mx2", "Mx2_1"],        # stage 0
         ["Emiss2", "Mx2_2"],     # stage 1
         ["pTmiss", "xF", <theta>]  # stage 2 (with fitted variables)
      ]
    The final pass is stage 3 => no new measurement, just final distribution.
    """
    print(f"⏳ Multi-stage processing for period='{period}' (analysis_type='{analysis_type}')...")

    period_code, trees = load_root_files(period)

    # Friendly names
    run_info_map = {
        "DVCS_Fa18_inb":  ("RGA Fa18 Inb DVCS",  "Fa18 Inb DVCS"),
        "DVCS_Fa18_out":  ("RGA Fa18 Out DVCS",  "Fa18 Out DVCS"),
        "DVCS_Sp19_inb":  ("RGA Sp19 Inb DVCS",  "Sp19 Inb DVCS"),
        "eppi0_Fa18_inb": ("RGA Fa18 Inb eppi0", "Fa18 Inb eppi0"),
        "eppi0_Fa18_out": ("RGA Fa18 Out eppi0", "Fa18 Out eppi0"),
        "eppi0_Sp19_inb": ("RGA Sp19 Inb eppi0", "Sp19 Inb eppi0")
    }
    run_info = run_info_map.get(period_code, (period_code, period_code))

    # 2) stage_vars
    # Stage 0 => Mx2, Mx2_1
    # Stage 1 => Emiss2, Mx2_2
    # Stage 2 => pTmiss, xF, plus the "theta_neutral_neutral" variable
    stage_vars = [
        ["Mx2", "Mx2_1"],
        ["Emiss2", "Mx2_2"],
        ["pTmiss", "xF"]  # We'll append the correct theta below
    ]
    # If dvcs => "theta_gamma_gamma", else => "theta_pi0_pi0"
    if analysis_type == "dvcs":
        stage_vars[-1].append("theta_gamma_gamma")
    else:
        stage_vars[-1].append("theta_pi0_pi0")
    #endif

    # 3) Loop topologies
    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        print(f"  🔄 Processing topology {topology}")

        # Single dictionary for final mu±3σ across all stages
        cumulative_cuts_dict = {"data": {}, "mc": {}}

        num_stages = len(stage_vars) + 1  # e.g. 3 + 1 = 4 total passes
        for stage_index in range(num_stages):
            # 1) Fill histograms
            data_hists, mc_hists = fill_stage_histograms(
                trees["data"], trees["mc"], topology, analysis_type,
                cumulative_cuts_dict, stage_index
            )

            # 2) Plot them
            cut_label = f"cut_{stage_index}"
            plot_title = f"{run_info[1]}"
            plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix=cut_label)

            # 3) If not the final pass, measure mu±3σ for stage_vars[stage_index]
            #    Now we do crystal-ball fits for pTmiss/theta to store in dictionary
            if stage_index < len(stage_vars):
                active_vars = stage_vars[stage_index]
                update_cuts_dict(data_hists, mc_hists, cumulative_cuts_dict, active_vars)
            #endif
        #endfor

        # 4) After final stage, save JSON once
        save_final_cuts(period_code, topology, output_dir, cumulative_cuts_dict)
    #endfor

    print(f"✅ Completed multi-stage for {period}\n")
#enddef


# 2) Fill
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
        return False

    # Fill data
    for event in data_tree:
        if not passes_topology(event):
            continue
        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, 0.0,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "data", "", topology
        ):
            continue
        if not passes_3sigma_cuts(event, False, cuts_dict):
            continue

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


# 3) get_hist_configs
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
            "Mx2_1":             (100, -1.5 , 1.5),
            "Mx2_2":             (100, 0, 3)
        }
    else:
        raise ValueError(f"Unrecognized analysis_type: {analysis_type}")
#enddef


# 4) update_cuts_dict
def update_cuts_dict(data_hists, mc_hists, cumulative_dict, active_vars):
    """
    For each var in active_vars, we *fit* a crystal ball for pTmiss & theta_...
    to get (mu, sigma), so the next stage's 3σ cut uses the fitted values.

    For other variables, we just store the raw histogram mean/std.
    """
    for var in active_vars:
        if var in ["theta_gamma_gamma", "theta_pi0_pi0", "pTmiss"]:
            # Fit crystal ball to get mu, sigma
            fit_mu_data, fit_sigma_data, _ = fit_crystal_ball(data_hists[var], var, True)
            fit_mu_mc,   fit_sigma_mc,   _ = fit_crystal_ball(mc_hists[var],   var, False)

            cumulative_dict["data"][var] = {"mean": fit_mu_data, "std": fit_sigma_data}
            cumulative_dict["mc"][var]   = {"mean": fit_mu_mc,   "std": fit_sigma_mc}
        else:
            # Use raw histogram
            d_mean = data_hists[var].GetMean()
            d_std  = data_hists[var].GetStdDev()
            m_mean = mc_hists[var].GetMean()
            m_std  = mc_hists[var].GetStdDev()

            cumulative_dict["data"][var] = {"mean": d_mean, "std": d_std}
            cumulative_dict["mc"][var]   = {"mean": m_mean, "std": m_std}
        #endif
    #endfor
#enddef


# 5) save_final_cuts
def save_final_cuts(period_code, topology, output_dir, cuts_dict):
    safe_topo = topology.replace("(", "").replace(")", "")
    filename = f"cuts_{period_code}_{safe_topo}_final.json"
    out_path = os.path.join(output_dir, filename)
    with open(out_path, "w") as f:
        json.dump(cuts_dict, f, indent=2)
    print(f"   ✅ Wrote final JSON => {out_path}")
#enddef


# 6) fit_crystal_ball (unchanged from before)
def fit_crystal_ball(hist, var_name, is_data=True):
    """
    Fit 'hist' with a 5-parameter Crystal Ball PDF.
    Returns (fit_mu, fit_sigma, fitFunc).
    """
    crystal_str = """
    function crystalBall(x, par) {
        double A    = par[0];
        double mean = par[1];
        double sigma= par[2];
        double alpha= par[3];
        double n    = par[4];
        double t = (x[0] - mean)/sigma;
        if (alpha < 0) t = -t;
        double absAlpha = fabs(alpha);

        if (t > -absAlpha) {
            return A * exp(-0.5 * t*t);
        } else {
            double nDivAbsA = n/absAlpha;
            double AA =  pow(nDivAbsA, n) * exp(-0.5*absAlpha*absAlpha);
            double B  = nDivAbsA - absAlpha;
            return A * AA * pow(B - t, -n);
        }
    }
    """
    ROOT.gInterpreter.ProcessLine(crystal_str)

    func_name = f"crystalBall_{var_name}_{'data' if is_data else 'mc'}"
    fcb = ROOT.TF1(func_name, "crystalBall", hist.GetXaxis().GetXmin(), hist.GetXaxis().GetXmax(), 5)

    A_guess = hist.GetMaximum()
    alpha_guess = 1.5
    n_guess = 2.0

    if var_name in ["theta_gamma_gamma", "theta_pi0_pi0"]:
        mu_guess    = 0.2
        sigma_guess = 0.02
    elif var_name == "pTmiss":
        mu_guess    = 0.05
        sigma_guess = 0.01
    else:
        mu_guess    = hist.GetMean()
        sigma_guess = hist.GetRMS() / 2.0
    #endif

    fcb.SetParameter(0, A_guess)
    fcb.SetParameter(1, mu_guess)
    fcb.SetParameter(2, sigma_guess)
    fcb.SetParameter(3, alpha_guess)
    fcb.SetParameter(4, n_guess)

    fit_opts = "R0Q"  # no drawing, quiet
    hist.Fit(fcb, fit_opts)

    fit_mu    = fcb.GetParameter(1)
    fit_sigma = fcb.GetParameter(2)

    return fit_mu, fit_sigma, fcb
#enddef


# 7) plot_results
def plot_results(data_hists, mc_hists, plot_title, topology, output_dir, suffix="cut_0"):
    """
    We do a crystal ball fit for pTmiss and theta_..., overlay it as a line,
    and in the legend show the fitted mu/sigma instead of raw means.
    """
    variables = list(data_hists.keys())
    canvas = ROOT.TCanvas("canvas", "", 2400, 1200)
    canvas.Divide(4, 2, 0.002, 0.002)

    base_legend = ROOT.TLegend(0.45, 0.75, 0.92, 0.89)
    base_legend.SetFillStyle(0)
    base_legend.SetTextFont(42)
    base_legend.SetBorderSize(1)
    base_legend.SetMargin(0.12)

    crystal_vars = ["theta_gamma_gamma", "theta_pi0_pi0", "pTmiss"]

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

        if dh.Integral() > 0:
            dh.Scale(1.0 / dh.Integral())
        if mh.Integral() > 0:
            mh.Scale(1.0 / mh.Integral())

        dh.GetYaxis().SetTitle("Normalized Counts")
        dh.GetYaxis().SetTitleOffset(1.4)
        dh.GetXaxis().SetTitle(format_label_name(var, "dvcs"))

        mh.GetYaxis().SetTitle("Normalized Counts")
        mh.GetYaxis().SetTitleOffset(1.4)
        mh.GetXaxis().SetTitle(format_label_name(var, "dvcs"))

        max_val = max(dh.GetMaximum(), mh.GetMaximum()) * 1.2
        dh.SetMaximum(max_val)
        mh.SetMaximum(max_val)

        if var in crystal_vars:
            # Fit data
            mu_data, sigma_data, fcb_data = fit_crystal_ball(dh, var, True)
            fcb_data.SetLineColor(ROOT.kBlue)
            fcb_data.SetLineStyle(1)
            fcb_data.SetLineWidth(2)

            # Fit MC
            mu_mc, sigma_mc, fcb_mc = fit_crystal_ball(mh, var, False)
            fcb_mc.SetLineColor(ROOT.kRed)
            fcb_mc.SetLineStyle(1)
            fcb_mc.SetLineWidth(2)

            dh.Draw("E1")
            mh.Draw("E1 SAME")
            fcb_data.Draw("SAME L")
            fcb_mc.Draw("SAME L")
        else:
            # Just raw hist means
            mu_data   = dh.GetMean()
            sigma_data= dh.GetStdDev()
            mu_mc     = mh.GetMean()
            sigma_mc  = mh.GetStdDev()

            dh.Draw("E1")
            mh.Draw("E1 SAME")
        #endif

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


# 8) combine_results (if you still want to combine final JSON outputs)
def combine_results(output_dir):
    import json

    combined = {}
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
            key_name = f"{period_code}_{topo}"
            try:
                with open(fpath, "r") as f:
                    combined[key_name] = json.load(f)
            except FileNotFoundError:
                print(f"⚠️ Missing file {fpath}")
        #endfor
    #endfor

    combined_path = os.path.join(output_dir, "combined_cuts.json")
    with open(combined_path, "w") as f:
        json.dump(combined, f, indent=2)
    print(f"✅ Wrote combined JSON => {combined_path}")
#enddef