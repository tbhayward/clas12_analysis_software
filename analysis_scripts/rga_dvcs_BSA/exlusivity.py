# exclusivity.py

import os
import json
import ROOT

# 1) Import cut logic from kin_cuts.py
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts

# 2) Import label utility (format_label_name) from label_formatting.py
from label_formatting import format_label_name

# 3) Import the ROOT loader from root_io.py
from root_io import load_root_files

# --------------------------------------------------------------------------------------
# process_events (first pass)
# --------------------------------------------------------------------------------------
def process_events(tree, topology, analysis_type, is_mc):
    """
    First pass: fill histograms with the 'placeholder' kinematic cuts.
    Returns: (hists, hists_loose)
      - hists:        All events that pass 'topology' check
      - hists_loose:  Subset that also pass apply_kinematic_cuts(...)
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
    hists_loose = {}

    # Create histograms
    for var in variables:
        hname        = f"{'mc' if is_mc else 'data'}_{var}"
        hname_loose  = f"{'mc' if is_mc else 'data'}_loose_{var}"
        hists[var]        = ROOT.TH1D(hname,       "", *hist_configs[var])
        hists_loose[var]  = ROOT.TH1D(hname_loose, "", *hist_configs[var])
    #endfor

    # Define topology condition
    if topology == "(FD,FD)":
        det_cond = lambda d1, d2: d1 == 1 and d2 == 1
    elif topology == "(CD,FD)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 1
    elif topology == "(CD,FT)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 0
    else:
        # No valid topology recognized
        return hists, hists_loose
    #endif

    # Loop over events
    for event in tree:
        d1 = event.detector1
        d2 = event.detector2
        if not det_cond(d1, d2):
            continue
        #endif

        # Basic kinematic cut
        cuts_passed = apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type=analysis_type,
            data_type=("mc" if is_mc else "data"),
            run_period="",
            topology=topology
        )

        # Fill histograms
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
# process_events_with_json_cuts (second pass)
# --------------------------------------------------------------------------------------
def process_events_with_json_cuts(tree, topology, analysis_type, is_mc, cuts_dict):
    """
    Second pass: apply both the same placeholder kinematic cuts AND mu+/-3*sigma from the JSON dict.
    Returns: histograms with only events that pass both sets of cuts.
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

    for var in variables:
        hname = f"{'mc' if is_mc else 'data'}_3sigma_{var}"
        hists[var] = ROOT.TH1D(hname, "", *hist_configs[var])
    #endfor

    # Define topology condition
    if topology == "(FD,FD)":
        det_cond = lambda d1, d2: d1 == 1 and d2 == 1
    elif topology == "(CD,FD)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 1
    elif topology == "(CD,FT)":
        det_cond = lambda d1, d2: d1 == 2 and d2 == 0
    else:
        return hists
    #endif

    # Loop over events
    for event in tree:
        d1 = event.detector1
        d2 = event.detector2
        if not det_cond(d1, d2):
            continue
        #endif

        # Basic kinematic cuts first
        basic_cuts_passed = apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type=analysis_type,
            data_type=("mc" if is_mc else "data"),
            run_period="",
            topology=topology
        )
        if not basic_cuts_passed:
            continue
        #endif

        # Additional mu+/-3*sigma cut
        if not passes_3sigma_cuts(event, is_mc, cuts_dict):
            continue
        #endif

        # Fill histograms
        for var in variables:
            val = getattr(event, var)
            hists[var].Fill(val)
        #endfor
    #endfor

    return hists
#enddef

# --------------------------------------------------------------------------------------
# plot_results
# --------------------------------------------------------------------------------------
def plot_results(data_hists, mc_hists, plot_title, topology, output_dir):
    """
    Compare data vs. MC histograms: normalize, set max, draw, show (mu, sigma).
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
        for h in (dh, mh):
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())
            h.GetYaxis().SetTitle("Normalized Counts")
            h.GetYaxis().SetTitleOffset(1.4)
            h.GetXaxis().SetTitle(format_label_name(var, "dvcs"))
        #endfor

        # Adjust maximum
        max_val = max(dh.GetMaximum(), mh.GetMaximum()) * 1.2
        dh.SetMaximum(max_val)
        mh.SetMaximum(max_val)

        # Means & sigmas
        mu_data   = dh.GetMean()
        sigma_data = dh.GetStdDev()
        mu_mc     = mh.GetMean()
        sigma_mc  = mh.GetStdDev()

        # Draw
        dh.Draw("E1")
        mh.Draw("E1 SAME")

        pad_leg = base_legend.Clone()
        pad_leg.AddEntry(dh, f"Data (#mu={mu_data:.3f}, #sigma={sigma_data:.3f})", "lep")
        pad_leg.AddEntry(mh, f"MC (#mu={mu_mc:.3f}, #sigma={sigma_mc:.3f})", "lep")
        pad_leg.Draw()
    #endfor

    # Save the canvas
    clean_title = plot_title.replace(" ", "_").replace("(", "").replace(")", "")
    outname = f"{clean_title}_{topology}_comparison.png"
    canvas.SaveAs(os.path.join(output_dir, outname))
    del canvas
#enddef

# --------------------------------------------------------------------------------------
# save_cuts
# --------------------------------------------------------------------------------------
def save_cuts(period, topology, output_dir, data, mc):
    """
    Writes (mu, std) for each variable of data and MC to a JSON file in output_dir.
    """
    cuts = {
        "data": {
            var: {"mean": data[var].GetMean(), "std": data[var].GetStdDev()}
            for var in data
        },
        "mc": {
            var: {"mean": mc[var].GetMean(), "std": mc[var].GetStdDev()}
            for var in mc
        }
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
# process_period
# --------------------------------------------------------------------------------------
def process_period(period, output_dir):
    """
    1) Load trees for data & MC.
    2) First pass: fill histograms + plot.
    3) Save JSON cuts for data, MC.
    4) Second pass: use JSON to do mu +/- 3*sigma cuts, re-plot.
    """
    print(f"⏳ Starting processing for {period}...")

    period_code, trees = load_root_files(period)
    run_info = {
        "Fa18_inb": ("RGA Fa18 Inb", "Fa18 Inb"),
        "Fa18_out": ("RGA Fa18 Out", "Fa18 Out"),
        "Sp19_inb": ("RGA Sp19 Inb", "Sp19 Inb")
    }[period_code]

    for topology in ["(FD,FD)", "(CD,FD)", "(CD,FT)"]:
        print(f"  🔄 Processing topology {topology}")

        # (a) First pass
        data_hists, data_loose = process_events(trees["data"], topology, "dvcs", False)
        mc_hists,   mc_loose   = process_events(trees["mc"],   topology, "dvcs", True)

        # (b) Plot first-pass results
        plot_title = f"{run_info[1]} DVCS"
        plot_results(data_hists, mc_hists, plot_title, topology, output_dir)

        # (c) Save JSON cuts using the "loose" histograms
        save_cuts(period_code, topology, output_dir, data_loose, mc_loose)

        # (d) Second pass: read JSON and apply 3σ cuts
        json_path = os.path.join(output_dir, f"cuts_{period_code}_{topology}.json")
        if not os.path.exists(json_path):
            print(f"   ⚠️ JSON cuts file not found, skipping second pass: {json_path}")
            continue
        #endif
        with open(json_path, "r") as infile:
            cuts_dict = json.load(infile)
        #endwith

        data_hists_3s = process_events_with_json_cuts(
            trees["data"], topology, "dvcs", False, cuts_dict
        )
        mc_hists_3s   = process_events_with_json_cuts(
            trees["mc"],   topology, "dvcs", True,  cuts_dict
        )

        # (e) Plot second-pass results
        plot_title_3s = f"{run_info[1]} DVCS (3sigma Cuts)"
        plot_results(data_hists_3s, mc_hists_3s, plot_title_3s, topology, output_dir)
    #endfor

    print(f"✅ Completed {period}\n")
#enddef

# --------------------------------------------------------------------------------------
# combine_results
# --------------------------------------------------------------------------------------
def combine_results(output_dir):
    """
    Combine all 'cuts_{period}_{topology}.json' into one 'combined_cuts.json'.
    """
    combined = {}
    for period in ["Fa18_inb", "Fa18_out", "Sp19_inb"]:
        for topology in ["FD_FD", "CD_FD", "CD_FT"]:
            path_ = os.path.join(output_dir, f"cuts_{period}_{topology}.json")
            try:
                with open(path_, "r") as f:
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