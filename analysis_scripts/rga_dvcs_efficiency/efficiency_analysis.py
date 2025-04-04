#!/usr/bin/env python3
import os
import glob
import json
import math
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from load_binning_scheme import load_binning_scheme

#############################################
# Integrated Efficiency Functions
#############################################

def load_cuts(period, topology):
    """
    Loads the final cuts dictionary from combined_cuts.json for a given (period, topology).
    Expect keys such as "DVCS_Fa18_inb_FD_FD", etc.
    """
    topo_clean = topology.replace("(", "").replace(")", "").replace(",", "_").strip()
    if "bkg" in period:
        dvcs_equiv = period.replace("eppi0_bkg", "DVCS")
        dictionary_key = f"{dvcs_equiv}_{topo_clean}"
    else:
        dictionary_key = f"{period}_{topo_clean}"
    combined_cuts_path = os.path.join("exclusivity", "combined_cuts.json")
    if not os.path.exists(combined_cuts_path):
        print(f"⚠️ {combined_cuts_path} does not exist; returning empty cuts dictionary.")
        return {}
    with open(combined_cuts_path, "r") as f:
        combined_cuts = json.load(f)
    if dictionary_key in combined_cuts:
        return combined_cuts[dictionary_key]
    else:
        print(f"⚠️ Key '{dictionary_key}' not found in {combined_cuts_path}.")
        print("Available keys:", list(combined_cuts.keys()))
        return {}
#enddef

def passes_exclusivity_cuts(event, cuts, topology, analysis_type, run_period):
    """
    Applies exclusivity cuts (for integrated analysis) to an event given a fixed topology.
    """
    if event.detector1 == 1 and event.detector2 == 1:
        event_topo = "(FD,FD)"
    elif event.detector1 == 2 and event.detector2 == 1:
        event_topo = "(CD,FD)"
    elif event.detector1 == 2 and event.detector2 == 0:
        event_topo = "(CD,FT)"
    else:
        return False
    if event_topo != topology:
        return False
    theta_val = getattr(event, "theta_gamma_gamma", None)
    if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, theta_val,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "data", run_period, topology
        ):
        return False
    if not passes_3sigma_cuts(event, False, cuts):
        return False
    return True
#enddef

def get_tree_entries_with_cuts(file_path, cuts, topology, analysis_type, run_period):
    """
    Opens the ROOT file and counts the number of entries in the "PhysicsEvents" tree
    that pass the exclusivity cuts.
    """
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: 'PhysicsEvents' tree not found in file {file_path}")
        f.Close()
        return 0
    count = 0
    for event in tree:
        try:
            if passes_exclusivity_cuts(event, cuts, topology, analysis_type, run_period):
                count += 1
        except Exception as e:
            print(f"Exception in applying exclusivity cuts: {e}")
            continue
    f.Close()
    return count
#enddef

def get_total_tree_entries(file_path):
    """
    Opens the ROOT file and returns the total number of entries in the "PhysicsEvents" tree.
    """
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: 'PhysicsEvents' tree not found in file {file_path}")
        f.Close()
        return 0
    n = tree.GetEntries()
    f.Close()
    return n
#enddef

def plot_normalized_efficiencies(output_dir):
    """
    Calculates and plots overall normalized efficiencies (integrated over kinematics)
    as a function of beam current for each run period.
    Efficiency is computed as (n_reco/n_gen) normalized to the 0 nA (nobkg) value.
    A forced linear fit through (0,1) is performed.
    """
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/efficiency_study"
    run_periods = {
        "rga_fa18_inb": "DVCS_Fa18_inb",
        "rga_fa18_out": "DVCS_Fa18_out",
        "rga_sp19_inb": "DVCS_Sp19_inb"
    }
    title_map = {
        "rga_fa18_inb": "RGA Fa18 Inb",
        "rga_fa18_out": "RGA Fa18 Out",
        "rga_sp19_inb": "RGA Sp19 Inb"
    }
    selected_topology = "(FD,FD)"
    analysis_type = "dvcs"
    for run_prefix, period_code in run_periods.items():
        pattern_gen = os.path.join(base_dir, f"gen_{run_prefix}_dvcs_*.root")
        pattern_reco = os.path.join(base_dir, f"{run_prefix}_dvcs_*.root")
        gen_files = glob.glob(pattern_gen)
        reco_files = glob.glob(pattern_reco)
        gen_dict = {}
        reco_dict = {}
        def extract_current(filename):
            base = os.path.basename(filename)
            parts = base.split("_")
            last_part = parts[-1].replace(".root", "")
            if "nobkg" in last_part:
                return 0
            else:
                try:
                    return int(last_part.replace("nA", ""))
                except:
                    return None
        for fpath in gen_files:
            cur = extract_current(fpath)
            if cur is not None:
                gen_dict[cur] = fpath
        for fpath in reco_files:
            cur = extract_current(fpath)
            if cur is not None:
                reco_dict[cur] = fpath
        if 0 not in gen_dict or 0 not in reco_dict:
            print(f"Missing nobkg files for run period {run_prefix}. Skipping.")
            continue
        cuts = load_cuts(period_code, selected_topology)
        currents = []
        efficiencies = []
        efficiency_errors = []
        n_gen_0 = get_total_tree_entries(gen_dict[0])
        n_reco_0 = get_tree_entries_with_cuts(reco_dict[0], cuts, selected_topology, analysis_type, period_code)
        if n_gen_0 == 0:
            print(f"Warning: 0 generated entries for normalization in run period {run_prefix}")
            norm_eff = 1.0
        else:
            norm_eff = (n_reco_0 / n_gen_0)
        for cur in sorted(gen_dict.keys()):
            if cur == 0:
                currents.append(0)
                efficiencies.append(1.0)
                efficiency_errors.append(0.0)
                continue
            n_gen = get_total_tree_entries(gen_dict[cur])
            n_reco = get_tree_entries_with_cuts(reco_dict[cur], cuts, selected_topology, analysis_type, period_code)
            if n_gen == 0:
                eff = 0
                err = 0
            else:
                eff = (n_reco / n_gen) / norm_eff
                err = (np.sqrt(n_reco) / n_gen) / norm_eff if n_reco > 0 else 0
            currents.append(cur)
            efficiencies.append(eff)
            efficiency_errors.append(err)
        currents_arr = np.array(currents)
        efficiencies_arr = np.array(efficiencies)
        weights2 = np.array([1/(e**2) if e > 0 else 1 for e in efficiency_errors])
        numerator = np.sum(weights2 * currents_arr * (efficiencies_arr - 1))
        denominator = np.sum(weights2 * currents_arr**2)
        m = numerator / denominator if denominator != 0 else 0
        sigma_m = np.sqrt(1 / denominator) if denominator != 0 else 0
        model = 1 + m * currents_arr
        chi2 = np.sum(weights2 * (efficiencies_arr - model)**2)
        ndf = len(currents_arr) - 1
        chi2_ndf = chi2 / ndf if ndf > 0 else 0
        fit_x = np.linspace(min(currents_arr), max(currents_arr), 100)
        fit_y = 1 + m * fit_x
        plt.figure()
        plt.errorbar(currents_arr, efficiencies_arr, yerr=efficiency_errors, fmt='ko', label='Data')
        plt.plot(fit_x, fit_y, 'r--', label=f"m = {m:.4f} ± {sigma_m:.4f}\n" + r"$\chi^{2}/ndf$ = " + f"{chi2_ndf:.2f}")
        plt.xlabel("current (nA)")
        plt.ylabel("normalized efficiency")
        plt.xlim(-5, 60)
        plt.ylim(0.7, 1.05)
        plt.legend(loc='upper right')
        plt.title(f"Normalized Efficiency for {title_map.get(run_prefix, run_prefix)}")
        output_path = os.path.join(output_dir, f"{run_prefix}_integrated.pdf")
        plt.savefig(output_path, format='pdf')
        plt.close()
        print(f"Saved plot for {run_prefix} to {output_path}")
#enddef

#############################################
# New Functions for Binned Efficiency Analysis (Efficient Version with Debug)
#############################################

# Global φ binning parameters are defined above.

def count_generated_in_bin(file_path, bin_idx, unique_xB_bins, unique_Q2_bins, unique_t_bins):
    """
    Opens the generated ROOT file and counts events falling into the bin defined by bin_idx.
    Uses event.x for xB, event.Q2 for Q², |event.t1| for t, and event.phi2 for φ.
    """
    if DEBUG:
        print(f"count_generated_in_bin: Opening generated file {file_path}")
    count = 0
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: 'PhysicsEvents' tree not found in file {file_path}")
        f.Close()
        return 0
    processed = 0
    for event in tree:
        processed += 1
        try:
            xB_val = float(event.x)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t1))
            phi_val = float(event.phi2)
        except Exception as e:
            continue
        i_xB = find_bin(xB_val, unique_xB_bins)
        i_Q2 = find_bin(Q2_val, unique_Q2_bins)
        i_t  = find_bin(t_val, unique_t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        # Correct i_phi if it equals N_PHI_BINS.
        if i_phi >= N_PHI_BINS:
            if DEBUG:
                print(f"count_generated_in_bin: Correcting i_phi from {i_phi} to {N_PHI_BINS - 1} for phi_val = {phi_val}")
            i_phi = N_PHI_BINS - 1
        if (i_xB, i_Q2, i_t, i_phi) == bin_idx:
            count += 1
        if DEBUG and processed % 10000 == 0:
            print(f"count_generated_in_bin: Processed {processed} events; current count = {count}")
    f.Close()
    if DEBUG:
        print(f"count_generated_in_bin: Finished file {file_path} with count = {count}")
    return count
#enddef

def count_reco_in_bin(file_path, bin_idx, unique_xB_bins, unique_Q2_bins, unique_t_bins,
                        cuts_dict, analysis_type, run_period):
    """
    Opens the reconstructed ROOT file and counts events falling into the bin defined by bin_idx.
    Only events that pass the exclusivity cuts (using the efficiency-specific function) are counted.
    """
    if DEBUG:
        print(f"count_reco_in_bin: Opening reconstructed file {file_path}")
    count = 0
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: 'PhysicsEvents' tree not found in file {file_path}")
        f.Close()
        return 0
    processed = 0
    for event in tree:
        processed += 1
        try:
            if not passes_exclusivity_cuts_eff(event, cuts_dict, analysis_type, run_period):
                continue
            xB_val = float(event.x)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t1))
            phi_val = float(event.phi2)
        except Exception as e:
            continue
        i_xB = find_bin(xB_val, unique_xB_bins)
        i_Q2 = find_bin(Q2_val, unique_Q2_bins)
        i_t  = find_bin(t_val, unique_t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        # Correct i_phi if it equals or exceeds N_PHI_BINS.
        if i_phi >= N_PHI_BINS:
            if DEBUG:
                print(f"count_reco_in_bin: Correcting i_phi from {i_phi} to {N_PHI_BINS - 1} for phi_val = {phi_val}")
            i_phi = N_PHI_BINS - 1
        if (i_xB, i_Q2, i_t, i_phi) == bin_idx:
            count += 1
        if DEBUG and processed % 10000 == 0:
            print(f"count_reco_in_bin: Processed {processed} events; current count = {count}")
    f.Close()
    if DEBUG:
        print(f"count_reco_in_bin: Finished file {file_path} with count = {count}")
    return count
#enddef

def calculate_binned_efficiencies_all(binning_scheme, output_dir):
    """
    Loops over all run periods and calculates binned efficiencies.
    For each run period, for each bin defined by unique xB, Q², t and φ,
    the efficiency is computed as (n_reco / total_generated) where total_generated
    is obtained from the generated ROOT file via GetEntries.
    A forced linear fit through (0,1) (model: y = m*x + 1) is then performed for each bin.
    The slope m and its uncertainty are saved in a JSON file (one per run period).
    """
    if DEBUG:
        print("calculate_binned_efficiencies_all: Starting calculation for all run periods.")
    RUN_PERIODS = {
        "rga_fa18_inb": "DVCS_Fa18_inb",
        "rga_fa18_out": "DVCS_Fa18_out",
        "rga_sp19_inb": "DVCS_Sp19_inb"
    }
    analysis_type = "dvcs"
    unique_xB_bins = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    unique_Q2_bins = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    unique_t_bins  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/efficiency_study"
    
    def extract_current(filename):
        base = os.path.basename(filename)
        parts = base.split("_")
        last_part = parts[-1].replace(".root", "")
        if "nobkg" in last_part:
            return 0
        else:
            try:
                return int(last_part.replace("nA", ""))
            except:
                return None
    #enddef
    
    for run_prefix, period_code in RUN_PERIODS.items():
        if DEBUG:
            print(f"\nProcessing run period {run_prefix} ({period_code})")
        # Load cuts for each topology.
        cuts_dict = {
            "(FD,FD)": load_cuts(period_code, "(FD,FD)"),
            "(CD,FD)": load_cuts(period_code, "(CD,FD)"),
            "(CD,FT)": load_cuts(period_code, "(CD,FT)")
        }
        if DEBUG:
            print(f"Loaded cuts for topologies: {list(cuts_dict.keys())}")
        pattern_gen = os.path.join(base_dir, f"gen_{run_prefix}_dvcs_*.root")
        pattern_reco = os.path.join(base_dir, f"{run_prefix}_dvcs_*.root")
        gen_files = glob.glob(pattern_gen)
        reco_files = glob.glob(pattern_reco)
        if DEBUG:
            print(f"Found {len(gen_files)} generated files and {len(reco_files)} reconstructed files for {run_prefix}")
        gen_dict = {extract_current(fp): fp for fp in gen_files if extract_current(fp) is not None}
        reco_dict = {extract_current(fp): fp for fp in reco_files if extract_current(fp) is not None}
        if 0 not in gen_dict or 0 not in reco_dict:
            print(f"Missing nobkg files for run period {run_prefix}. Skipping.")
            continue

        # Get total generated counts (using GetEntries, no event loop).
        total_gen = {}
        for cur, file_path in gen_dict.items():
            f_gen = ROOT.TFile.Open(file_path)
            if not f_gen or f_gen.IsZombie():
                print(f"Error opening generated file {file_path}")
                total_gen[cur] = 0
            else:
                tree_gen = f_gen.Get("PhysicsEvents")
                if tree_gen:
                    total_gen[cur] = tree_gen.GetEntries()
                else:
                    print(f"Generated tree not found in {file_path}")
                    total_gen[cur] = 0
            f_gen.Close()
        if DEBUG:
            print(f"Total generated counts: {total_gen}")

        # Loop through each reconstructed file once.
        rec_counts = {}  # rec_counts[current][bin_idx] = count
        for cur, file_path in reco_dict.items():
            if DEBUG:
                print(f"Processing reconstructed file for current {cur} nA: {file_path}")
            rec_counts[cur] = {}
            # Initialize counts for every possible bin.
            for i_xB in range(len(unique_xB_bins)):
                for i_Q2 in range(len(unique_Q2_bins)):
                    for i_t in range(len(unique_t_bins)):
                        for i_phi in range(N_PHI_BINS):
                            rec_counts[cur][(i_xB, i_Q2, i_t, i_phi)] = 0
            f_rec = ROOT.TFile.Open(file_path)
            if not f_rec or f_rec.IsZombie():
                print(f"Error opening reconstructed file {file_path}")
                continue
            tree_rec = f_rec.Get("PhysicsEvents")
            if not tree_rec:
                print(f"Reconstructed tree not found in {file_path}")
                f_rec.Close()
                continue
            processed = 0
            for event in tree_rec:
                processed += 1
                if processed % 10000 == 0 and DEBUG:
                    print(f"Processed {processed} events in current {cur}")
                if not passes_exclusivity_cuts_eff(event, cuts_dict, analysis_type, period_code):
                    continue
                try:
                    xB_val = float(event.x)
                    Q2_val = float(event.Q2)
                    t_val  = abs(float(event.t1))
                    phi_val = float(event.phi2)
                except Exception as e:
                    continue
                i_xB = find_bin(xB_val, unique_xB_bins)
                i_Q2 = find_bin(Q2_val, unique_Q2_bins)
                i_t  = find_bin(t_val, unique_t_bins)
                i_phi = np.digitize(phi_val, phi_edges) - 1
                if i_phi >= N_PHI_BINS:
                    if DEBUG:
                        print(f"Correcting i_phi from {i_phi} to {N_PHI_BINS-1} for phi_val = {phi_val}")
                    i_phi = N_PHI_BINS - 1
                if i_xB is None or i_Q2 is None or i_t is None or i_phi is None:
                    continue
                rec_counts[cur][(i_xB, i_Q2, i_t, i_phi)] += 1
            f_rec.Close()
            if DEBUG:
                print(f"Finished processing current {cur}. Total events processed: {processed}")
        # Build efficiency data for each bin.
        efficiency_data = {}
        for i_xB in range(len(unique_xB_bins)):
            for i_Q2 in range(len(unique_Q2_bins)):
                for i_t in range(len(unique_t_bins)):
                    for i_phi in range(N_PHI_BINS):
                        efficiency_data[(i_xB, i_Q2, i_t, i_phi)] = {"currents": [], "eff": [], "err": []}
        for cur in rec_counts.keys():
            for bin_idx in efficiency_data.keys():
                rec = rec_counts[cur].get(bin_idx, 0)
                tot = total_gen.get(cur, 0)
                if tot > 0:
                    eff = rec / tot
                    err = (np.sqrt(rec) / tot) if rec > 0 else 0
                else:
                    eff = 0
                    err = 0
                efficiency_data[bin_idx]["currents"].append(cur)
                efficiency_data[bin_idx]["eff"].append(eff)
                efficiency_data[bin_idx]["err"].append(err)
        # Normalize each bin's efficiencies by the 0 nA value.
        for bin_idx, data in efficiency_data.items():
            currents_arr = np.array(data["currents"])
            eff_arr = np.array(data["eff"])
            err_arr = np.array(data["err"])
            norm = 1.0
            zero_idx = np.where(currents_arr == 0)[0]
            if len(zero_idx) > 0:
                norm = eff_arr[zero_idx[0]]
            if norm == 0:
                norm = 1.0
            eff_arr = eff_arr / norm
            err_arr = err_arr / norm
            data["eff"] = eff_arr.tolist()
            data["err"] = err_arr.tolist()
        if DEBUG:
            print(f"Normalized efficiency data for {run_prefix}.")
        # For each bin, perform a forced linear fit through (0,1): model y = m*x + 1.
        fit_results = {}
        for bin_idx, data in efficiency_data.items():
            currents_arr = np.array(data["currents"])
            eff_arr = np.array(data["eff"])
            err_arr = np.array(data["err"])
            weights = np.array([1/(e**2) if e > 0 else 1 for e in err_arr])
            numerator = np.sum(weights * currents_arr * (eff_arr - 1))
            denominator = np.sum(weights * currents_arr**2)
            m = numerator / denominator if denominator != 0 else 0
            sigma_m = np.sqrt(1/denominator) if denominator != 0 else 0
            fit_results[bin_idx] = {"m": m, "m_err": sigma_m}
            if DEBUG:
                print(f"Bin {bin_idx}: Fit slope m = {m:.4f}, sigma_m = {sigma_m:.4f}")
        run_output_path = os.path.join(output_dir, f"efficiencies_{run_prefix}.json")
        fit_results_json = {str(k): v for k, v in fit_results.items()}
        with open(run_output_path, "w") as f:
            json.dump(fit_results_json, f, indent=2)
        print(f"Saved binned efficiency fit results for {run_prefix} to {run_output_path}")
#enddef

def plot_binned_efficiencies_all(binning_csv, output_dir, efficiencies_json_dir=None):
    """
    Reads the JSON files of binned efficiency fit results (one per run period) and the binning scheme,
    then creates a canvas for each overall xB bin. In each canvas, subplots are arranged by Q² (rows)
    and t (columns), and in each subplot the slope m (with error bars) versus φ (in degrees) is plotted.
    """
    if DEBUG:
        print("plot_binned_efficiencies_all: Loading binning scheme.")
    binning_scheme = load_binning_scheme(binning_csv)
    unique_xB_bins = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    unique_Q2_bins = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    unique_t_bins  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    RUN_PERIODS = {
        "rga_fa18_inb": "DVCS_Fa18_inb",
        "rga_fa18_out": "DVCS_Fa18_out",
        "rga_sp19_inb": "DVCS_Sp19_inb"
    }
    if efficiencies_json_dir is None:
        efficiencies_json_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    for run_prefix, period_code in RUN_PERIODS.items():
        json_path = os.path.join(efficiencies_json_dir, f"efficiencies_{run_prefix}.json")
        if not os.path.exists(json_path):
            print(f"plot_binned_efficiencies_all: JSON file for {run_prefix} not found, skipping plotting for this period.")
            continue
        if DEBUG:
            print(f"plot_binned_efficiencies_all: Loading JSON for {run_prefix} from {json_path}")
        with open(json_path, "r") as f:
            data_json = json.load(f)
        efficiencies = {}
        for key_str, val in data_json.items():
            key_tuple = tuple(int(x.strip()) for x in key_str.strip("()").split(","))
            efficiencies[key_tuple] = val
        for i_xB in range(len(unique_xB_bins)):
            nrows = len(unique_Q2_bins)
            ncols = len(unique_t_bins)
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(3 * ncols + 2, 3 * nrows + 2),
                                     squeeze=False)
            for i_Q2 in range(len(unique_Q2_bins)):
                for i_t in range(len(unique_t_bins)):
                    ax = axes[i_Q2, i_t]
                    m_vals = []
                    m_errs = []
                    phi_vals = []
                    for i_phi in range(N_PHI_BINS):
                        key = (i_xB, i_Q2, i_t, i_phi)
                        if key in efficiencies:
                            m_vals.append(efficiencies[key]["m"])
                            m_errs.append(efficiencies[key]["m_err"])
                        else:
                            m_vals.append(0)
                            m_errs.append(0)
                        phi_vals.append(phi_centers_deg[i_phi])
                    ax.errorbar(phi_vals, m_vals, yerr=m_errs, fmt='o', color='black', capsize=3)
                    ax.set_xlim(0, 360)
                    ax.set_xticks([0, 90, 180, 270, 360])
                    ax.grid(True, linestyle='--', alpha=0.5)
                    ax.set_title(f"Q² bin {i_Q2}, t bin {i_t}", fontsize=9)
                    if i_Q2 == nrows - 1:
                        ax.set_xlabel(r"$\phi\ (\deg)$", fontsize=9)
                    if i_t == 0:
                        ax.set_ylabel("Slope m", fontsize=9)
            fig.suptitle(f"Efficiency Slope for {period_code}, xB bin {i_xB}", fontsize=12)
            fig.tight_layout(rect=[0, 0, 1, 0.95])
            outpath = os.path.join(output_dir, f"efficiency_slope_{run_prefix}_xB_{i_xB}.png")
            plt.savefig(outpath, dpi=150)
            plt.close(fig)
            if DEBUG:
                print(f"plot_binned_efficiencies_all: Saved efficiency slope plot for {run_prefix}, xB bin {i_xB} to {outpath}")
#enddef