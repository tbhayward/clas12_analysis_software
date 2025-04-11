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
# Existing Functions for Integrated Efficiency
#############################################

def load_cuts(period, topology):
    """
    Loads the final cuts dictionary from combined_cuts.json for a given (period, topology).
    Expect keys like "DVCS_Fa18_inb_FD_FD", etc.
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
    Applies exclusivity cuts to an event (for integrated analysis) given a fixed topology.
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
    Opens the ROOT file and counts the number of entries in the tree that pass the exclusivity cuts.
    Assumes the tree is named 'PhysicsEvents'.
    """
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: tree not found in file {file_path}")
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
    Opens the ROOT file and returns the total number of entries in the tree.
    Assumes the tree is named 'PhysicsEvents'.
    """
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: tree not found in file {file_path}")
        f.Close()
        return 0
    n = tree.GetEntries()
    f.Close()
    return n
#enddef

def plot_mc_normalized_efficiencies(output_dir):
    """
    Calculates and plots overall normalized efficiencies as a function of beam current
    (integrated over kinematics) for each run period. Efficiency is defined as:
        (n_reco/n_gen) normalized to the 0 nA (nobkg) value.
    A forced linear fit through (0,1) is performed.
    """
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/mc_efficiency_study"
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

import math
import numpy as np
import os
import glob
import json
import ROOT
import matplotlib.pyplot as plt

# Global debug flag
DEBUG = False

# Global φ binning parameters.
N_PHI_BINS = 12
phi_edges = np.linspace(0, 2 * math.pi, N_PHI_BINS + 1)
phi_centers = (phi_edges[:-1] + phi_edges[1:]) / 2.0
phi_centers_deg = np.degrees(phi_centers)

def find_bin(value, bin_boundaries):
    global DEBUG
    """
    Given a value and a list of (min, max) boundaries, return the bin index (0-indexed)
    if the value falls in a bin; otherwise, return None.
    """
    for i, (low, high) in enumerate(bin_boundaries):
        if low <= value < high:
            if DEBUG:
                print(f"find_bin: value {value} is in bin {i} with boundaries ({low}, {high})")
            return i
    if DEBUG:
        print(f"find_bin: value {value} did not fall into any bin.")
    return None
#enddef

def passes_exclusivity_cuts_eff(event, cuts_dict, analysis_type, run_period):
    global DEBUG
    """
    Applies exclusivity cuts for binned efficiency analysis.
    Determines the event's topology from detector1 and detector2, then applies kinematic and 3σ cuts
    using the appropriate cuts dictionary from cuts_dict.
    """
    if event.detector1 == 1 and event.detector2 == 1:
        event_topo = "(FD,FD)"
    elif event.detector1 == 2 and event.detector2 == 1:
        event_topo = "(CD,FD)"
    elif event.detector1 == 2 and event.detector2 == 0:
        event_topo = "(CD,FT)"
    else:
        if DEBUG:
            print("passes_exclusivity_cuts_eff: Event failed detector topology check.")
        return False
    if event_topo not in cuts_dict:
        if DEBUG:
            print(f"passes_exclusivity_cuts_eff: No cuts found for topology {event_topo}.")
        return False
    cuts = cuts_dict[event_topo]
    theta_val = getattr(event, "theta_gamma_gamma", None)
    if theta_val is None:
        theta_val = getattr(event, "theta_pi0_pi0", None)
    if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, theta_val,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "data", run_period, event_topo
        ):
        if DEBUG:
            print("passes_exclusivity_cuts_eff: Event failed kinematic cuts.")
        return False
    if not passes_3sigma_cuts(event, False, cuts):
        if DEBUG:
            print("passes_exclusivity_cuts_eff: Event failed 3σ cuts.")
        return False
    return True
#enddef

def plot_data_normalized_efficiencies(output_dir):
    """
    Calculates and plots the normalized efficiency as a function of beam current
    (integrated over kinematics) for each run period for data.
    
    For each run, the raw efficiency is defined as:
         E_raw = counts / accumulated Faraday cup charge.
         
    The normalized efficiency is then given by:
         norm_eff = (E_raw) / (E_baseline)
    where E_baseline is the raw efficiency of the lowest-current run available 
    (assumed to be the baseline). Note that for some triggers (e.g. for FA18 out with
    v22_1_no_FT.cnf) a 5 nA run may not be available—in that case the minimum available
    current is used.
    
    A linear fit is performed on the normalized efficiency vs. beam current using a
    free intercept and slope.
    
    Parameters:
      output_dir (str): Directory where the output plots (in PDF format) will be saved.
    """
    import os
    import glob
    import numpy as np
    import matplotlib.pyplot as plt

    # Base directory for the data efficiency .root files.
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data_efficiency_study"
    
    # Mapping of run period identifiers to period codes for loading cuts and for plotting.
    # Note that these keys exactly match the start of the filename for each trigger
    run_periods = {
        "rga_fa18_inb_v22_0.cnf": "DVCS_Fa18_inb",                  # FA18 inb with v22_0.cnf trigger
        "rga_fa18_out_v22_0_no_dc_roads.cnf": "DVCS_Fa18_out",         # FA18 out with v22_0_no_dc_roads.cnf trigger
        "rga_fa18_out_v22_1_no_dc_roads.cnf": "DVCS_Fa18_out",         # FA18 out with v22_1_no_dc_roads.cnf trigger
        "rga_fa18_out_v22_1_no_FT.cnf": "DVCS_Fa18_out",               # FA18 out with v22_1_no_FT.cnf trigger
        "rga_sp19_inb": "DVCS_Sp19_inb"                                # SP19 inb remains the same trigger
    }
    
    # Title map for more friendly plot titles.
    title_map = {
        "rga_fa18_inb_v22_0.cnf": "RGA Fa18 Inb, v22_0.cnf",
        "rga_fa18_out_v22_0_no_dc_roads.cnf": "RGA Fa18 Out, v22_0_no_dc_roads.cnf",
        "rga_fa18_out_v22_1_no_dc_roads.cnf": "RGA Fa18 Out, v22_1_no_dc_roads.cnf",
        "rga_fa18_out_v22_1_no_FT.cnf": "RGA Fa18 Out, v22_1_no_FT.cnf",
        "rga_sp19_inb": "RGA Sp19 Inb"
    }
    
    # Mapping of accumulated Faraday cup charge (in nC) for each run period and beam current.
    # Each inner dictionary maps a beam current (nA) to its corresponding charge.
    charge_map = {
        "rga_fa18_inb_v22_0.cnf": {
            5: 93082,
            45: 1.55491e6,
            50: 2.32751e6,
            55: 2.00433e6
        },
        "rga_fa18_out_v22_0_no_dc_roads.cnf": {
            5: 56198.4,
            20: 116623.69,
            40: 2.23767e6
        },
        "rga_fa18_out_v22_1_no_dc_roads.cnf": {
            5: 79710.5,
            50: 1.82106e6
        },
        "rga_fa18_out_v22_1_no_FT.cnf": {
            50: 14739.981,
            70: 83940.98
        },
        "rga_sp19_inb": {
            5: 30860.979,
            10: 126675.44,
            50: 3.57217e6
        }
    }
    
    # Analysis parameters (common for all runs).
    selected_topology = "(FD,FD)"
    analysis_type = "dvcs"
    
    # Loop over each run period (trigger configuration).
    for run_prefix, period_code in run_periods.items():
        # Build the file pattern to collect all .root files for the current run period.
        # This will match files such as 'rga_fa18_inb_v22_0.cnf_dvcs_45nA.root' etc.
        pattern = os.path.join(base_dir, f"{run_prefix}_dvcs_*.root")
        file_list = glob.glob(pattern)
        
        # Dictionary to store files keyed by the extracted beam current (in nA)
        data_dict = {}
        
        # Function to extract the beam current from the filename.
        # Expects filenames like "rga_fa18_inb_v22_0.cnf_dvcs_40nA.root"
        def extract_current(filename):
            """
            Extracts the beam current (in nA) from the given filename.
            """
            base = os.path.basename(filename)
            parts = base.split("_")
            # The beam current information is assumed to be in the last part, e.g. "45nA.root"
            last_part = parts[-1].replace(".root", "")
            try:
                return int(last_part.replace("nA", ""))
            except:
                return None
        #enddef
        
        # Loop over the files found and populate data_dict with files that have valid current values.
        for fpath in file_list:
            current = extract_current(fpath)
            if current is not None:
                # Check if a corresponding accumulated charge exists in our mapping for this run period.
                if current in charge_map.get(run_prefix, {}):
                    data_dict[current] = fpath
                else:
                    print(f"Warning: No charge mapping for run period {run_prefix} at {current} nA. Skipping file {fpath}.")
            #endif
        #endfor
        
        # Check that we have at least one valid data file.
        if not data_dict:
            print(f"No data files found for run period {run_prefix}. Skipping.")
            continue
        #endif
        
        # Determine the baseline run using the lowest available beam current.
        # (For most FA18 runs this is expected to be 5 nA. For triggers missing a 5 nA run, the lowest available is used.)
        baseline_current = min(data_dict.keys())
        if baseline_current not in charge_map.get(run_prefix, {}):
            print(f"Baseline current {baseline_current} nA not found in charge mapping for run period {run_prefix}. Skipping.")
            continue
        #endif
        
        # Load the analysis cuts for the current period and selected topology.
        cuts = load_cuts(period_code, selected_topology)
        
        # Calculate the baseline efficiency from the baseline run.
        baseline_file = data_dict[baseline_current]
        baseline_counts = get_tree_entries_with_cuts(baseline_file, cuts, selected_topology, analysis_type, period_code)
        baseline_charge = charge_map[run_prefix][baseline_current]
        if baseline_charge == 0:
            print(f"Warning: Baseline charge is zero for run period {run_prefix} at {baseline_current} nA. Skipping.")
            continue
        #endif
        baseline_efficiency = baseline_counts / baseline_charge
        if baseline_efficiency == 0:
            print(f"Warning: Baseline efficiency is zero for run period {run_prefix} at {baseline_current} nA. Skipping.")
            continue
        #endif
        
        # Initialize lists to hold beam currents, normalized efficiencies, and their errors.
        currents = []
        norm_efficiencies = []
        norm_eff_errors = []
        
        # Loop over each beam current (sorted in ascending order) in the current run period.
        for current in sorted(data_dict.keys()):
            filename = data_dict[current]
            # Get the counts using the given cuts and analysis type.
            counts = get_tree_entries_with_cuts(filename, cuts, selected_topology, analysis_type, period_code)
            charge = charge_map[run_prefix].get(current, None)
            if charge is None or charge == 0:
                print(f"Warning: Invalid charge for run period {run_prefix} at {current} nA. Skipping this run.")
                continue
            #endif
            # Calculate the raw efficiency as counts per unit charge.
            raw_efficiency = counts / charge
            # Compute the normalized efficiency relative to the baseline efficiency.
            normalized_efficiency = raw_efficiency / baseline_efficiency
            # Propagate the Poisson uncertainty (sqrt(counts)) on counts.
            error = (np.sqrt(counts) / charge) / baseline_efficiency if counts > 0 else 0
            
            currents.append(current)
            norm_efficiencies.append(normalized_efficiency)
            norm_eff_errors.append(error)
        #endfor
        
        # Convert the lists to numpy arrays for ease of calculation.
        currents_arr = np.array(currents)
        norm_eff_arr = np.array(norm_efficiencies)
        
        # Define weights for the linear fit based on uncertainties (avoid division by zero).
        weights2 = np.array([1/(err**2) if err > 0 else 1 for err in norm_eff_errors])
        
        # --- Linear Fit (free intercept and slope) ---
        # Fit the model: norm_eff = b + m * current
        w = weights2
        x = currents_arr
        y = norm_eff_arr
        S = np.sum(w)
        Sx = np.sum(w * x)
        Sy = np.sum(w * y)
        Sxx = np.sum(w * x**2)
        Sxy = np.sum(w * x * y)
        Delta = S * Sxx - Sx**2

        if Delta != 0:
            m = (S * Sxy - Sx * Sy) / Delta
            b = (Sy * Sxx - Sx * Sxy) / Delta
            sigma_m = np.sqrt(S / Delta)
            sigma_b = np.sqrt(Sxx / Delta)
        else:
            m, b = 0, 0
            sigma_m, sigma_b = 0, 0
        
        # Calculate chi-squared for the fit.
        model = b + m * x
        chi2 = np.sum(w * (y - model)**2)
        ndf = len(x) - 2  # two free parameters: intercept and slope
        chi2_ndf = chi2 / ndf if ndf > 0 else 0
        
        # Prepare data to plot the fit line across the range of beam currents.
        fit_currents = np.linspace(min(x), max(x), 100)
        fit_norm_eff = b + m * fit_currents
        
        # --- Plotting ---
        plt.figure()
        plt.errorbar(currents_arr, norm_eff_arr, yerr=norm_eff_errors, fmt='ko', label='Data')
        plt.plot(fit_currents, fit_norm_eff, 'r--', 
                 label=f"m = {m:.4f} ± {sigma_m:.4f}\n" + r"$\chi^{2}/ndf$ = " + f"{chi2_ndf:.2f}")
        plt.xlabel("Current (nA)")
        plt.ylabel("Normalized Efficiency")
        plt.xlim(-5, 80)
        plt.ylim(0.0, 1.5)  # Adjust y-axis limits as needed.
        plt.legend(loc='upper right')
        # Use a friendly title if available; otherwise default to the run_prefix.
        plt.title(f"Normalized Efficiency for {title_map.get(run_prefix, run_prefix)}")
        
        # Save the plot as a PDF in the given output directory.
        output_path = os.path.join(output_dir, f"{run_prefix}_integrated.pdf")
        plt.savefig(output_path, format='pdf')
        plt.close()
        print(f"Saved plot for {run_prefix} to {output_path}")
    #endfor
#enddef


def plot_DAF_data_normalized_efficiencies(output_dir):
    """
    Calculates and plots the normalized efficiency as a function of beam current
    (integrated over kinematics) for each run period for data.
    
    For each run, the raw efficiency is defined as:
        E_raw = counts / accumulated Faraday cup charge.
        
    The normalized efficiency is then given by:
        norm_eff = (E_raw) / (E_baseline)
    where E_baseline is the raw efficiency of the 5 nA (lowest current) run.
    
    A linear fit is performed on the normalized efficiency vs. beam current using a
    forced fit through the baseline point (i.e. norm_eff = 1 at 5 nA).
    
    Parameters:
      output_dir (str): Directory where the output plots (in PDF format) will be saved.
    """
    # Base directory for the data efficiency .root files.
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/data_efficiency_study"
    
    # Mapping of run period identifiers to period codes and titles for plotting.
    run_periods = {
        "rga_fa18_inb_DAF": "DVCS_Fa18_inb",
        "rga_fa18_out_DAF": "DVCS_Fa18_out"
    }
    title_map = {
        "rga_fa18_inb_DAF": "RGA Fa18 Inb",
        "rga_fa18_out_DAF": "RGA Fa18 Out"
    }
    
    charge_map = {
        "rga_fa18_inb_DAF": {
            5: 159661.55,
            45: 1.25969e6,
        },
        "rga_fa18_out_DAF": {
            5: 166672.03,
            20: 177239.66,
            40: 3.6206e6,
        }
    }
    
    # Analysis parameters.
    selected_topology = "(FD,FD)"
    analysis_type = "dvcs"
    
    # Loop over each run period (e.g., Fall 2018 Inb, Fall 2018 Out, Spring 2019 Inb).
    for run_prefix, period_code in run_periods.items():
        # Define the file pattern to collect all .root files for the current run period.
        pattern = os.path.join(base_dir, f"{run_prefix}_dvcs_*.root")
        file_list = glob.glob(pattern)
        
        # Dictionary to store files keyed by the extracted beam current (in nA).
        data_dict = {}
        
        # Function to extract the beam current from the filename.
        # Expects filenames like "rga_fa18_inb_dvcs_40nA.root".
        def extract_current(filename):
            """
            Extracts the beam current (in nA) from the given filename.
            """
            base = os.path.basename(filename)
            parts = base.split("_")
            last_part = parts[-1].replace(".root", "")
            try:
                # Remove the trailing "nA" and convert to integer.
                return int(last_part.replace("nA", ""))
            except:
                return None
        #enddef
        
        # Loop over the files found and populate data_dict with files that have valid current values.
        for fpath in file_list:
            current = extract_current(fpath)
            if current is not None:
                # Check if a corresponding accumulated charge exists in our mapping for this run period.
                if current in charge_map.get(run_prefix, {}):
                    data_dict[current] = fpath
                else:
                    print(f"Warning: No charge mapping for run period {run_prefix} at {current} nA. Skipping file {fpath}.")
            #endif
        #endfor
        
        # Check that we have at least one valid data file.
        if not data_dict:
            print(f"No data files found for run period {run_prefix}. Skipping.")
            continue
        #endif
        
        # Determine the baseline run using the lowest available beam current.
        # Here we assume the baseline is at 5 nA.
        baseline_current = min(data_dict.keys())
        if baseline_current not in charge_map.get(run_prefix, {}):
            print(f"Baseline current {baseline_current} nA not found in charge mapping for run period {run_prefix}. Skipping.")
            continue
        #endif
        
        # Load the analysis cuts for the current period and selected topology.
        cuts = load_cuts(period_code, selected_topology)
        
        # Calculate the baseline efficiency from the baseline run.
        baseline_file = data_dict[baseline_current]
        baseline_counts = get_tree_entries_with_cuts(baseline_file, cuts, selected_topology, analysis_type, period_code)
        baseline_charge = charge_map[run_prefix][baseline_current]
        if baseline_charge == 0:
            print(f"Warning: Baseline charge is zero for run period {run_prefix} at {baseline_current} nA. Skipping.")
            continue
        #endif
        baseline_efficiency = baseline_counts / baseline_charge
        if baseline_efficiency == 0:
            print(f"Warning: Baseline efficiency is zero for run period {run_prefix} at {baseline_current} nA. Skipping.")
            continue
        #endif
        
        # Lists to hold beam currents, normalized efficiencies, and their errors for plotting and fitting.
        currents = []
        norm_efficiencies = []
        norm_eff_errors = []
        
        # Loop over each beam current (sorted from lowest to highest) in the current run period.
        for current in sorted(data_dict.keys()):
            filename = data_dict[current]
            counts = get_tree_entries_with_cuts(filename, cuts, selected_topology, analysis_type, period_code)
            charge = charge_map[run_prefix].get(current, None)
            if charge is None or charge == 0:
                print(f"Warning: Invalid charge for run period {run_prefix} at {current} nA. Skipping this run.")
                continue
            #endif
            # Calculate the raw efficiency as counts per unit charge.
            raw_efficiency = counts / charge
            # Compute the normalized efficiency relative to the baseline.
            normalized_efficiency = raw_efficiency / baseline_efficiency
            # Propagate the Poisson uncertainty (sqrt(counts)) on counts.
            error = (np.sqrt(counts) / charge) / baseline_efficiency if counts > 0 else 0
            
            currents.append(current)
            norm_efficiencies.append(normalized_efficiency)
            norm_eff_errors.append(error)
        #endfor
        
        # Convert lists to numpy arrays.
        currents_arr = np.array(currents)
        norm_eff_arr = np.array(norm_efficiencies)
        
        # Define weights for the linear fit (avoid division by zero).
        weights2 = np.array([1/(err**2) if err > 0 else 1 for err in norm_eff_errors])
        
                # --- Linear Fit (free intercept and slope) ---
        # Fit the model: norm_eff = b + m * current
        # Define weights for the fit using norm_eff_errors:
        w = weights2
        x = currents_arr
        y = norm_eff_arr
        S = np.sum(w)
        Sx = np.sum(w * x)
        Sy = np.sum(w * y)
        Sxx = np.sum(w * x**2)
        Sxy = np.sum(w * x * y)
        Delta = S * Sxx - Sx**2

        if Delta != 0:
            m = (S * Sxy - Sx * Sy) / Delta
            b = (Sy * Sxx - Sx * Sxy) / Delta
            sigma_m = np.sqrt(S / Delta)
            sigma_b = np.sqrt(Sxx / Delta)
        else:
            m, b = 0, 0
            sigma_m, sigma_b = 0, 0
        
        # Calculate chi-squared for the fit.
        model = b + m * x
        chi2 = np.sum(w * (y - model)**2)
        ndf = len(x) - 2  # two free parameters: slope and intercept
        chi2_ndf = chi2 / ndf if ndf > 0 else 0
        
        # Prepare data for plotting the fit line over the beam current range.
        fit_currents = np.linspace(min(x), max(x), 100)
        fit_norm_eff = b + m * fit_currents
        
        # --- Plotting ---
        plt.figure()
        plt.errorbar(currents_arr, norm_eff_arr, yerr=norm_eff_errors, fmt='ko', label='Data')
        plt.plot(fit_currents, fit_norm_eff, 'r--', 
                 label=f"m = {m:.4f} ± {sigma_m:.4f}\n" + r"$\chi^{2}/ndf$ = " + f"{chi2_ndf:.2f}")
        plt.xlabel("Current (nA)")
        plt.ylabel("Normalized Efficiency")
        plt.xlim(-5, 60)
        # plt.ylim(0.0, 1.5)  # Set y-axis limits as requested.
        plt.legend(loc='upper right')
        plt.title(f"Normalized Efficiency for {title_map.get(run_prefix, run_prefix)}")
        
        # Save the plot as a PDF in the output directory.
        output_path = os.path.join(output_dir, f"{run_prefix}_DAF_integrated.pdf")
        plt.savefig(output_path, format='pdf')
        plt.close()
        print(f"Saved plot for {run_prefix} to {output_path}")
    #endfor
#enddef

#############################################
# New Functions for Binned Efficiency Analysis (Efficient Version with Debug)
#############################################

# Global φ binning parameters are defined above.

def count_generated_in_bin(file_path, bin_idx, unique_xB_bins, unique_Q2_bins, unique_t_bins):
    global DEBUG
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
    global DEBUG
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
    global DEBUG
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

        # Get total generated counts using GetEntries.
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
    global DEBUG
    """
    Reads the JSON files of binned efficiency fit results (one per run period) and the binning scheme,
    then creates a canvas for each overall xB bin.
    Unlike before, each canvas will have only as many subplots (for Q² and t) as exist in that xB slice.
    In each subplot the slope m (with error bars) versus φ (in degrees) is plotted.
    The y–axis in every subplot is fixed from -0.1 to 0.1.
    """
    if DEBUG:
        print("plot_binned_efficiencies_all: Loading binning scheme.")
    binning_scheme = load_binning_scheme(binning_csv)
    # Determine overall unique xB, Q², and t boundaries.
    overall_unique_xB = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    overall_unique_Q2 = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    overall_unique_t = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    RUN_PERIODS = {
        "rga_fa18_inb": "DVCS_Fa18_inb",
        "rga_fa18_out": "DVCS_Fa18_out",
        "rga_sp19_inb": "DVCS_Sp19_inb"
    }
    if efficiencies_json_dir is None:
        efficiencies_json_dir = "output"
    os.makedirs(output_dir, exist_ok=True)
    
    # Loop over each run period.
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
        
        # For each overall xB bin, determine which Q² and t bins are present.
        for i_xB, xB_bound in enumerate(overall_unique_xB):
            # Collect Q² and t indices that appear for this xB.
            q2_indices = set()
            t_indices = set()
            for key in efficiencies.keys():
                if key[0] == i_xB:
                    q2_indices.add(key[1])
                    t_indices.add(key[2])
            q2_indices = sorted(list(q2_indices))
            t_indices = sorted(list(t_indices))
            if DEBUG:
                print(f"For run {run_prefix}, xB bin {i_xB}: Found Q² indices {q2_indices} and t indices {t_indices}")
            # If no Q² or t bins, skip.
            if len(q2_indices) == 0 or len(t_indices) == 0:
                continue
            nrows = len(q2_indices)
            ncols = len(t_indices)
            fig, axes = plt.subplots(nrows, ncols,
                                     figsize=(3 * ncols + 2, 3 * nrows + 2),
                                     squeeze=False)
            # Loop over the available Q² and t indices.
            for r, i_Q2 in enumerate(q2_indices):
                for c, i_t in enumerate(t_indices):
                    ax = axes[r, c]
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
                    ax.set_ylim(-0.1, 0.1)
                    ax.grid(True, linestyle='--', alpha=0.5)
                    ax.set_title(f"Q² idx {i_Q2}, t idx {i_t}", fontsize=9)
                    if r == nrows - 1:
                        ax.set_xlabel(r"$\phi\ (\deg)$", fontsize=9)
                    if c == 0:
                        ax.set_ylabel("Slope m", fontsize=9)
            fig.suptitle(f"Efficiency Slope for {period_code}, xB bin {i_xB}", fontsize=12)
            fig.tight_layout(rect=[0, 0, 1, 0.95])
            outpath = os.path.join(output_dir, f"efficiency_slope_{run_prefix}_xB_{i_xB}.png")
            plt.savefig(outpath, dpi=150)
            plt.close(fig)
            if DEBUG:
                print(f"plot_binned_efficiencies_all: Saved efficiency slope plot for {run_prefix}, xB bin {i_xB} to {outpath}")
#enddef