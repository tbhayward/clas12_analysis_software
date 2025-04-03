#!/usr/bin/env python3
import os
import glob
import json
import ROOT
import numpy as np
import matplotlib.pyplot as plt
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts

def load_cuts(period, topology):
    """
    Loads the final cuts dictionary from combined_cuts.json for a given (period, topology).

    If the period includes 'bkg', it is mapped to the corresponding DVCS period.
    """
    topo_clean = topology.replace("(", "").replace(")", "").replace(",", "_").strip()
    if "bkg" in period:
        dvcs_equiv = period.replace("eppi0_bkg", "DVCS")
        dictionary_key = f"{dvcs_equiv}_{topo_clean}"
    else:
        dictionary_key = f"{period}_{topo_clean}"
    #endelse

    combined_cuts_path = os.path.join("exclusivity", "combined_cuts.json")
    if not os.path.exists(combined_cuts_path):
        print(f"⚠️ {combined_cuts_path} does not exist; returning empty cuts dictionary.")
        return {}
    #endif

    with open(combined_cuts_path, "r") as f:
        combined_cuts = json.load(f)
    #endwith

    if dictionary_key in combined_cuts:
        return combined_cuts[dictionary_key]
    else:
        print(f"⚠️ Key '{dictionary_key}' not found in {combined_cuts_path}.")
        print("Available keys:", list(combined_cuts.keys()))
        return {}
    #endif
#enddef

def passes_exclusivity_cuts(event, cuts, topology, analysis_type, run_period):
    """
    Applies the exclusivity cuts to an event based on the provided cuts dictionary.
    It first determines the event’s topology from detector1 and detector2,
    then applies kinematic cuts and finally the 3σ cuts.
    """
    # Determine event topology based on detector1 and detector2.
    if event.detector1 == 1 and event.detector2 == 1:
        event_topo = "(FD,FD)"
    elif event.detector1 == 2 and event.detector2 == 1:
        event_topo = "(CD,FD)"
    elif event.detector1 == 2 and event.detector2 == 0:
        event_topo = "(CD,FT)"
    else:
        return False
    #endif

    # Only proceed if event topology matches the selected topology.
    if event_topo != topology:
        return False
    #endif

    # Get theta value (for DVCS we expect theta_gamma_gamma).
    theta_val = getattr(event, "theta_gamma_gamma", None)

    # Apply kinematic cuts.
    if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, theta_val,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF,
            analysis_type, "data", run_period, topology
        ):
        return False
    #endif

    # Apply 3σ exclusivity cuts.
    if not passes_3sigma_cuts(event, False, cuts):
        return False
    #endif

    return True
#enddef

def get_tree_entries_with_cuts(file_path, cuts, topology, analysis_type, run_period):
    """
    Opens the ROOT file and counts the number of entries in the tree
    that pass the exclusivity cuts based on the provided cuts dictionary.
    Assumes the tree is named 'tree'. Adjust if necessary.
    """
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    #endif
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: tree not found in file {file_path}")
        f.Close()
        return 0
    #endif

    count = 0
    for event in tree:
        try:
            if passes_exclusivity_cuts(event, cuts, topology, analysis_type, run_period):
                count += 1
            #endif
        except Exception as e:
            print(f"Exception in applying exclusivity cuts: {e}")
            continue
        #endfor
    #endfor
    f.Close()
    return count
#enddef

def get_total_tree_entries(file_path):
    """
    Opens the ROOT file and returns the total number of entries in the tree.
    Assumes the tree is named 'tree'. Adjust if necessary.
    """
    f = ROOT.TFile.Open(file_path)
    if not f or f.IsZombie():
        print(f"Error opening file {file_path}")
        return 0
    #endif
    tree = f.Get("PhysicsEvents")
    if not tree:
        print(f"Error: tree not found in file {file_path}")
        f.Close()
        return 0
    #endif
    n = tree.GetEntries()
    f.Close()
    return n
#enddef

def plot_normalized_efficiencies(output_dir):
    """
    Calculates and plots normalized efficiencies as a function of beam current
    for each run period. For each current x nA, the efficiency is calculated as:

        efficiency = (n_reco(x nA)/n_gen(x nA)) / (n_reco(0 nA)/n_gen(0 nA))

    The results are then fitted to a linear polynomial.
    Three separate canvases (one per run period) are saved as PDFs.
    """
    # Define base directory where the ROOT files are stored.
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/efficiency_study"
    
    # Define the run periods with their file prefix mapping.
    run_periods = {
        "rga_fa18_inb": "DVCS_Fa18_inb",
        "rga_fa18_out": "DVCS_Fa18_out",
        "rga_sp19_inb": "DVCS_Sp19_inb"
    }
    
    # Define the selected topology and analysis type.
    selected_topology = "(FD,FD)"
    analysis_type = "dvcs"
    
    # Loop over each run period.
    for run_prefix, period_code in run_periods.items():
        # Find all generated and reconstructed ROOT files for this run period.
        pattern_gen = os.path.join(base_dir, f"gen_{run_prefix}_dvcs_*.root")
        pattern_reco = os.path.join(base_dir, f"{run_prefix}_dvcs_*.root")
        gen_files = glob.glob(pattern_gen)
        reco_files = glob.glob(pattern_reco)
        
        # Create dictionaries mapping current (in nA) to file path for both gen and reco.
        gen_dict = {}
        reco_dict = {}
        
        def extract_current(filename):
            """
            Extracts the beam current from the filename.
            Expected pattern: ..._{current}nA.root, where current is a number
            or 'nobkg' (which is mapped to 0).
            """
            base = os.path.basename(filename)
            parts = base.split("_")
            # Last part is expected to be like '45nA.root' or 'nobkg.root'
            last_part = parts[-1]
            last_part = last_part.replace(".root", "")
            if "nobkg" in last_part:
                return 0
            else:
                try:
                    current_str = last_part.replace("nA", "")
                    return int(current_str)
                except:
                    return None
                #endif
            #endif
        #enddef

        for fpath in gen_files:
            cur = extract_current(fpath)
            if cur is not None:
                gen_dict[cur] = fpath
            #endif
        #endfor

        for fpath in reco_files:
            cur = extract_current(fpath)
            if cur is not None:
                reco_dict[cur] = fpath
            #endif
        #endfor

        # Check that we have a 0 nA (nobkg) file for normalization.
        if 0 not in gen_dict or 0 not in reco_dict:
            print(f"Missing nobkg files for run period {run_prefix}. Skipping.")
            continue
        #endif
        
        # Load the exclusivity cuts for this period.
        cuts = load_cuts(period_code, selected_topology)
        
        # Calculate efficiency for each current.
        currents = []
        efficiencies = []
        efficiency_errors = []
        
        # Calculate normalization factor from the nobkg files (0 nA).
        n_gen_0 = get_total_tree_entries(gen_dict[0])
        n_reco_0 = get_tree_entries_with_cuts(reco_dict[0], cuts, selected_topology, analysis_type, period_code)
        if n_gen_0 == 0:
            print(f"Warning: 0 generated entries for normalization in run period {run_prefix}")
            norm_eff = 1.0
        else:
            norm_eff = (n_reco_0 / n_gen_0)
        #endif

        # Loop over all currents (ensuring 0 is included for normalization).
        for cur in sorted(gen_dict.keys()):
            if cur == 0:
                currents.append(0)
                efficiencies.append(1.0)  # by definition, normalized to itself
                efficiency_errors.append(0.0)
                continue
            #endif
            n_gen = get_total_tree_entries(gen_dict[cur])
            n_reco = get_tree_entries_with_cuts(reco_dict[cur], cuts, selected_topology, analysis_type, period_code)
            
            if n_gen == 0:
                eff = 0
                err = 0
            else:
                eff = (n_reco / n_gen) / norm_eff
                # Estimate error from Poisson statistics (assuming error in n_reco dominates).
                if n_reco > 0:
                    err = (np.sqrt(n_reco) / n_gen) / norm_eff
                else:
                    err = 0
                #endif
            #endif
            currents.append(cur)
            efficiencies.append(eff)
            efficiency_errors.append(err)
        #endfor
        
        # Convert lists to numpy arrays for fitting.
        currents_arr = np.array(currents)
        efficiencies_arr = np.array(efficiencies)
        weights = np.array([1/err if err > 0 else 1 for err in efficiency_errors])
        
        # Perform a linear fit: efficiency = a * current + b.
        coeffs, cov = np.polyfit(currents_arr, efficiencies_arr, 1, w=weights, cov=True)
        a, b = coeffs
        
        # Generate fitted line data.
        fit_x = np.linspace(min(currents_arr), max(currents_arr), 100)
        fit_y = a * fit_x + b
        
        # Create the plot.
        plt.figure()
        plt.errorbar(currents_arr, efficiencies_arr, yerr=efficiency_errors, fmt='ko', label='Data')  # black points with error bars
        plt.plot(fit_x, fit_y, 'r--', label=f'Fit: y = {a:.3e}x + {b:.3e}')  # dashed red line
        plt.xlabel("current (nA)")
        plt.ylabel("normalized efficiency")
        plt.legend(loc='upper right')
        plt.title(f"Normalized Efficiency for {run_prefix}")
        
        # Save the plot as a PDF.
        output_path = os.path.join(output_dir, f"{run_prefix}_integrated.pdf")
        plt.savefig(output_path, format='pdf')
        plt.close()
        print(f"Saved plot for {run_prefix} to {output_path}")
    #endfor
#enddef