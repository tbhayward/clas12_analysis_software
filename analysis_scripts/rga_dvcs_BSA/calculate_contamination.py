import os
import json
import math
import numpy as np
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from root_io import load_root_files

# Define the number of φ bins and their edges.
N_PHI_BINS = 12
phi_edges = np.linspace(0, 2 * math.pi, N_PHI_BINS + 1)

def find_bin(value, bin_boundaries):
    """
    Given a value and a list of (min, max) boundaries, return the bin index (0-indexed)
    if the value falls in a bin; otherwise, return None.
    """
    for i, (low, high) in enumerate(bin_boundaries):
        if low <= value < high:
            return i
    return None

import os
import json

def load_cuts(period, topology):
    """
    Loads the final cuts dictionary from combined_cuts.json for a given (period, topology).

    If the period includes 'bkg' (e.g. 'eppi0_bkg_Sp19_inb'), we map it to the corresponding
    'DVCS' period (e.g. 'DVCS_Sp19_inb') so that the same 3σ cuts are used.
    """
    # Clean up the topology string for the JSON key.
    topo_clean = topology.replace("(", "").replace(")", "").replace(",", "_").strip()

    # If this period is a 'bkg' sample, map it to the corresponding DVCS period.
    if "bkg" in period:
        # Example: "eppi0_bkg_Sp19_inb" -> "DVCS_Sp19_inb"
        dvcs_equiv = period.replace("eppi0_bkg", "DVCS")
        dictionary_key = f"{dvcs_equiv}_{topo_clean}"
    else:
        dictionary_key = f"{period}_{topo_clean}"

    # Load combined_cuts.json
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

    if not cuts_dict:
        print(f"[WARNING] No cuts found for key '{dictionary_key}'. Is combined_cuts.json correct?")


def calculate_contamination(period, topology, analysis_type, binning_scheme):
    """
    Calculate the 4D contamination for a given DVCS period and topology.
    
    This function loads the DVCS trees, and the corresponding π⁰ trees (from
    eppi0 and eppi0_bkg samples), applies the kinematic and 3σ cuts, fills 4D bins,
    and then computes the contamination.
    
    Returns a dictionary keyed by a 4‐tuple (i_xB, i_Q2, i_t, i_phi) with contamination info.
    """
    print(f"[calculate_contamination] Beginning contamination calculation for {period}, topology {topology}, analysis {analysis_type}")
    print(f"DEBUG: Looking for key: {period}_{topology}")

    _, dvcs_trees = load_root_files(period)

    # Check that the 'data' key exists.
    if "data" not in dvcs_trees:
        raise KeyError(f"[ERROR] DVCS trees for period {period} do not contain 'data'.")

    # For DVCS, we use:
    #   DVCS data: dvcs_trees["data"]
    # For π⁰ samples we need to build the corresponding period names.
    pi0_exp_period = period.replace("DVCS", "eppi0")
    pi0_reco_period = period.replace("DVCS", "eppi0")
    pi0_bkg_period  = period.replace("DVCS", "eppi0_bkg")
    
    _, pi0_exp_trees = load_root_files(pi0_exp_period)
    _, pi0_reco_trees = load_root_files(pi0_reco_period)
    pi0_bkg_trees = {}
    try:
        _, pi0_bkg_trees = load_root_files(pi0_bkg_period)
    except ValueError:
        print(f"[WARNING] No background MC found for {pi0_bkg_period}, skipping.")
        pi0_bkg_trees = {}

    _, pi0_bkg_trees = load_root_files(pi0_bkg_period)
    
    # Load the cuts dictionary (using DVCS cuts even for bkg, if needed).
    cuts_dict = load_cuts(period, topology)

    # In calculate_contamination.py, add to topology check:
    fd_fd_count = cd_fd_count = cd_ft_count = 0
    for event in pi0_bkg_trees["mc"]:
        if event.detector1 == 1 and event.detector2 == 1:
            fd_fd_count +=1
        elif event.detector1 == 2 and event.detector2 == 1:
            cd_fd_count +=1
        elif event.detector1 == 2 and event.detector2 == 0:
            cd_ft_count +=1

    print(f"\n[TOPOLOGY] {pi0_bkg_period} raw MC counts:")
    print(f"  (FD,FD): {fd_fd_count}")
    print(f"  (CD,FD): {cd_fd_count}")
    print(f"  (CD,FT): {cd_ft_count}")

    # After loading cuts_dict in calculate_contamination.py:
    print(f"\n=== Cuts for {period} {topology} ===")
    for var, vals in cuts_dict.get("data", {}).items():
        print(f"{var}: μ={vals['mean']:.3f} ±3σ→({vals['mean']-3*vals['std']:.3f}, {vals['mean']+3*vals['std']:.3f})")
    
    # Build bin boundaries from the binning scheme.
    xB_bins = [(b.xBmin, b.xBmax) for b in binning_scheme]
    Q2_bins = [(b.Q2min, b.Q2max) for b in binning_scheme]
    t_bins  = [(b.tmin, b.tmax) for b in binning_scheme]
    

    # Extract the unique bin boundaries from the binning scheme.
    unique_xB_bins = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    unique_Q2_bins = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    unique_t_bins  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))
    # Initialize a results dictionary for every 4D bin.
    results = {}
    for i_xB in range(len(unique_xB_bins)):
        for i_Q2 in range(len(unique_Q2_bins)):
            for i_t in range(len(unique_t_bins)):
                for i_phi in range(N_PHI_BINS):
                    results[(i_xB, i_Q2, i_t, i_phi)] = {
                        'N_data': 0,
                        'N_pi0_mc': 0,
                        'N_pi0_exp': 0,
                        'N_pi0_reco': 0
                    }

    
    # --- Count DVCS data events ---
    count = 0
    for event in dvcs_trees["data"]:
        # if count >= 10000:
        #     break
        # count += 1
        try:
            if not apply_kinematic_cuts(
                event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
                event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                event.pTmiss, event.xF,
                analysis_type, "data", "", topology
            ):
                continue
            if not passes_3sigma_cuts(event, False, cuts_dict):
                continue
        except Exception as e:
            print(f"[calculate_contamination] Exception in DVCS data cut: {e}")
            continue

        try:
            xB_val = float(event.x)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t1))
            phi_val = float(event.phi2)
            # phi_val = float(event.phi1)
        except Exception as e:
            continue
        i_xB = find_bin(xB_val, unique_xB_bins)
        i_Q2 = find_bin(Q2_val, unique_Q2_bins)
        i_t  = find_bin(t_val, unique_t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_data'] += 1



    # --- Count π⁰ misidentification events from eppi0_bkg MC ---
    # In the π⁰ MC event loop (pi0_bkg_trees["mc"] section):
    total_events = 0
    passed_kinematic = 0
    passed_3sigma = 0
    for event in pi0_bkg_trees["mc"]:
        total_events += 1
        # if count >= 10000:
        #     break
        # count += 1
        try:
            if not apply_kinematic_cuts(
                event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
                event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                event.pTmiss, event.xF,
                analysis_type, "mc", "", topology
            ):
                continue
            passed_kinematic += 1
            if not passes_3sigma_cuts(event, True, cuts_dict):
                continue
            passed_3sigma += 1
        except Exception as e:
            continue

        try:
            xB_val = float(event.x)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t1))
            phi_val = float(event.phi2)
            # phi_val = float(event.phi1)
        except Exception as e:
            continue
        i_xB = find_bin(xB_val, unique_xB_bins)
        i_Q2 = find_bin(Q2_val, unique_Q2_bins)
        i_t  = find_bin(t_val, unique_t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_pi0_mc'] += 1

    print(f"[CUTS] {pi0_bkg_period} {topology}:")
    print(f"  Total MC: {total_events}")
    print(f"  Passed kinematic: {passed_kinematic} ({passed_kinematic/total_events:.1%})")
    print(f"  Passed 3σ: {passed_3sigma} ({passed_3sigma/passed_kinematic:.1%})")

    # --- Count π⁰ experimental events from eppi0 data ---
    count = 0
    for event in pi0_exp_trees.get("data", []):
        # if count >= 10000:
        #     break
        # count += 1
        try:
            if not apply_kinematic_cuts(
                event.t1, event.open_angle_ep2, event.theta_pi0_pi0,
                event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                event.pTmiss, event.xF,
                "eppi0", "data", "", topology
            ):
                continue
            if not passes_3sigma_cuts(event, False, cuts_dict):
                continue
        except Exception as e:
            continue

        try:
            xB_val = float(event.x)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t1))
            phi_val = float(event.phi2)
            # phi_val = float(event.phi1)
        except Exception as e:
            continue
        i_xB = find_bin(xB_val, unique_xB_bins)
        i_Q2 = find_bin(Q2_val, unique_Q2_bins)
        i_t  = find_bin(t_val, unique_t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_pi0_exp'] += 1

    # --- Count π⁰ reconstructed events from eppi0 MC ---
    count = 0
    for event in pi0_exp_trees.get("mc", []):
        # if count >= 10000:
        #     break
        # count += 1
        try:
            if not apply_kinematic_cuts(
                event.t1, event.open_angle_ep2, event.theta_pi0_pi0,
                event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                event.pTmiss, event.xF,
                "eppi0", "mc", "", topology
            ):
                continue
            if not passes_3sigma_cuts(event, True, cuts_dict):
                continue
        except Exception as e:
            continue

        try:
            xB_val = float(event.x)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t1))
            phi_val = float(event.phi2)
            # phi_val = float(event.phi1)
        except Exception as e:
            continue
        i_xB = find_bin(xB_val, unique_xB_bins)
        i_Q2 = find_bin(Q2_val, unique_Q2_bins)
        i_t  = find_bin(t_val, unique_t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_pi0_reco'] += 1

    # --- Compute contamination in each 4D bin ---
    for key, counts in results.items():
        N_data = counts['N_data']
        if N_data == 0 or counts['N_pi0_reco'] == 0:
            counts['c_i'] = 0.0
            counts['c_i_err'] = 0.0
        else:
            ratio = counts['N_pi0_exp'] / counts['N_pi0_reco']
            c_i = counts['N_pi0_mc'] * ratio / N_data
            rel_pi0_mc = 1 / math.sqrt(counts['N_pi0_mc']) if counts['N_pi0_mc'] > 0 else 0
            rel_pi0_exp = 1 / math.sqrt(counts['N_pi0_exp']) if counts['N_pi0_exp'] > 0 else 0
            rel_pi0_reco = 1 / math.sqrt(counts['N_pi0_reco']) if counts['N_pi0_reco'] > 0 else 0
            rel_data = 1 / math.sqrt(N_data)
            # For the ratio (N_pi0_exp/N_pi0_reco) we sum the relative uncertainties.
            rel_ratio = math.sqrt(rel_pi0_exp**2 + rel_pi0_reco**2)
            rel_err = math.sqrt(rel_pi0_mc**2 + rel_ratio**2 + rel_data**2)
            c_i_err = c_i * rel_err
            counts['c_i'] = c_i
            counts['c_i_err'] = c_i_err

    return results