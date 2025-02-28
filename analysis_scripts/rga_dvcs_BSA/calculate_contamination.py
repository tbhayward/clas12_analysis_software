#!/usr/bin/env python3
"""
Module: calculate_contamination.py

This module defines a function to compute the π⁰ contamination (cᵢ) in DVCS analyses
using a 4D binning scheme:
  - xB, Q², |t|: boundaries provided by the CSV (loaded via load_binning_scheme)
  - φ: 12 equally spaced bins between 0 and 2π

For each 4D bin the contamination is defined as:

    cᵢ = (N_pi0_mc) * [ (N_pi0_exp) / (N_pi0_reco) ] / (N_data)

where
  - N_data is the DVCS data event count (from DVCS_* data),
  - N_pi0_mc is the π⁰ misidentification count (from eppi0_bkg_* MC),
  - N_pi0_exp is the experimental π⁰ count (from eppi0_* data),
  - N_pi0_reco is the reconstructed π⁰ count (from eppi0_* MC).

The final cuts (e.g. mu ± 2.5σ) are read in from the combined cuts JSON produced by exclusivity.
Uncertainties are propagated assuming Poisson statistics.
"""

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

def load_cuts(period, topology):
    """
    Loads the final cuts dictionary from combined_cuts.json for a given (period, topology).

    If the period includes 'bkg' (e.g. 'eppi0_bkg_Sp19_inb'), we REPLACE that
    with the corresponding 'DVCS' period (so we reuse the DVCS dictionary).
    That way, the bkg sample uses the same 3σ cuts as DVCS, avoiding KeyError('data').
    """
    # 1) Clean up the topology string for the JSON key.
    #    We remove parentheses, commas, AND all internal spaces:
    topo_clean = (topology
                  .replace("(", "")
                  .replace(")", "")
                  .replace(",", "_")
                  .replace(" ", ""))  # Remove ALL spaces

    # 2) If this period is a 'bkg' sample, re-map it to the DVCS version:
    #    e.g. "eppi0_bkg_Sp19_inb" -> "DVCS_Sp19_inb"
    if "bkg" in period:
        dvcs_equiv = period.replace("eppi0_bkg", "DVCS")  
        # or whatever transformation your naming scheme requires
        dictionary_key = f"{dvcs_equiv}_{topo_clean}"
    else:
        # Normal case
        dictionary_key = f"{period}_{topo_clean}"

    # Debug print to confirm exactly what we are looking for:
    print(f"DEBUG load_cuts => looking up key: '{dictionary_key}'")

    # 3) Load the combined_cuts.json
    combined_cuts_path = os.path.join("exclusivity", "combined_cuts.json")
    cuts_dict = {}

    if os.path.exists(combined_cuts_path):
        with open(combined_cuts_path, "r") as f:
            combined_cuts = json.load(f)

        if dictionary_key in combined_cuts:
            cuts_dict = combined_cuts[dictionary_key]
        else:
            # If the key isn't found, we return an empty dict => 3σ cuts will fail
            print(f"⚠️ Key '{dictionary_key}' not found in {combined_cuts_path}; "
                  f"available keys: {list(combined_cuts.keys())}\n"
                  f"Returning empty cuts_dict.")
    else:
        print(f"⚠️ {combined_cuts_path} does not exist, returning empty cuts_dict.")

    return cuts_dict

def calculate_contamination(period, topology, analysis_type, binning_scheme):
    """
    Calculate the 4D contamination in the DVCS analysis for a given period.

    Parameters:
      period: DVCS period (e.g., "DVCS_Fa18_inb")
      topology: a string representing the event topology, e.g. "(FD,FD)"
      analysis_type: should be "dvcs"
      binning_scheme: list of namedtuples (loaded from the CSV via load_binning_scheme)
                      with fields: xBmin, xBmax, Q2min, Q2max, tmin, tmax

    Returns:
      A dictionary with keys being a 4-tuple of bin indices: (i_xB, i_Q2, i_t, i_phi)
      and values a dictionary containing:
          {
              'N_data':     <int>,   # DVCS data count
              'N_pi0_mc':   <int>,   # π⁰ misidentification count (from eppi0_bkg MC)
              'N_pi0_exp':  <int>,   # π⁰ experimental count (from eppi0 data)
              'N_pi0_reco': <int>,   # π⁰ reconstructed count (from eppi0 MC)
              'c_i':        <float>, # contamination in this bin
              'c_i_err':    <float>  # uncertainty on contamination
          }
    """
    print("beginning of calculate_contamination")
    # Load the DVCS trees.
    _, dvcs_trees = load_root_files(period)
    
    # For DVCS, use:
    #   Data: from DVCS_* data.
    # Now, for the π⁰ quantities we use fixed mapping:
    #   Experimental π⁰: from eppi0_* data (period with "DVCS" replaced by "eppi0")
    #   Reconstructed π⁰: from eppi0_* MC (same as experimental sample for this example)
    #   Misidentified π⁰: from eppi0_bkg_* MC (period with "DVCS" replaced by "eppi0_bkg")
    pi0_exp_period = period.replace("DVCS", "eppi0")
    pi0_reco_period = period.replace("DVCS", "eppi0")
    pi0_bkg_period  = period.replace("DVCS", "eppi0_bkg")
    
    _, pi0_exp_trees = load_root_files(pi0_exp_period)
    _, pi0_reco_trees = load_root_files(pi0_reco_period)
    _, pi0_bkg_trees  = load_root_files(pi0_bkg_period)


    # Load the cuts from the exclusivity processing.
    cuts_dict = load_cuts(period, topology)
    
    # Build bin boundaries for xB, Q², and t from the binning scheme.
    xB_bins = [(b.xBmin, b.xBmax) for b in binning_scheme]
    Q2_bins = [(b.Q2min, b.Q2max) for b in binning_scheme]
    t_bins  = [(b.tmin, b.tmax) for b in binning_scheme]
    
    # Initialize the results dictionary for every 4D bin.
    results = {}
    for i_xB in range(len(xB_bins)):
        for i_Q2 in range(len(Q2_bins)):
            for i_t in range(len(t_bins)):
                for i_phi in range(N_PHI_BINS):
                    results[(i_xB, i_Q2, i_t, i_phi)] = {
                        'N_data': 0,
                        'N_pi0_mc': 0,
                        'N_pi0_exp': 0,
                        'N_pi0_reco': 0
                    }
    
    # --- Count DVCS data events ---
    for event in dvcs_trees["data"]:
        if not apply_kinematic_cuts(event.t1, event.open_angle_ep2, 0.0,
                                    event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                                    event.pTmiss, event.xF,
                                    analysis_type, "data", "", topology):
            continue
        if not passes_3sigma_cuts(event, False, cuts_dict):
            continue
        try:
            xB_val = float(event.xB)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t))
            phi_val = float(event.phi2)
        except AttributeError:
            continue
        i_xB = find_bin(xB_val, xB_bins)
        i_Q2 = find_bin(Q2_val, Q2_bins)
        i_t  = find_bin(t_val, t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_data'] += 1

    # --- Count π⁰ misidentification events from eppi0_bkg MC ---
    for event in pi0_bkg_trees["mc"]:
        if not apply_kinematic_cuts(event.t1, event.open_angle_ep2, 0.0,
                                    event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                                    event.pTmiss, event.xF,
                                    analysis_type, "mc", "", topology):
            continue
        if not passes_3sigma_cuts(event, True, cuts_dict):
            continue
        try:
            xB_val = float(event.xB)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t))
            phi_val = float(event.phi2)
        except AttributeError:
            continue
        i_xB = find_bin(xB_val, xB_bins)
        i_Q2 = find_bin(Q2_val, Q2_bins)
        i_t  = find_bin(t_val, t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_pi0_mc'] += 1

    # --- Count π⁰ experimental events from eppi0 data ---
    for event in pi0_exp_trees["data"]:
        if not apply_kinematic_cuts(event.t1, event.open_angle_ep2, 0.0,
                                    event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                                    event.pTmiss, event.xF,
                                    "eppi0", "data", "", topology):
            continue
        if not passes_3sigma_cuts(event, False, cuts_dict):
            continue
        try:
            xB_val = float(event.xB)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t))
            phi_val = float(event.phi2)
        except AttributeError:
            continue
        i_xB = find_bin(xB_val, xB_bins)
        i_Q2 = find_bin(Q2_val, Q2_bins)
        i_t  = find_bin(t_val, t_bins)
        i_phi = np.digitize(phi_val, phi_edges) - 1
        if i_xB is None or i_Q2 is None or i_t is None or i_phi is None or i_phi < 0 or i_phi >= N_PHI_BINS:
            continue
        results[(i_xB, i_Q2, i_t, i_phi)]['N_pi0_exp'] += 1

    # --- Count π⁰ reconstructed events from eppi0 MC ---
    for event in pi0_exp_trees["mc"]:
        if not apply_kinematic_cuts(event.t1, event.open_angle_ep2, 0.0,
                                    event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                                    event.pTmiss, event.xF,
                                    "eppi0", "mc", "", topology):
            continue
        if not passes_3sigma_cuts(event, True, cuts_dict):
            continue
        try:
            xB_val = float(event.xB)
            Q2_val = float(event.Q2)
            t_val  = abs(float(event.t))
            phi_val = float(event.phi2)
        except AttributeError:
            continue
        i_xB = find_bin(xB_val, xB_bins)
        i_Q2 = find_bin(Q2_val, Q2_bins)
        i_t  = find_bin(t_val, t_bins)
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
            rel_ratio = rel_pi0_exp + rel_pi0_reco
            rel_err = math.sqrt(rel_pi0_mc**2 + rel_ratio**2 + rel_data**2)
            c_i_err = c_i * rel_err
            counts['c_i'] = c_i
            counts['c_i_err'] = c_i_err

    return results