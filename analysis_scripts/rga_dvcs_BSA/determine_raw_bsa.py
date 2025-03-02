import os
import json
import math
import numpy as np
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from root_io import load_root_files
from load_binning_scheme import load_binning_scheme

def determine_topology(event):
    """Determine topology based on detector1/detector2 values"""
    if event.detector1 == 1 and event.detector2 == 1:
        return "(FD,FD)"
    elif event.detector1 == 2 and event.detector2 == 1:
        return "(CD,FD)"
    elif event.detector1 == 2 and event.detector2 == 0:
        return "(CD,FT)"
    return None

def load_cuts(period, topology):
    """Load cuts from combined_cuts.json with proper topology handling"""
    topo_clean = topology.replace("(", "").replace(")", "").replace(",", "_").strip()
    
    if "bkg" in period:
        dvcs_period = period.replace("eppi0_bkg", "DVCS")
        json_key = f"{dvcs_period}_{topo_clean}"
    else:
        json_key = f"{period}_{topo_clean}"

    cuts_path = os.path.join("exclusivity", "combined_cuts.json")
    with open(cuts_path) as f:
        return json.load(f).get(json_key, {})

def calculate_raw_bsa(period, channel, binning_csv, output_dir):
    """Calculate BSA with beam polarization scaling and validity flags"""
    binning_scheme = load_binning_scheme(binning_csv)
    unique_xB = sorted({(b.xBmin, b.xBmax) for b in binning_scheme})
    unique_Q2 = sorted({(b.Q2min, b.Q2max) for b in binning_scheme})
    unique_t = sorted({(b.tmin, b.tmax) for b in binning_scheme})
    phi_bins = np.linspace(0, 2*math.pi, 13)

    # Initialize storage for counts and beam polarization
    results = {(i_xB, i_Q2, i_t, i_phi): {"N_plus":0, "N_minus":0} 
              for i_xB in range(len(unique_xB)) for i_Q2 in range(len(unique_Q2))
              for i_t in range(len(unique_t)) for i_phi in range(12)}
    
    beam_pol_sum = 0.0
    beam_pol_count = 0

    try:
        _, trees = load_root_files(period)
        data_tree = trees.get("data", [])
    except Exception as e:
        print(f"Error loading {period}: {e}")
        return {}

    # First pass: collect beam polarization data and counts
    for event in data_tree:
        topology = determine_topology(event)
        if not topology: continue
        
        cuts_dict = load_cuts(period, topology)
        if not cuts_dict: continue

        theta_var = event.theta_gamma_gamma if channel == "dvcs" else event.theta_pi0_pi0
        if not apply_kinematic_cuts(
            event.t1, event.open_angle_ep2, theta_var,
            event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
            event.pTmiss, event.xF, "dvcs" if channel == "dvcs" else "eppi0", 
            "data", "", topology
        ) or not passes_3sigma_cuts(event, False, cuts_dict):
            continue

        # try:
            # Collect beam polarization data
        beam_pol = event.beam_pol
        beam_pol_sum += beam_pol
        beam_pol_count += 1

        # Bin the event
        xB, Q2, t, phi = event.x, event.Q2, abs(event.t1), event.phi2
        i_xB = next(i for i, (lo,hi) in enumerate(unique_xB) if lo <= xB < hi)
        i_Q2 = next(i for i, (lo,hi) in enumerate(unique_Q2) if lo <= Q2 < hi)
        i_t = next(i for i, (lo,hi) in enumerate(unique_t) if lo <= t < hi)
        i_phi = np.digitize(phi, phi_bins) - 1
        
        if 0 <= i_phi < 12:
            key = (i_xB, i_Q2, i_t, i_phi)
            if event.helicity > 0:
                results[key]["N_plus"] += 1
            else:
                results[key]["N_minus"] += 1
        # except Exception as e:
        #     print(f"Error processing event: {e}")
        #     continue

    # Calculate average beam polarization
    if beam_pol_count == 0:
        print(f"No valid events found for {period} {channel}")
        return {}
        
    beam_pol_avg = beam_pol_sum / beam_pol_count
    print(f"Calculated beam polarization for {period} {channel}: {beam_pol_avg:.3f}")

    # Second pass: calculate scaled BSA values
    valid_results = {}
    for key, counts in results.items():
        total = counts["N_plus"] + counts["N_minus"]
        if total == 0: continue
        
        # Calculate raw BSA and error
        bsa_raw = (counts["N_plus"] - counts["N_minus"]) / total
        bsa_err_raw = math.sqrt((4 * counts["N_plus"] * counts["N_minus"]) / (total**3))
        
        # Apply beam polarization scaling
        try:
            bsa = bsa_raw / beam_pol_avg
            bsa_err = bsa_err_raw / beam_pol_avg
        except ZeroDivisionError:
            print(f"Invalid beam polarization {beam_pol_avg} for {period}")
            bsa = 0
            bsa_err = 0

        valid_results[key] = {
            "bsa": round(bsa, 5),
            "bsa_err": round(bsa_err, 5),
            "valid": True,
            "beam_pol": round(beam_pol_avg, 3),
            "N_total": total
        }

    # Save results
    os.makedirs(output_dir, exist_ok=True)
    out_path = os.path.join(output_dir, f"raw_bsa_{channel}_{period}.json")
    with open(out_path, "w") as f:
        json.dump({str(k): v for k, v in valid_results.items()}, f, indent=2)
        
    return valid_results