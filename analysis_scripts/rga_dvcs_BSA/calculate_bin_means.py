# calculate_bin_means.py

import os
import json
import math
import numpy as np

# We assume you'll use the same structure for reading ROOT files 
# and applying cuts as in your other scripts.
from root_io import load_root_files
from kin_cuts import apply_kinematic_cuts, passes_3sigma_cuts
from load_binning_scheme import Binning

# Number of phi bins and their edges
N_PHI_BINS = 12
phi_edges = np.linspace(0, 2 * math.pi, N_PHI_BINS + 1)

def find_bin(value, bin_boundaries):
    """
    Returns the bin index (0-based) for 'value' given a list of (min, max) boundaries,
    or None if the value does not fall in any of the bins.
    """
    for i, (low, high) in enumerate(bin_boundaries):
        if low <= value < high:
            return i
        #endif
    #endfor
    return None
#endif

def calculate_bin_means(dvcs_periods, topologies, analysis_type, binning_scheme, output_json):
    """
    Calculates the mean xB, Q2, -t, and phi in each 4D bin for ALL given dvcs_periods 
    and topologies, combining them into a single global set of bin averages.
    
    Arguments:
      dvcs_periods  : list of strings (e.g. ["DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb"])
      topologies    : list of topologies (e.g. ["(FD,FD)", "(CD,FD)", "(CD,FT)"])
      analysis_type : typically "dvcs", but you can generalize
      binning_scheme: a list of Binning(xBmin, xBmax, Q2min, Q2max, tmin, tmax)
      output_json   : path to JSON file where global bin-averaged kinematics are saved

    The function will:
      1. Build one set of accumulators for every 4D bin.
      2. Loop over each period + each topology, load the DVCS data, apply cuts,
         and accumulate the total sums of xB, Q2, (-t), and phi in each bin.
      3. At the end, compute the mean for each bin and write to a single JSON.
    """

    print(f"[calculate_bin_means] Combining all periods: {dvcs_periods}")
    print(f"                    Combining topologies: {topologies}")

    # Build sets of unique bin boundaries from your binning scheme
    unique_xB_bins = sorted(set((b.xBmin, b.xBmax) for b in binning_scheme))
    unique_Q2_bins = sorted(set((b.Q2min, b.Q2max) for b in binning_scheme))
    unique_t_bins  = sorted(set((b.tmin, b.tmax) for b in binning_scheme))

    # 1) Prepare accumulators for sums and counts in each bin
    bin_accumulators = {}
    for i_xB in range(len(unique_xB_bins)):
        for i_Q2 in range(len(unique_Q2_bins)):
            for i_t in range(len(unique_t_bins)):
                for i_phi in range(N_PHI_BINS):
                    bin_accumulators[(i_xB, i_Q2, i_t, i_phi)] = {
                        'sum_xB':  0.0,
                        'sum_Q2':  0.0,
                        'sum_t':   0.0,  # We'll accumulate |t| (or -t) directly
                        'sum_phi': 0.0,
                        'count':   0
                    }
                #endfor
            #endfor
        #endfor
    #endfor

    # 2) Loop over each period and each topology, load the relevant data, apply cuts, accumulate sums
    for period in dvcs_periods:
        for topo in topologies:
            print(f"\n[calculate_bin_means] Processing {period}, topology={topo}")

            # Load DVCS trees
            _, dvcs_trees = load_root_files(period)
            if "data" not in dvcs_trees:
                raise KeyError(f"[ERROR] DVCS trees for period '{period}' do not contain 'data'.")

            # 3) Loop over DVCS data, apply cuts, and fill accumulators
            for event in dvcs_trees["data"]:
                try:
                    # Kinematic cuts
                    if not apply_kinematic_cuts(
                        event.t1, event.open_angle_ep2, event.theta_gamma_gamma,
                        event.Emiss2, event.Mx2, event.Mx2_1, event.Mx2_2,
                        event.pTmiss, event.xF,
                        analysis_type, "data", "", topo
                    ):
                        continue
                    #endif

                    # 3σ cuts (you can load the real dictionary if needed)
                    if not passes_3sigma_cuts(event, False, {}):
                        continue
                    #endif
                except Exception as e:
                    print(f"[calculate_bin_means] Exception during cuts on {period}: {e}")
                    continue
                #endtry

                # Sort the event into a 4D bin
                try:
                    xB_val   = float(event.x)
                    Q2_val   = float(event.Q2)
                    t_val    = abs(float(event.t1))  # if t1 is negative, use abs
                    phi_val  = float(event.phi2)
                    # phi_val  = float(event.Delta_phi)
                except Exception:
                    # Skip if we can't parse these
                    continue
                #endtry

                i_xB  = find_bin(xB_val,  unique_xB_bins)
                i_Q2  = find_bin(Q2_val,  unique_Q2_bins)
                i_t   = find_bin(t_val,   unique_t_bins)
                i_phi = np.digitize(phi_val, phi_edges) - 1

                # Make sure it actually fell in a valid bin
                if (i_xB is None) or (i_Q2 is None) or (i_t is None) or (i_phi < 0 or i_phi >= N_PHI_BINS):
                    continue
                #endif

                # Accumulate sums
                key = (i_xB, i_Q2, i_t, i_phi)
                bin_accumulators[key]['sum_xB']  += xB_val
                bin_accumulators[key]['sum_Q2']  += Q2_val
                bin_accumulators[key]['sum_t']   += t_val
                bin_accumulators[key]['sum_phi'] += phi_val
                bin_accumulators[key]['count']   += 1
            #endfor (dvcs_trees["data"])
        #endfor (topologies)
    #endfor (dvcs_periods)

    # 4) Compute means and store in a dictionary
    bin_means = {}
    for key, sums in bin_accumulators.items():
        count = sums['count']
        if count > 0:
            i_xB, i_Q2, i_t, i_phi = key
            xB_avg  = sums['sum_xB']  / count
            Q2_avg  = sums['sum_Q2']  / count
            t_avg   = sums['sum_t']   / count
            phi_avg = sums['sum_phi'] / count

            bin_means[str(key)] = {
                "xB_avg":  xB_avg,
                "Q2_avg":  Q2_avg,
                "t_avg":   t_avg,    # this is average of |t| if t1 < 0 in your data
                "phi_avg": phi_avg
            }
        #endif
    #endfor

    # 5) Write out to JSON
    with open(output_json, "w") as f:
        json.dump(bin_means, f, indent=2)

    print(f"\n[calculate_bin_means] Saved GLOBAL bin-averaged kinematics over all periods to {output_json}")
    return bin_means
#endif