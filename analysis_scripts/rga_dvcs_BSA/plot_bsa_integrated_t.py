#!/usr/bin/env python3

import os
import json
import numpy as np

def load_combined_bsa_json(json_filepath):
    """Loads combined BSA data from JSON, returning keys as integer tuples."""
    with open(json_filepath) as f:
        data = json.load(f)
        combined_data = {
            tuple(map(int, k.strip("()").replace(" ", "").split(","))): v
            for k, v in data.items()
        }
    return combined_data

def integrate_t_bins(input_json, output_json):
    combined_data = load_combined_bsa(input_json)

    # Use a dictionary to group data by (xB, Q2, phi) combinations
    integrated_data = {}

    # Loop over the original data and group by (xB, Q2, phi)
    grouped_data = {}
    for (xB_idx, Q2_idx, t_idx, phi_idx), values in combined_data.items():
        key = (xB_idx, Q2_idx, phi_idx) = (xB, Q2, phi) = (xB_idx, Q2_idx, phi_idx) = (xB, Q2, phi) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = (xB_idx, Q2_idx, phi_idx) = bin_key[0], bin_key[1], bin_key[3]

    integrated_results = {}

    # Organize data by xB, Q2, phi (integrate out the t dimension)
    bin_groups = {}
    for bin_key, values in combined_data.items():
        xB_idx, Q2_idx, _, phi_idx = bin_key
        key_3d = (xB_idx, Q2_idx, phi_idx)
        if key not in bin_groups:
            bin_groups[key] = []
        bin_groups[key].append((values["bsa"], values["bsa_err"]))

    for key, measurements in bin_groups.items():
        bsa_vals, bsa_errs = zip(*measurements)

        # Propagate uncertainties via weighted average
        weights = [1 / (err ** 2) for err in bsa_errs]
        total_weight = sum(weights)

        if total_weight > 0:
            combined_bsa = sum(w * val for w, val in zip(weights, bsa_vals)) / total_weight
            combined_err = np.sqrt(1 / total_weight)
        else:
            combined_bsa = np.mean(bsa_vals)
            combined_err = np.std(bsa_vals) if len(bsa_vals) > 1 else bsa_errs[0]

        integrated_results[key] = {
            "bsa": round(combined_bsa, 5),
            "bsa_err": round(combined_err, 5),
            "n_points": len(bsa_vals),
            "valid": True
        }

    # Save results
    with open(output_json, 'w') as f:
        json.dump({str(k): v for k, v in integrated_results.items()}, f, indent=2)

    print(f"Integrated BSA (t-integrated) results saved to: {output_json}")

# Example usage if running as standalone script
if __name__ == "__main__":
    combined_bsa_json = os.path.join("final_results", "combined_bsa.json")
    output_integrated_json = os.path.join("final_results", "integrated_bsa.json")
    
    integrate_t_bins(combined_bsa_json, output_integrated_json)