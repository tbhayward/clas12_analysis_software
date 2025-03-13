#!/usr/bin/env python3

import os
import json
import numpy as np


def load_combined_bsa_json(json_filepath):
    """Loads combined BSA data from JSON, returning keys as integer tuples."""
    with open(json_filepath) as f:
        data = json.load(f)
        return {tuple(map(int, k.strip("()").split(", "))): v for k, v in data.items()}


def integrate_t_bins(input_json, output_json):
    combined_data = load_combined_bsa_json(input_json)

    # Use a dictionary to group data by (xB, Q2, phi) combinations
    integrated_data = {}

    # Loop over the original data and group by (xB, Q2, phi)
    bin_groups = {}
    for bin_key, values in combined_data.items():
        xB_idx, Q2_idx, _, phi_idx = bin_key
        key_3d = (xB_idx, Q2_idx, phi_idx)

        if key_3d not in bin_groups:
            bin_groups[key_3d] = []

        bin_groups[key_3d].append(values)

    # Now compute the combined BSA for each group
    for key_3d, measurements in bin_groups.items():
        bsa_vals = [m["bsa"] for m in measurements]
        bsa_errs = [m["bsa_err"] for m in measurements]

        weights = [1 / (err ** 2) for err in bsa_errs]
        total_weight = sum(weights)

        if total_weight > 0:
            combined_bsa = sum(w * val for w, val in zip(weights, bsa_vals)) / total_weight
            combined_err = np.sqrt(1 / total_weight)
        else:
            combined_bsa = np.mean(bsa_vals)
            combined_err = np.std(bsa_vals) if len(bsa_vals) > 1 else 0

        integrated_data[key_3d] = {
            "bsa": round(combined_bsa, 5),
            "bsa_err": round(combined_err, 5),
            "n_points": len(bin_groups[key_3d]),
            "valid": True
        }

    # Save results
    with open(output_json, 'w') as f:
        json.dump({str(k): v for k, v in integrated_data.items()}, f, indent=2)


if __name__ == "__main__":
    combined_bsa_json = os.path.join("final_results", "combined_bsa.json")
    output_integrated_json = os.path.join("final_results", "combined_bsa_integrated_t.json")
    integrate_t_bins(combined_bsa_json, output_integrated_json)