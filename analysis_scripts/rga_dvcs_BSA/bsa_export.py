# bsa_export.py

import os
import json
import math

def export_bsa_to_text(periods,
                       bin_means_json,
                       final_results_dir,
                       output_file):
    """
    periods            : list of DVCS‐period strings
    bin_means_json     : path to global bin‐means JSON
    final_results_dir  : path to 'final_results' dir (contains the adjusted_bsa_*.json and will hold the .txt)
    output_file        : filename for the combined .txt (e.g. "bsa_all_periods.txt")
    """
    # 1) load the per‐bin mean kinematics
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    # 2) ensure final_results folder exists
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # 3) fix our desired chronological order
    desired_order = [
        "DVCS_Sp18_out",
        "DVCS_Sp18_inb",
        "DVCS_Fa18_inb",
        "DVCS_Fa18_out",
        "DVCS_Sp19_inb",
    ]
    # only keep those present in `periods`
    period_list = [p for p in desired_order if p in periods]

    # 4) beam energies lookup
    beam_energies = {
        "DVCS_Fa18_inb":  10.604,
        "DVCS_Fa18_out":  10.604,
        "DVCS_Sp19_inb":  10.200,
        "DVCS_Sp18_inb":  10.594,
        "DVCS_Sp18_out":  10.594,
    }

    # helper to turn a key‐string into a tuple of ints
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").replace(" ", "").split(",")))

    # 5) write the single combined text file
    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")

        for period in period_list:
            Eb = beam_energies.get(period)
            if Eb is None:
                print(f"[export] no Eb defined for {period}, skipping")
                continue

            json_path = os.path.join(final_results_dir,
                                     f"adjusted_bsa_{period}.json")
            try:
                with open(json_path) as f:
                    adj = json.load(f)
            except FileNotFoundError:
                print(f"[export] missing {json_path}, skipping")
                continue

            # sort the bin‐keys by their integer tuples
            for key_str in sorted(adj.keys(), key=key_to_tuple):
                vals = adj[key_str]
                if not vals.get("valid", True):
                    continue
                mean = bin_means.get(key_str)
                if mean is None:
                    continue

                # convert phi from radians to degrees
                phi_deg = math.degrees(mean["phi_avg"])
                q2      = mean["Q2_avg"]
                xb      = mean["xB_avg"]
                t_val   = -mean["t_avg"]       # write –t
                A       = vals["bsa"]
                sigmaA  = vals["bsa_err"]

                # round everything to 3 decimals
                out.write(
                    f"{phi_deg:.3f} "
                    f"{q2:.3f} "
                    f"{xb:.3f} "
                    f"{t_val:.3f} "
                    f"{Eb:.3f} "
                    f"{A:.3f} "
                    f"{sigmaA:.3f}\n"
                )

    print(f"[export] wrote combined file → {output_path}")