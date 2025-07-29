# bsa_export.py

import os
import json

def export_bsa_to_text(periods,
                       bin_means_json,
                       final_results_dir,
                       output_file):
    """
    periods            : list of strings, e.g. ["DVCS_Fa18_inb", …]
    bin_means_json     : path to the global bin‐means JSON
    final_results_dir  : directory containing adjusted_bsa_{period}.json
    output_file        : full path to the single .txt you want written

    Produces one text file with a header and then one line per (period,bin):
      # phi(rad) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA
    """
    # load the per‐bin mean kinematics
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    # beam energies by period
    beam_energies = {
        "DVCS_Fa18_inb": 10.604,
        "DVCS_Fa18_out": 10.604,
        "DVCS_Sp19_inb": 10.200,
        "DVCS_Sp18_inb": 10.594,
        "DVCS_Sp18_out": 10.594,
    }

    # ensure parent folder exists
    parent = os.path.dirname(output_file)
    if parent:
        os.makedirs(parent, exist_ok=True)

    with open(output_file, "w") as out:
        out.write("# phi(rad) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")

        for period in periods:
            Eb = beam_energies.get(period)
            if Eb is None:
                print(f"[export] no Eb defined for {period}, skipping")
                continue

            # load that period's adjusted BSA
            fname = f"adjusted_bsa_{period}.json"
            path  = os.path.join(final_results_dir, fname)
            try:
                with open(path) as f:
                    adj = json.load(f)
            except FileNotFoundError:
                print(f"[export] missing {path}, skipping")
                continue

            for key, vals in adj.items():
                if not vals.get("valid", True):
                    continue
                mean = bin_means.get(key)
                if mean is None:
                    continue

                φ   = mean["phi_avg"]
                Q2  = mean["Q2_avg"]
                xB  = mean["xB_avg"]
                t   = -mean["t_avg"]        # write –t
                A   = vals["bsa"]
                σA  = vals["bsa_err"]

                out.write(f"{φ} {Q2} {xB} {t} {Eb:.3f} {A} {σA}\n")

    print(f"[export] wrote combined file → {output_file}")