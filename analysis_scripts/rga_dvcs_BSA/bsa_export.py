# bsa_export.py

import os
import json
import math

def export_bsa_to_text(periods,
                       bin_means_json,
                       final_results_dir,
                       output_file):
    """
    (unchanged)
    Writes one line per period+bin into final_results/output_file
    """
    # … your existing code here, unchanged …

def export_bsa_grouped_to_text(bin_means_json,
                               final_results_dir,
                               output_file):
    """
    Combines:
      - SP18_out & SP18_inb → “Sp18”
      - FA18_inb & FA18_out → “Fa18”
      - SP19_inb           → “Sp19”

    For each bin‐index present in all periods of a group, does a
    weighted average of A (and σA) with weights = 1/σA², and uses the
    same global φ, Q², xB, t from bin_means_json.  Writes one file
    final_results/output_file with three blocks in order:
      Sp18 combined
      Fa18 combined
      Sp19_inb
    """
    # 1) load bin‐means
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    # 2) ensure output dir
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # 3) grouping definitions in the order you requested
    group_defs = [
        ("Sp18", ["DVCS_Sp18_out", "DVCS_Sp18_inb"]),
        ("Fa18", ["DVCS_Fa18_inb", "DVCS_Fa18_out"]),
        ("Sp19", ["DVCS_Sp19_inb"]),
    ]

    # 4) beam energies (they’re identical within each pair)
    beam_energies = {
        "DVCS_Fa18_inb": 10.604,
        "DVCS_Fa18_out": 10.604,
        "DVCS_Sp19_inb": 10.200,
        "DVCS_Sp18_inb": 10.594,
        "DVCS_Sp18_out": 10.594,
    }

    # helper to sort the JSON keys
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").replace(" ", "").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")

        for group_name, periods in group_defs:
            # load each period’s adjusted_bsa JSON
            adjs = []
            for p in periods:
                path = os.path.join(final_results_dir,
                                    f"adjusted_bsa_{p}.json")
                try:
                    with open(path) as f:
                        adjs.append(json.load(f))
                except FileNotFoundError:
                    print(f"[grouped export] missing {path}, skipping group {group_name}")
                    adjs = []
                    break
            if not adjs:
                continue

            # find the keys common to all periods in this group
            common = set(adjs[0].keys())
            for d in adjs[1:]:
                common &= set(d.keys())
            if not common:
                continue

            # pick Eb from the first period
            Eb = beam_energies[periods[0]]

            # separator/comment for this block
            out.write(f"\n# {group_name} combined\n")

            for key_str in sorted(common, key=key_to_tuple):
                # gather A_i, σA_i for each period
                As     = []
                sigAs  = []
                valid  = True
                for d in adjs:
                    v = d[key_str]
                    if not v.get("valid", True):
                        valid = False
                        break
                    As.append(v["bsa"])
                    sigAs.append(v["bsa_err"])
                if not valid:
                    continue

                # weighted average of A
                weights = [1.0/(s**2) for s in sigAs]
                Wsum    = sum(weights)
                A_avg   = sum(w*a for w,a in zip(weights, As)) / Wsum
                sigA    = math.sqrt(1.0/Wsum)

                # kinematics from global bin_means
                m = bin_means.get(key_str)
                if m is None:
                    continue
                phi_deg = math.degrees(m["phi_avg"])
                q2      = m["Q2_avg"]
                xb      = m["xB_avg"]
                t_val   = -m["t_avg"]

                # write rounded to 3 decimals
                out.write(
                    f"{phi_deg:.3f} "
                    f"{q2:.3f} "
                    f"{xb:.3f} "
                    f"{t_val:.3f} "
                    f"{Eb:.3f} "
                    f"{A_avg:.3f} "
                    f"{sigA:.3f}\n"
                )

    print(f"[grouped export] wrote → {output_path}")