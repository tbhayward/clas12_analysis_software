# bsa_export.py

import os
import json
import math


def export_bsa_to_text(periods,
                       bin_means_json,
                       final_results_dir,
                       output_file):
    """
    periods            : list of DVCS-period strings
    bin_means_json     : path to global bin-means JSON
    final_results_dir  : path to 'final_results' dir (contains the adjusted_bsa_*.json and will hold the .txt)
    output_file        : filename for the combined .txt (e.g. "bsa_all_periods.txt")
    """
    # 1) load the per-bin mean kinematics
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

    # helper to turn a key-string into a tuple of ints
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").replace(" ", "").split(",")))

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

            # sort the bin-keys by their integer tuples
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
                t_val   = -mean["t_avg"]       # write -t
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


def export_bsa_grouped_to_text(bin_means_json,
                               final_results_dir,
                               output_file):
    """
    Writes one combined text file with no extra headers.  Combines:
      - SP18_out & SP18_inb → "Sp18"
      - FA18_inb & FA18_out → "Fa18"
      - SP19_inb           → "Sp19"

    Uses weighted averages of A with weights=1/σA², if σA=0 only the periods with σA>0 are used.
    """
    # 1) load bin-means
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    # 2) ensure final_results folder exists
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # grouping definitions
    group_defs = [
        ["DVCS_Sp18_out", "DVCS_Sp18_inb"],
        ["DVCS_Fa18_inb", "DVCS_Fa18_out"],
        ["DVCS_Sp19_inb"],
    ]

    # beam energies
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

        for periods in group_defs:
            # load each period’s adjusted_bsa JSON
            adjs = []
            for p in periods:
                path = os.path.join(final_results_dir,
                                    f"adjusted_bsa_{p}.json")
                try:
                    with open(path) as f:
                        adjs.append(json.load(f))
                except FileNotFoundError:
                    adjs = []
                    break
            if not adjs:
                continue

            # find keys common to all
            common = set(adjs[0].keys())
            for d in adjs[1:]:
                common &= set(d.keys())
            if not common:
                continue

            # pick Eb from first period
            Eb = beam_energies[periods[0]]

            for key_str in sorted(common, key=key_to_tuple):
                # gather values
                As = []
                sigAs = []
                for d in adjs:
                    v = d[key_str]
                    if not v.get("valid", True):
                        As = []
                        break
                    As.append(v["bsa"])
                    sigAs.append(v["bsa_err"])
                if not As:
                    continue

                # filter zero-error
                good = [(A, s) for A, s in zip(As, sigAs) if s > 0]
                if not good:
                    continue
                As, sigAs = zip(*good)

                # weighted average
                weights = [1.0/(s*s) for s in sigAs]
                Wsum = sum(weights)
                A_avg = sum(w*a for w,a in zip(weights, As)) / Wsum
                sigA = math.sqrt(1.0 / Wsum)

                # kinematics
                m = bin_means.get(key_str)
                if m is None:
                    continue
                phi_deg = math.degrees(m["phi_avg"])
                q2 = m["Q2_avg"]
                xb = m["xB_avg"]
                t_val = -m["t_avg"]

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