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
    final_results_dir  : path to 'final_results' dir
    output_file        : filename for the combined .txt (e.g. "bsa_all_periods.txt")

    Writes a single file with header and then one line per period+bin.
    Rounds all numeric columns to two decimal places.
    """
    # load global bin kinematics
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    # ensure output directory exists
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # desired chronological order
    desired_order = [
        "DVCS_Sp18_out",
        "DVCS_Sp18_inb",
        "DVCS_Fa18_inb",
        "DVCS_Fa18_out",
        "DVCS_Sp19_inb",
    ]
    period_list = [p for p in desired_order if p in periods]

    # beam energies lookup
    beam_energies = {
        "DVCS_Fa18_inb":  10.604,
        "DVCS_Fa18_out":  10.604,
        "DVCS_Sp19_inb":  10.200,
        "DVCS_Sp18_inb":  10.594,
        "DVCS_Sp18_out":  10.594,
    }

    # helper to sort bin keys
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").replace(" ", "").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for period in period_list:
            Eb = beam_energies.get(period)
            if Eb is None:
                continue

            path = os.path.join(final_results_dir, f"adjusted_bsa_{period}.json")
            try:
                with open(path) as f:
                    adj = json.load(f)
            except FileNotFoundError:
                continue

            for key_str in sorted(adj.keys(), key=key_to_tuple):
                vals = adj[key_str]
                if not vals.get("valid", True):
                    continue
                mean = bin_means.get(key_str)
                if mean is None:
                    continue

                # compute values
                phi_deg = math.degrees(mean["phi_avg"])
                q2      = mean["Q2_avg"]
                xb      = mean["xB_avg"]
                t_val   = -mean["t_avg"]
                A       = vals["bsa"]
                sigmaA  = vals["bsa_err"]

                # write rounded to two decimals
                out.write(
                    f"{phi_deg:.2f} "
                    f"{q2:.2f} "
                    f"{xb:.2f} "
                    f"{t_val:.2f} "
                    f"{Eb:.2f} "
                    f"{A:.2f} "
                    f"{sigmaA:.2f}\n"
                )
    print(f"[export] wrote combined file → {output_path}")



def export_bsa_grouped_to_text(bin_means_json,
                               final_results_dir,
                               output_file):
    """
    Combines:
      - SP18_out & SP18_inb → first block
      - FA18_inb & FA18_out → second block
      - SP19_inb           → third block

    Writes one file with no extra headers, rounded to two decimals.
    Uses weighted averages of A (weights=1/σA²), skipping σA<=0 periods.
    """
    with open(bin_means_json) as f:
        bin_means = json.load(f)
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    group_defs = [
        ["DVCS_Sp18_out", "DVCS_Sp18_inb"],
        ["DVCS_Fa18_inb", "DVCS_Fa18_out"],
        ["DVCS_Sp19_inb"],
    ]
    beam_energies = {
        "DVCS_Fa18_inb": 10.604,
        "DVCS_Fa18_out": 10.604,
        "DVCS_Sp19_inb": 10.200,
        "DVCS_Sp18_inb": 10.594,
        "DVCS_Sp18_out": 10.594,
    }
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").replace(" ", "").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")

        for periods in group_defs:
            adjs = []
            for p in periods:
                path = os.path.join(final_results_dir, f"adjusted_bsa_{p}.json")
                try:
                    with open(path) as f:
                        adjs.append(json.load(f))
                except FileNotFoundError:
                    adjs = []
                    break
            if not adjs:
                continue

            common = set(adjs[0].keys())
            for d in adjs[1:]:
                common &= set(d.keys())
            if not common:
                continue

            Eb = beam_energies[periods[0]]
            for key_str in sorted(common, key=key_to_tuple):
                As, sigAs = [], []
                for d in adjs:
                    v = d[key_str]
                    if not v.get("valid", True):
                        As = []
                        break
                    As.append(v["bsa"])
                    sigAs.append(v["bsa_err"])
                if not As:
                    continue
                # filter and average
                good = [(A,s) for A,s in zip(As,sigAs) if s>0]
                if not good:
                    continue
                As, sigAs = zip(*good)
                weights = [1.0/(s*s) for s in sigAs]
                W = sum(weights)
                A_avg = sum(w*a for w,a in zip(weights,As)) / W
                sigA = math.sqrt(1.0 / W)

                m = bin_means.get(key_str)
                if m is None:
                    continue
                phi_deg = math.degrees(m["phi_avg"])
                q2      = m["Q2_avg"]
                xb      = m["xB_avg"]
                t_val   = -m["t_avg"]

                out.write(
                    f"{phi_deg:.2f} "
                    f"{q2:.2f} "
                    f"{xb:.2f} "
                    f"{t_val:.2f} "
                    f"{Eb:.2f} "
                    f"{A_avg:.2f} "
                    f"{sigA:.2f}\n"
                )
    print(f"[grouped export] wrote → {output_path}")


def export_bsa_spfa_combined_to_text(bin_means_json,
                                      final_results_dir,
                                      output_file):
    """
    Combines SP18 (out+inb) and FA18 (inb+out) into one group,
    treating Eb=10.60 GeV, then writes SP19_inb after.

    Weighted averages for A, sigma; rounds to two decimals.
    """
    with open(bin_means_json) as f:
        bin_means = json.load(f)
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # two groups: [Sp18+Fa18], [Sp19]
    groups = [
        ["DVCS_Sp18_out", "DVCS_Sp18_inb", "DVCS_Fa18_inb", "DVCS_Fa18_out"],
        ["DVCS_Sp19_inb"],
    ]
    # Eb for first combined group
    Eb_group1 = 10.60
    # lookup for Sp19
    beam_energies = {"DVCS_Sp19_inb": 10.200}
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").replace(" ", "").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for idx, periods in enumerate(groups):
            # load JSONs
            adjs = []
            for p in periods:
                path = os.path.join(final_results_dir, f"adjusted_bsa_{p}.json")
                try:
                    with open(path) as f:
                        adjs.append(json.load(f))
                except FileNotFoundError:
                    adjs = []
                    break
            if not adjs:
                continue

            common = set(adjs[0].keys())
            for d in adjs[1:]:
                common &= set(d.keys())
            if not common:
                continue

            Eb = Eb_group1 if idx == 0 else beam_energies.get(periods[0], 0)

            for key_str in sorted(common, key=key_to_tuple):
                As, sigAs = [], []
                for d in adjs:
                    v = d[key_str]
                    if not v.get("valid", True):
                        As = []
                        break
                    As.append(v["bsa"])
                    sigAs.append(v["bsa_err"])
                if not As:
                    continue
                good = [(A,s) for A,s in zip(As,sigAs) if s>0]
                if not good:
                    continue
                As, sigAs = zip(*good)
                weights = [1.0/(s*s) for s in sigAs]
                W = sum(weights)
                A_avg = sum(w*a for w,a in zip(weights,As)) / W
                sigA = math.sqrt(1.0 / W)

                m = bin_means.get(key_str)
                if m is None:
                    continue
                phi_deg = math.degrees(m["phi_avg"])
                q2      = m["Q2_avg"]
                xb      = m["xB_avg"]
                t_val   = -m["t_avg"]

                out.write(
                    f"{phi_deg:.2f} "
                    f"{q2:.2f} "
                    f"{xb:.2f} "
                    f"{t_val:.2f} "
                    f"{Eb:.2f} "
                    f"{A_avg:.2f} "
                    f"{sigA:.2f}\n"
                )
    print(f"[merged export] wrote → {output_path}")