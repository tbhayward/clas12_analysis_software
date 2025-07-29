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
    output_file        : filename for the combined .txt

    Writes one file with header plus one line per period+bin,
    two-decimal rounding, skips |A|>1 or |σA|>1.
    """
    # 1) load global bin kinematics
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    # 2) ensure output directory exists
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # 3) define desired period order
    desired_order = [
        "DVCS_Sp18_out",
        "DVCS_Sp18_inb",
        "DVCS_Fa18_inb",
        "DVCS_Fa18_out",
        "DVCS_Sp19_inb",
    ]
    period_list = [p for p in desired_order if p in periods]

    # 4) beam energies lookup
    beam_energies = {
        "DVCS_Fa18_inb":  10.604,
        "DVCS_Fa18_out":  10.604,
        "DVCS_Sp19_inb":  10.200,
        "DVCS_Sp18_inb":  10.594,
        "DVCS_Sp18_out":  10.594,
    }

    # helper to sort keys
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").split(",")))

    # 5) write header and data
    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for period in period_list:
            Eb = beam_energies.get(period)
            if Eb is None:
                continue
            json_path = os.path.join(final_results_dir, f"adjusted_bsa_{period}.json")
            try:
                with open(json_path) as f:
                    data = json.load(f)
            except FileNotFoundError:
                continue

            for key in sorted(data.keys(), key=key_to_tuple):
                record = data[key]
                if not record.get("valid", True):
                    continue
                A = record.get("bsa")
                sigmaA = record.get("bsa_err")
                # skip if outside [-1,1]
                if A is None or sigmaA is None or abs(A) > 1 or abs(sigmaA) > 1:
                    continue

                mean = bin_means.get(key)
                if mean is None:
                    continue
                phi_deg = math.degrees(mean.get("phi_avg", 0))
                q2      = mean.get("Q2_avg", 0)
                xb      = mean.get("xB_avg", 0)
                t_val   = -mean.get("t_avg", 0)

                out.write(
                    f"{phi_deg:.2f} {q2:.2f} {xb:.2f} {t_val:.2f} {Eb:.2f} "
                    f"{A:.2f} {sigmaA:.2f}\n"
                )
    print(f"[export_bsa_to_text] wrote → {output_path}")



def export_bsa_grouped_to_text(bin_means_json,
                               final_results_dir,
                               output_file):
    """
    Groups:
      - SP18_out & SP18_inb
      - FA18_inb & FA18_out
      - SP19_inb

    Uses union of all bins within each group; if a bin
    appears in only one period, it is still output.
    Weighted average of A over available periods, weights=1/sigmaA^2.
    Skips entries with sigmaA<=0 or |A|>1 or |sigmaA|>1.
    Rounds to two decimals.
    """
    with open(bin_means_json) as f:
        bin_means = json.load(f)
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    groups = [
        ["DVCS_Sp18_out", "DVCS_Sp18_inb"],
        ["DVCS_Fa18_inb", "DVCS_Fa18_out"],
        ["DVCS_Sp19_inb"],
    ]
    beam_energies = {
        "DVCS_Sp18_out": 10.594,
        "DVCS_Sp18_inb": 10.594,
        "DVCS_Fa18_inb": 10.604,
        "DVCS_Fa18_out": 10.604,
        "DVCS_Sp19_inb": 10.200,
    }
    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for periods in groups:
            # load all JSONs in this group
            adjs = []
            for p in periods:
                pth = os.path.join(final_results_dir, f"adjusted_bsa_{p}.json")
                try:
                    with open(pth) as f:
                        adjs.append(json.load(f))
                except FileNotFoundError:
                    continue
            if not adjs:
                continue

            # union of all keys
            keys = set().union(*(d.keys() for d in adjs))
            Eb = beam_energies.get(periods[0], 0)

            for key in sorted(keys, key=key_to_tuple):
                # collect available A, sigmaA
                As = []
                sigAs = []
                for d in adjs:
                    rec = d.get(key)
                    if not rec or not rec.get("valid", True):
                        continue
                    a = rec.get("bsa")
                    sa = rec.get("bsa_err")
                    # require sa>0 and within bounds
                    if sa and sa>0 and abs(a)<=1 and abs(sa)<=1:
                        As.append(a)
                        sigAs.append(sa)
                if not As:
                    continue

                # weighted average
                weights = [1.0/(s*s) for s in sigAs]
                Wsum = sum(weights)
                A_avg = sum(w*a for w,a in zip(weights, As)) / Wsum
                sigmaA = math.sqrt(1.0 / Wsum)

                mean = bin_means.get(key)
                if not mean:
                    continue
                phi_deg = math.degrees(mean.get("phi_avg", 0))
                q2      = mean.get("Q2_avg", 0)
                xb      = mean.get("xB_avg", 0)
                t_val   = -mean.get("t_avg", 0)

                out.write(
                    f"{phi_deg:.2f} {q2:.2f} {xb:.2f} {t_val:.2f} {Eb:.2f} "
                    f"{A_avg:.2f} {sigmaA:.2f}\n"
                )
    print(f"[export_bsa_grouped_to_text] wrote → {output_path}")



def export_bsa_spfa_combined_to_text(bin_means_json,
                                      final_results_dir,
                                      output_file):
    """
    Combines SP18 & FA18 into a single group (Eb=10.60 GeV),
    then SP19_inb separately.

    Union of keys across all periods; averages available
    values with weights=1/sigmaA^2; skips invalid.
    Rounds to two decimals, skip |A|>1 or |σA|>1.
    """
    with open(bin_means_json) as f:
        bin_means = json.load(f)
    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    groups = [
        ["DVCS_Sp18_out", "DVCS_Sp18_inb",
         "DVCS_Fa18_inb", "DVCS_Fa18_out"],
        ["DVCS_Sp19_inb"],
    ]
    Eb_values = [10.60, None]
    beam_energies = {"DVCS_Sp19_inb": 10.200}

    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for idx, periods in enumerate(groups):
            adjs = []
            for p in periods:
                pth = os.path.join(final_results_dir, f"adjusted_bsa_{p}.json")
                try:
                    with open(pth) as f:
                        adjs.append(json.load(f))
                except FileNotFoundError:
                    continue
            if not adjs:
                continue

            keys = set().union(*(d.keys() for d in adjs))
            Eb = Eb_values[idx] if Eb_values[idx] else beam_energies.get(periods[0],0)

            for key in sorted(keys, key=key_to_tuple):
                As = []
                sigAs = []
                for d in adjs:
                    rec = d.get(key)
                    if not rec or not rec.get("valid", True):
                        continue
                    a = rec.get("bsa")
                    sa = rec.get("bsa_err")
                    if sa and sa>0 and abs(a)<=1 and abs(sa)<=1:
                        As.append(a)
                        sigAs.append(sa)
                if not As:
                    continue

                weights = [1.0/(s*s) for s in sigAs]
                Wsum = sum(weights)
                A_avg = sum(w*a for w,a in zip(weights, As)) / Wsum
                sigmaA = math.sqrt(1.0 / Wsum)

                mean = bin_means.get(key)
                if not mean:
                    continue
                phi_deg = math.degrees(mean.get("phi_avg", 0))
                q2      = mean.get("Q2_avg", 0)
                xb      = mean.get("xB_avg", 0)
                t_val   = -mean.get("t_avg", 0)

                out.write(
                    f"{phi_deg:.2f} {q2:.2f} {xb:.2f} {t_val:.2f} {Eb:.2f} "
                    f"{A_avg:.2f} {sigmaA:.2f}\n"
                )
    print(f"[export_bsa_spfa_combined_to_text] wrote → {output_path}")