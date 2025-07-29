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
    # load bin means
    with open(bin_means_json) as f:
        bin_means = json.load(f)

    os.makedirs(final_results_dir, exist_ok=True)
    output_path = os.path.join(final_results_dir, output_file)

    # periods in desired order
    desired = [
        "DVCS_Sp18_out", "DVCS_Sp18_inb",
        "DVCS_Fa18_inb", "DVCS_Fa18_out",
        "DVCS_Sp19_inb"
    ]
    period_list = [p for p in desired if p in periods]

    # beam energies
    beam_energies = {
        "DVCS_Fa18_inb":  10.604,
        "DVCS_Fa18_out":  10.604,
        "DVCS_Sp19_inb":  10.200,
        "DVCS_Sp18_inb":  10.594,
        "DVCS_Sp18_out":  10.594,
    }

    def key_to_tuple(k):
        return tuple(map(int, k.strip("()").split(",")))

    with open(output_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for period in period_list:
            Eb = beam_energies.get(period)
            if Eb is None:
                continue
            path = os.path.join(final_results_dir, f"adjusted_bsa_{period}.json")
            try:
                with open(path) as f:
                    data = json.load(f)
            except FileNotFoundError:
                continue

            for key in sorted(data.keys(), key=key_to_tuple):
                vals = data[key]
                if not vals.get("valid", True):
                    continue
                A = vals["bsa"]
                sigmaA = vals["bsa_err"]
                # skip out-of-range
                if abs(A) > 1 or abs(sigmaA) > 1:
                    continue
                m = bin_means.get(key)
                if m is None:
                    continue
                phi_deg = math.degrees(m["phi_avg"])
                q2 = m["Q2_avg"]
                xb = m["xB_avg"]
                t_val = -m["t_avg"]
                out.write(
                    f"{phi_deg:.2f} {q2:.2f} {xb:.2f} {t_val:.2f} {Eb:.2f} {A:.2f} {sigmaA:.2f}\n"
                )
    print(f"[export] wrote → {output_path}")


def export_bsa_grouped_to_text(bin_means_json,
                               final_results_dir,
                               output_file):
    """
    Groups: [Sp18_out,Sp18_inb], [Fa18_inb,Fa18_out], [Sp19_inb].
    Union of keys; for each key average over available periods,
    skipping zero or invalid entries, two-decimal rounding,
    skip |A|>1 or |σA|>1.
    """
    with open(bin_means_json) as f:
        bin_means = json.load(f)
    os.makedirs(final_results_dir, exist_ok=True)
    out_path = os.path.join(final_results_dir, output_file)

    groups = [
        ["DVCS_Sp18_out","DVCS_Sp18_inb"],
        ["DVCS_Fa18_inb","DVCS_Fa18_out"],
        ["DVCS_Sp19_inb"]
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

    with open(out_path, "w") as out:
        out.write("# phi(deg) q2(GeV2) xb t(GeV2) Eb(GeV) A sigA\n")
        for periods in groups:
            # load JSONs
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
            # union keys
            keys = set().union(*[d.keys() for d in adjs])
            Eb = beam_energies.get(periods[0], 0)
            for key in sorted(keys, key=key_to_tuple):
                As, sigAs = [], []
                for d in adjs:
                    if key in d and d[key].get("valid",True):
                        a = d[key]["bsa"]
                        sa= d[key]["bsa_err"]
                        if sa>0 and abs(a)<=1 and abs(sa)<=1:
                            As.append(a); sigAs.append(sa)
                if not As:
                    continue
                # weighted average
                ws = [1/(s*s) for s in sigAs]
                W = sum(ws)
                A_avg = sum(w*a for w,a in zip(ws,As))/W
                sigmaA = math.sqrt(1.0/W)
                m = bin_means.get(key)
                if not m:
                    continue
                phi_deg = math.degrees(m["phi_avg"])
                q2 = m["Q2_avg"]; xb = m["xB_avg"]
                t_val = -m["t_avg"]
                out.write(
                    f"{phi_deg:.2f} {q2:.2f} {xb:.2f} {t_val:.2f} {Eb:.2f} 