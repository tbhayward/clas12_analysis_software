#!/usr/bin/env python3
import os
import numpy as np

# ---- paths (adjust if needed) ----
RUNINFO = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"
PT_FILE = "output/Pt_by_run.txt"

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

# ---- 1) Build runnum -> period map using NH3 sections in RUNINFO ----
PERIOD_KEYS = {
    "rgc su22 nh3": "RGC_Su22",
    "rgc fa22 nh3": "RGC_Fa22",
    "rgc sp23 nh3": "RGC_Sp23",
}

run_to_period = {}
current_period = None

with open(RUNINFO) as f:
    for raw in f:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            key = line.lower().strip("# ").strip()
            current_period = PERIOD_KEYS.get(key, None)
            continue
        if current_period is None:
            continue  # only keep NH3 sections we recognize
        parts = [p.strip() for p in line.split(",")]
        if not parts or not parts[0].isdigit():
            continue
        runnum = int(parts[0])
        run_to_period[runnum] = current_period

# ---- 2) Read DIS results (Pt_by_run.txt), grab Pt_avg, avg_sig (stat), avg_sys (sys) ----
# Format: Run Pt_GRV sigma_GRV sys_GRV Pt_ABD sigma_ABD sys_ABD Pt_avg avg_sig avg_sys
by_period = {"RGC_Su22": [], "RGC_Fa22": [], "RGC_Sp23": []}  # list of (run, pt_avg, sig_stat, sys)

with open(PT_FILE) as f:
    for raw in f:
        line = raw.strip()
        if (not line) or line.startswith("#") or line.lower().startswith("run"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        try:
            run     = int(parts[0])
            pt_avg  = float(parts[7])
            sig_stat= float(parts[8])
            sys_avg = float(parts[9])
        except ValueError:
            continue
        period = run_to_period.get(run)
        if period is not None:
            by_period[period].append((run, pt_avg, sig_stat, sys_avg))

# ---- 3) Helpers for weighted mean/stat and sys mean ----
def weighted_mean_and_stat(entries):
    """
    entries: list of (run, y, sigma_stat, sigma_sys)
    Returns: n_used, mean_w, stat_unc, sys_mean
    - mean_w is the weighted mean of y with w=1/sigma_stat^2
    - stat_unc = sqrt(1/sum w)
    - sys_mean = arithmetic mean of sigma_sys over the used entries
    """
    if not entries:
        return 0, float("nan"), float("nan"), float("nan")

    # Keep only entries with positive, finite statistical uncertainty
    ys, sigs_stat, sigs_sys = [], [], []
    for _, y, s_stat, s_sys in entries:
        if np.isfinite(s_stat) and s_stat > 0:
            ys.append(y)
            sigs_stat.append(s_stat)
            sigs_sys.append(s_sys)

    n_used = len(ys)
    if n_used == 0:
        return 0, float("nan"), float("nan"), float("nan")

    ys        = np.array(ys, dtype=float)
    sigs_stat = np.array(sigs_stat, dtype=float)
    sigs_sys  = np.array(sigs_sys, dtype=float)

    w = 1.0 / (sigs_stat**2)
    sumw = np.sum(w)
    if not np.isfinite(sumw) or sumw <= 0:
        return 0, float("nan"), float("nan"), float("nan")

    mean_w   = float(np.sum(w * ys) / sumw)
    stat_unc = float(np.sqrt(1.0 / sumw))
    sys_mean = float(np.mean(sigs_sys)) if sigs_sys.size > 0 else float("nan")
    return n_used, mean_w, stat_unc, sys_mean

# ---- 4) For each period, split into pos/neg and compute stats ----
rows = []
header = ("Period\t"
          "N_pos\tMean_Pos\tStat_Pos\tSys_Pos\t"
          "N_neg\tMean_Neg\tStat_Neg\tSys_Neg")
rows.append(header)

for period in ["RGC_Su22", "RGC_Fa22", "RGC_Sp23"]:
    entries = by_period.get(period, [])

    pos_entries = [(r, y, s, u) for (r, y, s, u) in entries if y > 0]
    neg_entries = [(r, y, s, u) for (r, y, s, u) in entries if y < 0]

    npos, mpos, spos, upos = weighted_mean_and_stat(pos_entries)
    nneg, mneg, sneg, uneg = weighted_mean_and_stat(neg_entries)

    rows.append(
        f"{period}\t"
        f"{npos}\t{(mpos if np.isfinite(mpos) else float('nan')):.4f}\t{(spos if np.isfinite(spos) else float('nan')):.4f}\t{(upos if np.isfinite(upos) else float('nan')):.4f}\t"
        f"{nneg}\t{(mneg if np.isfinite(mneg) else float('nan')):.4f}\t{(sneg if np.isfinite(sneg) else float('nan')):.4f}\t{(uneg if np.isfinite(uneg) else float('nan')):.4f}"
    )

# ---- 5) Save and print ----
out_path = "output/period_pt_sign_means.txt"
with open(out_path, "w") as out:
    for ln in rows:
        out.write(ln + "\n")

print("\n".join(rows))
print(f"[INFO] Wrote {out_path}")