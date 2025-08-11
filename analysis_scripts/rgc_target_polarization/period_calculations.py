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
by_period = {"RGC_Su22": [], "RGC_Fa22": [], "RGC_Sp23": []}  # list of (run, pt_avg, sig_stat, sys_abs)

with open(PT_FILE) as f:
    for raw in f:
        line = raw.strip()
        if (not line) or line.startswith("#") or line.lower().startswith("run"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        try:
            run      = int(parts[0])
            pt_avg   = float(parts[7])
            sig_stat = float(parts[8])
            sys_abs  = float(parts[9])
        except ValueError:
            continue
        period = run_to_period.get(run)
        if period is not None:
            by_period[period].append((run, pt_avg, sig_stat, sys_abs))

# ---- 3) Helpers for weighted mean/stat and relative->absolute sys ----
def weighted_mean_and_stat_with_rel_sys(entries):
    """
    entries: list of (run, y, sigma_stat, sigma_sys_abs)

    Returns:
      n_used         : number of runs used for the weighted mean/stat
      mean_w         : weighted mean of y with weights w = 1/sigma_stat^2
      stat_unc       : sqrt(1 / sum w)
      sys_abs_from_rel: absolute sys built as (mean of per-run relative sys) * |mean_w|

    Notes:
      - Relative sys for a run is defined as sigma_sys_abs / |y|.
      - We average those relative sys values over runs where |y| > 0 and both y, sigma_sys_abs are finite.
      - Final period sys is that mean-relative multiplied by |mean_w|.
    """
    if not entries:
        return 0, float("nan"), float("nan"), float("nan")

    # For stat-weighted mean:
    ys, sigs_stat = [], []
    # For relative sys averaging:
    rel_sys_vals = []

    for _, y, s_stat, s_sys_abs in entries:
        # collect for stat-weighted mean if stat is finite and > 0
        if np.isfinite(s_stat) and s_stat > 0 and np.isfinite(y):
            ys.append(float(y))
            sigs_stat.append(float(s_stat))
        # collect relative sys if finite and |y| > 0
        if np.isfinite(y) and np.isfinite(s_sys_abs) and abs(y) > 0.0:
            rel_sys_vals.append(float(s_sys_abs) / abs(float(y)))

    n_used = len(ys)
    if n_used == 0:
        return 0, float("nan"), float("nan"), float("nan")

    ys        = np.array(ys, dtype=float)
    sigs_stat = np.array(sigs_stat, dtype=float)

    w = 1.0 / (sigs_stat**2)
    sumw = np.sum(w)
    if not np.isfinite(sumw) or sumw <= 0:
        return 0, float("nan"), float("nan"), float("nan")

    mean_w   = float(np.sum(w * ys) / sumw)
    stat_unc = float(np.sqrt(1.0 / sumw))

    # Mean of relative sys across qualifying runs
    if len(rel_sys_vals) == 0:
        sys_abs_from_rel = float("nan")
    else:
        rel_sys_mean = float(np.mean(rel_sys_vals))
        sys_abs_from_rel = rel_sys_mean * abs(mean_w)

    return n_used, mean_w, stat_unc, sys_abs_from_rel

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

    npos, mpos, spos, upos_abs = weighted_mean_and_stat_with_rel_sys(pos_entries)
    nneg, mneg, sneg, uneg_abs = weighted_mean_and_stat_with_rel_sys(neg_entries)

    def fmt(x):
        return f"{x:.4f}" if np.isfinite(x) else "nan"

    rows.append(
        f"{period}\t"
        f"{npos}\t{fmt(mpos)}\t{fmt(spos)}\t{fmt(upos_abs)}\t"
        f"{nneg}\t{fmt(mneg)}\t{fmt(sneg)}\t{fmt(uneg_abs)}"
    )

# ---- 5) Save and print ----
out_path = "output/period_pt_sign_means.txt"
with open(out_path, "w") as out:
    for ln in rows:
        out.write(ln + "\n")

print("\n".join(rows))
print(f"[INFO] Wrote {out_path}")