#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from matplotlib.patches import Patch

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

# ===============================
# Utilities
# ===============================

def plot_pos_neg(ax, x, y, yerr, color, label, xoffset=0.0, alpha=1.0, zorder=None):
    """Plot |y| with filled markers for y>=0 and empty markers for y<0 on a given axes."""
    x = np.asarray(x) + xoffset
    y = np.asarray(y)
    yerr = None if yerr is None else np.asarray(yerr)

    pos = y >= 0
    neg = ~pos
    y_abs = np.abs(y)

    if pos.any():
        ax.errorbar(x[pos], y_abs[pos], yerr=None if yerr is None else yerr[pos],
                    fmt='o', color=color, label=label, alpha=alpha, zorder=zorder)
    if neg.any():
        ax.errorbar(x[neg], y_abs[neg], yerr=None if yerr is None else yerr[neg],
                    fmt='o', color=color, mfc='none', label=None, alpha=alpha, zorder=zorder)

def segment_bands_by_value(run_x, mean, err):
    """
    Given arrays (run_x, mean, err) for NH3 elastic values,
    return list of segments where (mean, err) is constant.
    Each segment is (x_start, x_end, y_low, y_high), with y_low/high using |mean| ± err.
    Segmentation is only by change in (mean, err), not by run contiguity.
    """
    bands = []
    if len(run_x) == 0:
        return bands
    start = 0
    current_val = mean[0]
    current_err = err[0]
    for i in range(1, len(run_x)):
        if (mean[i] != current_val) or (err[i] != current_err):
            x0 = run_x[start]
            x1 = run_x[i-1]
            m  = abs(current_val)
            e  = current_err
            bands.append((x0, x1, max(0.0, m - e), m + e))
            start = i
            current_val = mean[i]
            current_err = err[i]
    # last segment
    x0 = run_x[start]
    x1 = run_x[-1]
    m  = abs(current_val)
    e  = current_err
    bands.append((x0, x1, max(0.0, m - e), m + e))
    return bands

def draw_bands(ax, bands, color='green', alpha=0.2, label=None):
    """
    Draw a list of (x0, x1, y0, y1) bands as filled rectangles.
    Adds a legend entry once if 'label' provided.
    """
    first = True
    for (x0, x1, y0, y1) in bands:
        ax.fill_between([x0, x1], [y0, y0], [y1, y1],
                        color=color, alpha=alpha, label=label if first else None)
        first = False

# ===============================
# 1) Model curves plot
# ===============================

def ALL_GRV(x):
    return 0.00823729 + 1.62853*x - 1.38493*x**2 + 1.07047*x**3 - 0.747653*x**4

def ALL_ABD(x):
    return 0.0558035 + 1.23137*x - 1.05596*x**2 + 1.95783*x**3 - 1.22263*x**4

x_vals = np.linspace(0.0, 0.7, 500)
y_grv  = ALL_GRV(x_vals)
y_abd  = ALL_ABD(x_vals)

plt.figure(figsize=(8,6))
plt.plot(x_vals, y_grv, color='blue', label="GRSV [hep-ph] 0011215v1")
plt.plot(x_vals, y_abd, color='red',  label="ABDY [hep-ph] 0705.1553")
plt.axhline(0, color='grey', linestyle='--', linewidth=1)
plt.xlabel(r"$x_{B}$", fontsize=14)
plt.ylabel(r"$F_{LL}/F_{UU}$", fontsize=14)
plt.ylim(-0.2, 0.80)
plt.xlim(0.0, 0.7)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/models_plot.pdf")
plt.close()

print("[INFO] Saved output/models_plot.pdf")

# ===============================
# 2) Load run-by-run extraction data from C++ output
# New format:
# Run Pt_GRV sigma_GRV sys_GRV Pt_ABD sigma_ABD sys_ABD Pt_avg avg_sig avg_sys
# ===============================
extractions_file = "output/Pt_by_run.txt"
runnum, Pt_grv, s_grv, sys_grv, Pt_abd, s_abd, sys_abd, Pt_avg, avg_sig, avg_sys = ([] for _ in range(10))

with open(extractions_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("Run"):
            continue
        parts = line.split()
        if len(parts) < 10:
            continue
        try:
            runnum.append(int(parts[0]))
            Pt_grv.append(float(parts[1]))
            s_grv.append(float(parts[2]))
            sys_grv.append(float(parts[3]))
            Pt_abd.append(float(parts[4]))
            s_abd.append(float(parts[5]))
            sys_abd.append(float(parts[6]))
            Pt_avg.append(float(parts[7]))
            avg_sig.append(float(parts[8]))
            avg_sys.append(float(parts[9]))
        except ValueError:
            continue

runnum  = np.array(runnum)
Pt_grv  = np.array(Pt_grv);   s_grv  = np.array(s_grv);   sys_grv = np.array(sys_grv)
Pt_abd  = np.array(Pt_abd);   s_abd  = np.array(s_abd);   sys_abd = np.array(sys_abd)
Pt_avg  = np.array(Pt_avg);   avg_sig= np.array(avg_sig); avg_sys = np.array(avg_sys)

# sort by runnum (for consistent left/right mapping)
order = np.argsort(runnum)
runnum   = runnum[order]
Pt_grv   = Pt_grv[order];  s_grv   = s_grv[order];  sys_grv = sys_grv[order]
Pt_abd   = Pt_abd[order];  s_abd   = s_abd[order];  sys_abd = sys_abd[order]
Pt_avg   = Pt_avg[order];  avg_sig = avg_sig[order]; avg_sys = avg_sys[order]

xmin = runnum.min() - 10
xmax = runnum.max() + 10
run_index = np.arange(1, len(runnum)+1)
run_index_scaled = run_index * 5  # scaled x-axis for visibility

# ===============================
# 3) 1x2: Pt per run (GRV + ABD), with stat and stat+sys (per-model)
# ===============================
tot_grv = np.sqrt(s_grv**2 + sys_grv**2)
tot_abd = np.sqrt(s_abd**2 + sys_abd**2)

fig, axes = plt.subplots(1, 2, figsize=(16,6), sharey=True)

# Left: vs runnum (offset ±0.005)
ax = axes[0]
plot_pos_neg(ax, runnum, Pt_grv, s_grv,  color='blue', label="GRSV [stat]",      xoffset=-0.005)
plot_pos_neg(ax, runnum, Pt_grv, tot_grv, color='blue', label="GRSV [stat+sys]", xoffset=-0.005, alpha=0.4)
plot_pos_neg(ax, runnum, Pt_abd, s_abd,  color='red',  label="ABDY [stat]",      xoffset=+0.005)
plot_pos_neg(ax, runnum, Pt_abd, tot_abd, color='red', label="ABDY [stat+sys]",  xoffset=+0.005, alpha=0.4)
ax.set_xlabel("runnum", fontsize=14)
ax.set_ylabel(r"$P_{t}$", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(xmin, xmax)
ax.legend(fontsize=10)

# Right: vs run index × 5 (offsets multiplied by 5 -> ±1.25)
ax = axes[1]
plot_pos_neg(ax, run_index_scaled, Pt_grv, s_grv,  color='blue', label="GRSV [stat]",      xoffset=-1.25)
plot_pos_neg(ax, run_index_scaled, Pt_grv, tot_grv, color='blue', label="GRSV [stat+sys]", xoffset=-1.25, alpha=0.4)
plot_pos_neg(ax, run_index_scaled, Pt_abd, s_abd,  color='red',  label="ABDY [stat]",      xoffset=+1.25)
plot_pos_neg(ax, run_index_scaled, Pt_abd, tot_abd, color='red', label="ABDY [stat+sys]",  xoffset=+1.25, alpha=0.4)
ax.set_xlabel("run index × 5", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(0, run_index_scaled.max() + 5)

plt.tight_layout()
plt.savefig("output/model_extractions.pdf")
plt.close()
print("[INFO] Saved output/model_extractions.pdf")

# ===============================
# 4) 1x2: Average Pt with stat and stat+sys (avg_sys from file)
# ===============================
stat_plus_sys = np.sqrt(avg_sig**2 + avg_sys**2)

fig, axes = plt.subplots(1, 2, figsize=(16,6), sharey=True)

# Left: vs runnum
ax = axes[0]
plot_pos_neg(ax, runnum, Pt_avg, avg_sig,       color='black', label="Avg $P_{t}$ (stat)")
plot_pos_neg(ax, runnum, Pt_avg, stat_plus_sys, color='black', label="Avg $P_{t}$ (stat+sys)", alpha=0.4)
ax.set_xlabel("runnum", fontsize=14)
ax.set_ylabel(r"$P_{t}$", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(xmin, xmax)
ax.legend(fontsize=10)

# Right: vs run index × 5 (offsets multiplied by 5 where used; here no offset needed)
ax = axes[1]
plot_pos_neg(ax, run_index_scaled, Pt_avg, avg_sig,       color='black', label="Avg $P_{t}$ (stat)")
plot_pos_neg(ax, run_index_scaled, Pt_avg, stat_plus_sys, color='black', label="Avg $P_{t}$ (stat+sys)", alpha=0.4)
ax.set_xlabel("run index × 5", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(0, run_index_scaled.max() + 5)

plt.tight_layout()
plt.savefig("output/avg_pt_plot.pdf")
plt.close()
print("[INFO] Saved output/avg_pt_plot.pdf")

# ===============================
# 5) 1x2: Method comparison (NH3-only from CSV) with Elastic bands
# ===============================
run_csv_file = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

run_csv_nums = []
pol_csv      = []
pol_csv_err  = []

current_section_is_nh3 = False
with open(run_csv_file, newline='') as csvfile:
    for raw in csvfile:
        line = raw.strip()
        if not line:
            continue
        if line.startswith("#"):
            header = line.lower()
            current_section_is_nh3 = ("nh3" in header)
            continue
        if not current_section_is_nh3:
            continue
        row = line.split(",")
        if len(row) < 3:
            continue
        try:
            run_csv_nums.append(int(row[0]))
            pol_csv.append(float(row[-2]))      # mean polarization (Elastic method)
            pol_csv_err.append(float(row[-1]))  # its stat uncertainty
        except ValueError:
            continue

run_csv_nums = np.array(run_csv_nums)
pol_csv      = np.array(pol_csv)
pol_csv_err  = np.array(pol_csv_err)

# Sort Elastic arrays by run
if run_csv_nums.size > 0:
    order_np = np.argsort(run_csv_nums)
    run_csv_nums = run_csv_nums[order_np]
    pol_csv      = pol_csv[order_np]
    pol_csv_err  = pol_csv_err[order_np]
run_index_np = np.arange(1, len(run_csv_nums)+1) if run_csv_nums.size > 0 else np.array([])
run_index_np_scaled = run_index_np * 5

# Build band segments where (pol, err) is constant
elastic_bands_runnum = segment_bands_by_value(run_csv_nums, pol_csv, pol_csv_err)
elastic_bands_index  = segment_bands_by_value(run_index_np_scaled, pol_csv, pol_csv_err)

fig, axes = plt.subplots(1, 2, figsize=(16,6), sharey=True)

# Left: vs runnum (DIS points; Elastic band)
ax = axes[0]
plot_pos_neg(ax, runnum, Pt_avg, avg_sig,       color='orange', label="DIS (TBH)", xoffset=-0.005)
plot_pos_neg(ax, runnum, Pt_avg, stat_plus_sys, color='orange', label=None,       xoffset=-0.005, alpha=0.4)
if elastic_bands_runnum:
    draw_bands(ax, elastic_bands_runnum, color='green', alpha=0.2, label="Elastic (NP)")
ax.set_xlabel("runnum", fontsize=14)
ax.set_ylabel(r"$|P_{t}|$", fontsize=14)
ax.set_ylim(0.0, 1.2)
if run_csv_nums.size > 0:
    ax.set_xlim(min(xmin, run_csv_nums.min()-10), max(xmax, run_csv_nums.max()+10))
else:
    ax.set_xlim(xmin, xmax)
ax.legend(fontsize=10)

# Right: vs run index × 5 (offsets multiplied by 5 -> DIS xoffset -1.25; Elastic band across index)
ax = axes[1]
plot_pos_neg(ax, run_index_scaled, Pt_avg, avg_sig,       color='orange', label="DIS (TBH)", xoffset=-1.25)
plot_pos_neg(ax, run_index_scaled, Pt_avg, stat_plus_sys, color='orange', label=None,       xoffset=-1.25, alpha=0.4)
if elastic_bands_index:
    draw_bands(ax, elastic_bands_index, color='green', alpha=0.2, label="Elastic (NP)")
ax.set_xlabel("run index × 5", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(0, max(run_index_scaled.max(), run_index_np_scaled.max() if run_index_np.size > 0 else 0) + 5)

plt.tight_layout()
plt.savefig("output/method_comparison.pdf")
plt.close()
print("[INFO] Saved output/method_comparison.pdf")