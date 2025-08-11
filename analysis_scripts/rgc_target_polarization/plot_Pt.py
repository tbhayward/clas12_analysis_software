#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

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
# (Format: Run Pt_GRV sigma_GRV Pt_ABD sigma_ABD Pt_avg avg_sig avg_sys)
# ===============================
extractions_file = "output/Pt_by_run.txt"
runnum, Pt_grv, s_grv, Pt_abd, s_abd, Pt_avg, avg_sig, avg_sys = [], [], [], [], [], [], [], []

with open(extractions_file) as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith("#") or line.startswith("Run"):
            continue
        parts = line.split()
        if len(parts) < 8:
            continue
        try:
            runnum.append(int(parts[0]))
            Pt_grv.append(float(parts[1]))
            s_grv.append(float(parts[2]))
            Pt_abd.append(float(parts[3]))
            s_abd.append(float(parts[4]))
            Pt_avg.append(float(parts[5]))
            avg_sig.append(float(parts[6]))
            avg_sys.append(float(parts[7]))
        except ValueError:
            continue

runnum = np.array(runnum)
Pt_grv = np.array(Pt_grv)
s_grv  = np.array(s_grv)
Pt_abd = np.array(Pt_abd)
s_abd  = np.array(s_abd)
Pt_avg = np.array(Pt_avg)
avg_sig= np.array(avg_sig)
avg_sys= np.array(avg_sys)

# sort by runnum (for consistent left/right mapping)
order = np.argsort(runnum)
runnum   = runnum[order]
Pt_grv   = Pt_grv[order]
s_grv    = s_grv[order]
Pt_abd   = Pt_abd[order]
s_abd    = s_abd[order]
Pt_avg   = Pt_avg[order]
avg_sig  = avg_sig[order]
avg_sys  = avg_sys[order]

xmin = runnum.min() - 10
xmax = runnum.max() + 10
run_index = np.arange(1, len(runnum)+1)

# ===============================
# 3) 1x2: Pt per run (GRV + ABD), with stat and stat+sys
# ===============================
tot_grv = np.sqrt(s_grv**2 + avg_sys**2)
tot_abd = np.sqrt(s_abd**2 + avg_sys**2)

fig, axes = plt.subplots(1, 2, figsize=(16,6), sharey=True)

# Left: vs runnum
ax = axes[0]
plot_pos_neg(ax, runnum, Pt_grv, s_grv, color='blue', label="GRSV [stat]", xoffset=-0.005)
plot_pos_neg(ax, runnum, Pt_grv, tot_grv, color='blue', label="GRSV [stat+sys]", xoffset=-0.005, alpha=0.4)
plot_pos_neg(ax, runnum, Pt_abd, s_abd, color='red', label="ABDY [stat]", xoffset=+0.005)
plot_pos_neg(ax, runnum, Pt_abd, tot_abd, color='red', label="ABDY [stat+sys]", xoffset=+0.005, alpha=0.4)
ax.set_xlabel("runnum", fontsize=14)
ax.set_ylabel(r"$P_{t}$", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(xmin, xmax)
ax.legend(fontsize=10)

# Right: vs run index
ax = axes[1]
plot_pos_neg(ax, run_index, Pt_grv, s_grv, color='blue', label="GRSV [stat]", xoffset=-0.05)
plot_pos_neg(ax, run_index, Pt_grv, tot_grv, color='blue', label="GRSV [stat+sys]", xoffset=-0.05, alpha=0.4)
plot_pos_neg(ax, run_index, Pt_abd, s_abd, color='red', label="ABDY [stat]", xoffset=+0.05)
plot_pos_neg(ax, run_index, Pt_abd, tot_abd, color='red', label="ABDY [stat+sys]", xoffset=+0.05, alpha=0.4)
ax.set_xlabel("run index", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(0.5, len(run_index)+0.5)

plt.tight_layout()
plt.savefig("output/model_extractions.pdf")
plt.close()
print("[INFO] Saved output/model_extractions.pdf")

# ===============================
# 4) 1x2: Average Pt with stat and stat+sys
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

# Right: vs run index
ax = axes[1]
plot_pos_neg(ax, run_index, Pt_avg, avg_sig,       color='black', label="Avg $P_{t}$ (stat)")
plot_pos_neg(ax, run_index, Pt_avg, stat_plus_sys, color='black', label="Avg $P_{t}$ (stat+sys)", alpha=0.4)
ax.set_xlabel("run index", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(0.5, len(run_index)+0.5)

plt.tight_layout()
plt.savefig("output/avg_pt_plot.pdf")
plt.close()
print("[INFO] Saved output/avg_pt_plot.pdf")

# ===============================
# 5) 1x2: Method comparison (NH3-only from CSV)
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
            pol_csv.append(float(row[-2]))      # polarization column
            pol_csv_err.append(float(row[-1]))  # stat err column
        except ValueError:
            continue

run_csv_nums = np.array(run_csv_nums)
pol_csv      = np.array(pol_csv)
pol_csv_err  = np.array(pol_csv_err)

# Sort Elastic arrays by run
order_np = np.argsort(run_csv_nums)
run_csv_nums = run_csv_nums[order_np]
pol_csv      = pol_csv[order_np]
pol_csv_err  = pol_csv_err[order_np]
run_index_np = np.arange(1, len(run_csv_nums)+1)

fig, axes = plt.subplots(1, 2, figsize=(16,6), sharey=True)

# Left: vs runnum (offset ±0.005)
ax = axes[0]
plot_pos_neg(ax, runnum, Pt_avg, avg_sig,           color='orange', label="DIS (TBH)", xoffset=-0.005)
plot_pos_neg(ax, runnum, Pt_avg, stat_plus_sys,     color='orange', label=None,       xoffset=-0.005, alpha=0.4)
plot_pos_neg(ax, run_csv_nums, pol_csv, pol_csv_err, color='green',  label="Elastic (NP)", xoffset=+0.005)
ax.set_xlabel("runnum", fontsize=14)
ax.set_ylabel(r"$P_{t}$", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(min(xmin, run_csv_nums.min()-10 if len(run_csv_nums) else xmin),
            max(xmax, run_csv_nums.max()+10 if len(run_csv_nums) else xmax))
ax.legend(fontsize=10)

# Right: vs run index (independent indices, offset ±0.05)
ax = axes[1]
plot_pos_neg(ax, run_index,   Pt_avg, avg_sig,       color='orange', label="DIS (TBH)", xoffset=-0.05)
plot_pos_neg(ax, run_index,   Pt_avg, stat_plus_sys, color='orange', label=None,       xoffset=-0.05, alpha=0.4)
plot_pos_neg(ax, run_index_np, pol_csv, pol_csv_err, color='green',  label="Elastic (NP)", xoffset=+0.05)
ax.set_xlabel("run index", fontsize=14)
ax.set_ylim(0.0, 1.2)
ax.set_xlim(0.5, max(len(run_index), len(run_index_np)) + 0.5)

plt.tight_layout()
plt.savefig("output/method_comparison.pdf")
plt.close()
print("[INFO] Saved output/method_comparison.pdf")