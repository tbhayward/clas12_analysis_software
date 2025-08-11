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

def plot_pos_neg(x, y, yerr, color, label, xoffset=0.0, alpha=1.0, zorder=None):
    """Plot |y| with filled markers for y>=0 and empty markers for y<0."""
    x = np.asarray(x) + xoffset
    y = np.asarray(y)
    yerr = None if yerr is None else np.asarray(yerr)

    pos = y >= 0
    neg = ~pos
    y_abs = np.abs(y)

    if pos.any():
        plt.errorbar(x[pos], y_abs[pos], yerr=None if yerr is None else yerr[pos],
                     fmt='o', color=color, label=label, alpha=alpha, zorder=zorder)
    if neg.any():
        plt.errorbar(x[neg], y_abs[neg], yerr=None if yerr is None else yerr[neg],
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
plt.xlim(0.0, 0.7)  # lock to requested range
plt.legend(fontsize=11)  # slightly smaller legend text
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

xmin = runnum.min() - 10
xmax = runnum.max() + 10

# ===============================
# 3) Plot Pt per run (GRV + ABD) with ±0.005 horizontal offset
#     y from 0 to 1.2, negative values drawn as empty circles
# ===============================
plt.figure(figsize=(10,6))
plot_pos_neg(runnum, Pt_grv, s_grv, color='blue', label="GRSV [hep-ph] 0011215v1", xoffset=-0.005)
plot_pos_neg(runnum, Pt_abd, s_abd, color='red',  label="ABDY [hep-ph] 0705.1553", xoffset=+0.005)
plt.xlabel("runnum", fontsize=14)
plt.ylabel(r"$P_{t}$", fontsize=14)
plt.ylim(0.0, 1.2)
plt.xlim(xmin, xmax)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/model_extractions.pdf")
plt.close()

print("[INFO] Saved output/model_extractions.pdf")

# ===============================
# 4) Plot average Pt with stat and stat+sys errors (y from 0 → 1.2)
#     negative values as empty circles
# ===============================
stat_plus_sys = np.sqrt(avg_sig**2 + avg_sys**2)

plt.figure(figsize=(10,6))
# Stat bars (solid)
plot_pos_neg(runnum, Pt_avg, avg_sig, color='black', label="Avg $P_{t}$ (stat)")
# Stat+sys bars (semi-transparent, no duplicate legend)
plot_pos_neg(runnum, Pt_avg, stat_plus_sys, color='black', label=None, alpha=0.4)
plt.xlabel("runnum", fontsize=14)
plt.ylabel(r"$P_{t}$", fontsize=14)
plt.ylim(0.0, 1.2)
plt.xlim(xmin, xmax)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/avg_pt_plot.pdf")
plt.close()

print("[INFO] Saved output/avg_pt_plot.pdf")

# ===============================
# 5) Method comparison plot (NH3-only from CSV), y from 0 → 1.2
#     negative values as empty circles; ±0.005 horizontal offset
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
            # Only true for NH3 sections in any RGC period
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
            pol_csv.append(float(row[-2]))  # polarization column
            pol_csv_err.append(float(row[-1]))  # stat err column
        except ValueError:
            continue

run_csv_nums = np.array(run_csv_nums)
pol_csv      = np.array(pol_csv)
pol_csv_err  = np.array(pol_csv_err)

plt.figure(figsize=(10,6))
# DIS (TBH) at -0.005
plot_pos_neg(runnum, Pt_avg, avg_sig, color='orange', label="DIS (TBH)", xoffset=-0.005)
# Elastic (NP) at +0.005
plot_pos_neg(run_csv_nums, pol_csv, pol_csv_err, color='green', label="Elastic (NP)", xoffset=+0.005)
plt.xlabel("runnum", fontsize=14)
plt.ylabel(r"$P_{t}$", fontsize=14)
plt.ylim(0.0, 1.2)
plt.xlim(xmin, xmax)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/method_comparison.pdf")
plt.close()

print("[INFO] Saved output/method_comparison.pdf")