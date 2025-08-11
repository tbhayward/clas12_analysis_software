#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import csv
import os

# Ensure output directory exists
os.makedirs("output", exist_ok=True)

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
# 3) Plot Pt per run (GRV + ABD) with small horizontal offset
# ===============================
plt.figure(figsize=(10,6))
plt.errorbar(runnum - 0.005, Pt_grv, yerr=s_grv, fmt='o', color='blue', label="GRSV [hep-ph] 0011215v1")
plt.errorbar(runnum + 0.005, Pt_abd, yerr=s_abd, fmt='o', color='red',  label="ABDY [hep-ph] 0705.1553")
plt.xlabel("runnum", fontsize=14)
plt.ylabel(r"$P_{t}$", fontsize=14)
plt.ylim(-1, 1)
plt.xlim(xmin, xmax)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/model_extractions.pdf")
plt.close()

print("[INFO] Saved output/model_extractions.pdf")

# ===============================
# 4) Plot average Pt with stat and stat+sys errors
# ===============================
stat_plus_sys = np.sqrt(avg_sig**2 + avg_sys**2)

plt.figure(figsize=(10,6))
plt.errorbar(runnum, Pt_avg, yerr=avg_sig,        fmt='o', color='black', label="Avg $P_{t}$ (stat)")
plt.errorbar(runnum, Pt_avg, yerr=stat_plus_sys,  fmt='o', color='black', alpha=0.4, label="Avg $P_{t}$ (stat+sys)")
plt.xlabel("runnum", fontsize=14)
plt.ylabel(r"$P_{t}$", fontsize=14)
plt.ylim(-1, 1)
plt.xlim(xmin, xmax)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/avg_pt_plot.pdf")
plt.close()

print("[INFO] Saved output/avg_pt_plot.pdf")

# ===============================
# 5) Method comparison plot
# ===============================
# Load clas12_run_info.csv (runnum, ..., pol_s, pol_e) – we use last two cols
run_csv_file = "/u/home/thayward/clas12_analysis_software/analysis_scripts/asymmetry_extraction/imports/clas12_run_info.csv"

run_csv_nums = []
pol_csv      = []
pol_csv_err  = []

with open(run_csv_file, newline='') as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        if not row or row[0].startswith("#"):
            continue
        try:
            run_csv_nums.append(int(row[0]))
            pol_csv.append(float(row[-2]))   # polarization column
            pol_csv_err.append(float(row[-1]))  # stat err column
        except ValueError:
            continue

run_csv_nums = np.array(run_csv_nums)
pol_csv      = np.array(pol_csv)
pol_csv_err  = np.array(pol_csv_err)

plt.figure(figsize=(10,6))
plt.errorbar(runnum,       Pt_avg,  yerr=avg_sig,      fmt='o', color='orange', label="DIS (TBH)")
plt.errorbar(run_csv_nums, pol_csv, yerr=pol_csv_err,  fmt='o', color='green',  label="Elastic (NP)")
plt.xlabel("runnum", fontsize=14)
plt.ylabel(r"$P_{t}$", fontsize=14)
plt.ylim(-1, 1)
plt.xlim(xmin, xmax)
plt.legend(fontsize=11)
plt.tight_layout()
plt.savefig("output/method_comparison.pdf")
plt.close()

print("[INFO] Saved output/method_comparison.pdf")