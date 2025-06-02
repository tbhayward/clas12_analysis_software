#!/usr/bin/env python3

import sys
import os
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------------------------------------------------------
# Parse command-line argument for run period
# -----------------------------------------------------------------------------
if len(sys.argv) < 2:
    print("Usage: python plot_spin_asymmetry.py <run_period_number>")
    print("  where <run_period_number> is 1 (RGCSu22), 2 (RGCFa22), or 3 (RGCSp23)")
    sys.exit(1)
# endif

period = sys.argv[1]

if period not in ("1", "2", "3"):
    print(f"Invalid run period '{period}'. Must be '1', '2', or '3'.")
    sys.exit(1)
# endif

# -----------------------------------------------------------------------------
# Define data for RGCSu22 (period == "1"); placeholders for others
# -----------------------------------------------------------------------------
if period == "1":
    run_period_name = "RGCSu22"

    # Each entry: [xB, central_value, uncertainty]
    enpichi2FitsALUsinphi_RGCSu22 = [
        [0.093969228, 1.121102910, 0.103682047],
        [0.168352864, 0.075343496, 0.011124971],
        [0.255015669, 0.110872309, 0.007668837],
        [0.348460186, 0.136797553, 0.008004000],
        [0.441445130, 0.122327908, 0.010538986],
        [0.535418362, 0.095188929, 0.018040875]
    ]

    enpichi2FitsAULsinphi_RGCSu22 = [
        [0.093969228, -0.338857579, 0.194599332],
        [0.168352864, -0.051397484, 0.010823311],
        [0.255015669, -0.003922942, 0.005337434],
        [0.348460186, 0.035992215, 0.005302819],
        [0.441445130, 0.039496051, 0.006732558],
        [0.535418362, 0.031684185, 0.011741457]
    ]

    enpichi2FitsAULsin2phi_RGCSu22 = [
        [0.093969228, -0.933145689, 0.505068832],
        [0.168352864, -0.020302810, 0.022454560],
        [0.255015669, -0.068209117, 0.012126491],
        [0.348460186, -0.086572841, 0.011778306],
        [0.441445130, -0.051343918, 0.014025767],
        [0.535418362, -0.036623704, 0.024819806]
    ]

    enpichi2FitsALL_RGCSu22 = [
        [0.093969228, -0.226411668, 0.322623368],
        [0.168352864, 0.191460695, 0.023617359],
        [0.255015669, 0.302080680, 0.022686092],
        [0.348460186, 0.434065099, 0.029778569],
        [0.441445130, 0.522671691, 0.036871963],
        [0.535418362, 0.634253359, 0.048909323]
    ]

    enpichi2FitsALLcosphi_RGCSu22 = [
        [0.093969228, -0.622700113, 0.579882138],
        [0.168352864, 0.187035564, 0.037593926],
        [0.255015669, 0.067043637, 0.035191562],
        [0.348460186, -0.028525676, 0.044653070],
        [0.441445130, -0.024721894, 0.056791830],
        [0.535418362, -0.041561041, 0.074395483]
    ]
# endif

elif period == "2":
    print("Data for RGCFa22 not yet available. Please add enpichi2Fits... arrays.")
    sys.exit(1)
# endif

elif period == "3":
    print("Data for RGCSp23 not yet available. Please add enpichi2Fits... arrays.")
    sys.exit(1)
# endif

# -----------------------------------------------------------------------------
# Convert lists to numpy arrays
# -----------------------------------------------------------------------------
# ALU sinφ
x_ALUsinphi = np.array([row[0] for row in enpichi2FitsALUsinphi_RGCSu22])
y_ALUsinphi = np.array([row[1] for row in enpichi2FitsALUsinphi_RGCSu22])
err_ALUsinphi = np.array([row[2] for row in enpichi2FitsALUsinphi_RGCSu22])

# AUL sinφ
x_AULsinphi = np.array([row[0] for row in enpichi2FitsAULsinphi_RGCSu22])
y_AULsinphi = np.array([row[1] for row in enpichi2FitsAULsinphi_RGCSu22])
err_AULsinphi = np.array([row[2] for row in enpichi2FitsAULsinphi_RGCSu22])

# AUL sin2φ
x_AULsin2phi = np.array([row[0] for row in enpichi2FitsAULsin2phi_RGCSu22])
y_AULsin2phi = np.array([row[1] for row in enpichi2FitsAULsin2phi_RGCSu22])
err_AULsin2phi = np.array([row[2] for row in enpichi2FitsAULsin2phi_RGCSu22])

# ALL n=0
x_ALL_n0 = np.array([row[0] for row in enpichi2FitsALL_RGCSu22])
y_ALL_n0 = np.array([row[1] for row in enpichi2FitsALL_RGCSu22])
err_ALL_n0 = np.array([row[2] for row in enpichi2FitsALL_RGCSu22])

# ALL n=1 (cosφ)
x_ALLcosphi = np.array([row[0] for row in enpichi2FitsALLcosphi_RGCSu22])
y_ALLcosphi = np.array([row[1] for row in enpichi2FitsALLcosphi_RGCSu22])
err_ALLcosphi = np.array([row[2] for row in enpichi2FitsALLcosphi_RGCSu22])

# -----------------------------------------------------------------------------
# Create output directory if it does not exist
# -----------------------------------------------------------------------------
out_dir = "output"
if not os.path.isdir(out_dir):
    os.makedirs(out_dir, exist_ok=True)
# endif

# -----------------------------------------------------------------------------
# Plotting
# -----------------------------------------------------------------------------
plt.figure(figsize=(15, 5))

# Subplot 1: F_LU^{sinφ} / F_UU
ax1 = plt.subplot(1, 3, 1)
ax1.errorbar(
    x_ALUsinphi,
    y_ALUsinphi,
    yerr=err_ALUsinphi,
    fmt="o",
    color="tab:blue",
    ecolor="tab:gray",
    capsize=3,
    label=None
)
ax1.set_xlim(0, 0.7)
ax1.set_ylim(-0.2, 0.2)
ax1.set_xlabel(r"$x_{B}$")
ax1.set_ylabel(r"$F_{LU}^{\sin\phi}/F_{UU}$")
ax1.grid(True, linestyle="--", alpha=0.6)

# Subplot 2: F_UL^{sin n φ} / F_UU, n=1,2
ax2 = plt.subplot(1, 3, 2, sharex=ax1)
ax2.errorbar(
    x_AULsinphi,
    y_AULsinphi,
    yerr=err_AULsinphi,
    fmt="s",
    color="tab:orange",
    ecolor="tab:gray",
    capsize=3,
    label="n=1"
)
ax2.errorbar(
    x_AULsin2phi,
    y_AULsin2phi,
    yerr=err_AULsin2phi,
    fmt="^",
    color="tab:green",
    ecolor="tab:gray",
    capsize=3,
    label="n=2"
)
ax2.set_xlim(0, 0.7)
ax2.set_xlabel(r"$x_{B}$")
ax2.set_ylabel(r"$F_{UL}^{\sin\,n\phi}/F_{UU}$")
ax2.legend(frameon=False)
ax2.grid(True, linestyle="--", alpha=0.6)

# Subplot 3: F_LL^{cos n φ} / F_UU, n=0,1
ax3 = plt.subplot(1, 3, 3, sharex=ax1)
ax3.errorbar(
    x_ALL_n0,
    y_ALL_n0,
    yerr=err_ALL_n0,
    fmt="D",
    color="tab:red",
    ecolor="tab:gray",
    capsize=3,
    label="n=0"
)
ax3.errorbar(
    x_ALLcosphi,
    y_ALLcosphi,
    yerr=err_ALLcosphi,
    fmt="v",
    color="tab:purple",
    ecolor="tab:gray",
    capsize=3,
    label="n=1"
)
ax3.set_xlim(0, 0.7)
ax3.set_ylim(-0.8, 0.8)
ax3.set_xlabel(r"$x_{B}$")
ax3.set_ylabel(r"$F_{LL}^{\cos\,n\phi}/F_{UU}$")
ax3.legend(frameon=False)
ax3.grid(True, linestyle="--", alpha=0.6)

plt.tight_layout()

# Save the figure
output_filename = os.path.join(out_dir, f"rgc_enpi+_{run_period_name}.pdf")
plt.savefig(output_filename)

print(f"Plot saved to '{output_filename}'")