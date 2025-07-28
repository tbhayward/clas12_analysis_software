import numpy as np
import matplotlib.pyplot as plt
import os

# Create output directory if it doesn't exist
os.makedirs('output/plots', exist_ok=True)

# Read data from file
data_file = '/u/home/thayward/clas12_analysis_software/analysis_scripts/higher_level_dvcs_cross_section/imports/binned_rga_prl_ALU.txt'
data = []
with open(data_file, 'r') as f:
    for line in f:
        if 'Point' in line:
            parts = line.split(',')
            # Extract relevant values
            xB = float(parts[1].split('=')[1])
            Q2 = float(parts[2].split('=')[1])
            t_minus = float(parts[3].split('=')[1])
            Amp = float(parts[4].split('=')[1])
            sigma = float(parts[5].split('=')[1])
            data.append((xB, Q2, t_minus, Amp, sigma))

# Convert to numpy array and compute xi
data = np.array(data)
xB = data[:, 0]
xi = xB / (2 - xB)
Q2 = data[:, 1]
t_minus = data[:, 2]
Amp = data[:, 3]
sigma = data[:, 4]

# Define binning scheme
xi_bins = [(0.04, 0.1), (0.1, 0.2), (0.2, 0.45)]
Q2_bins = [(1, 2), (2, 4), (4, 8)]
t_bins = [(0.1, 0.3), (0.3, 0.6), (0.6, 3.0)]

# Create figure with 3 subplots
fig, axs = plt.subplots(1, 3, figsize=(18, 6))
plt.subplots_adjust(wspace=0.3)

# Color and marker scheme for bins
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
markers = ['o', 's', '^']

# Plot 1: Amp vs xi (grouped by Q2 and t)
for i, q_bin in enumerate(Q2_bins):
    for j, t_bin in enumerate(t_bins):
        mask = (Q2 >= q_bin[0]) & (Q2 < q_bin[1]) & (t_minus >= t_bin[0]) & (t_minus < t_bin[1])
        if np.any(mask):
            label = f'Q²: {q_bin[0]}-{q_bin[1]} GeV²\n-t: {t_bin[0]}-{t_bin[1]} GeV²'
            axs[0].errorbar(xi[mask], Amp[mask], yerr=sigma[mask], 
                            fmt=markers[j], color=colors[i], mec='k', ms=8, capsize=4, label=label)

axs[0].set_xlabel('ξ', fontsize=14)
axs[0].set_ylabel('$A_{LU}$', fontsize=14)
axs[0].grid(True, linestyle='--', alpha=0.7)
axs[0].legend(loc='upper right', fontsize=9, framealpha=0.95)

# Plot 2: Amp vs -t (grouped by xi and Q2)
for i, x_bin in enumerate(xi_bins):
    for j, q_bin in enumerate(Q2_bins):
        mask = (xi >= x_bin[0]) & (xi < x_bin[1]) & (Q2 >= q_bin[0]) & (Q2 < q_bin[1])
        if np.any(mask):
            label = f'ξ: {x_bin[0]:.2f}-{x_bin[1]:.2f}\nQ²: {q_bin[0]}-{q_bin[1]} GeV²'
            axs[1].errorbar(t_minus[mask], Amp[mask], yerr=sigma[mask], 
                            fmt=markers[j], color=colors[i], mec='k', ms=8, capsize=4, label=label)

axs[1].set_xlabel('$-t$ (GeV$^{2}$)', fontsize=14)
axs[1].grid(True, linestyle='--', alpha=0.7)
axs[1].legend(loc='upper right', fontsize=9, framealpha=0.95)

# Plot 3: Amp vs Q² (grouped by xi and t)
for i, x_bin in enumerate(xi_bins):
    for j, t_bin in enumerate(t_bins):
        mask = (xi >= x_bin[0]) & (xi < x_bin[1]) & (t_minus >= t_bin[0]) & (t_minus < t_bin[1])
        if np.any(mask):
            label = f'ξ: {x_bin[0]:.2f}-{x_bin[1]:.2f}\n-t: {t_bin[0]}-{t_bin[1]} GeV²'
            axs[2].errorbar(Q2[mask], Amp[mask], yerr=sigma[mask], 
                            fmt=markers[j], color=colors[i], mec='k', ms=8, capsize=4, label=label)

axs[2].set_xlabel('$Q^{2}$ (GeV$^{2}$)', fontsize=14)
axs[2].grid(True, linestyle='--', alpha=0.7)
axs[2].legend(loc='upper right', fontsize=9, framealpha=0.95)

# Set common y-axis limits
all_amps = np.concatenate([Amp + sigma, Amp - sigma])
y_min, y_max = np.min(all_amps), np.max(all_amps)
y_padding = 0.1 * (y_max - y_min)
for ax in axs:
    ax.set_ylim(y_min - y_padding, y_max + y_padding)

plt.tight_layout()
plt.savefig('output/plots/binned_rga_results.pdf', bbox_inches='tight')
plt.close()