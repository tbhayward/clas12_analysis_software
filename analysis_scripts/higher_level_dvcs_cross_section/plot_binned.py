import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages

# Create output directory
os.makedirs('output/plots', exist_ok=True)

# Read data from file
data_file = '/u/home/thayward/clas12_analysis_software/analysis_scripts/higher_level_dvcs_cross_section/imports/binned_rga_prl_ALU.txt'
data = []
with open(data_file, 'r') as f:
    for line in f:
        if 'Point' in line:
            parts = line.split(',')
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

# Define bin ranges
xi_bins = [(0.04, 0.1), (0.1, 0.2), (0.2, 0.45)]
Q2_bins = [(1, 2), (2, 4), (4, 8)]
t_bins = [(0.1, 0.3), (0.3, 0.6), (0.6, 3.0)]

# Function to create subplots for each bin combination
def create_subplots(x_var, y_var, fixed_var1_bins, fixed_var2_bins, 
                   x_label, fixed_var1_name, fixed_var2_name):
    fig, axs = plt.subplots(len(fixed_var1_bins), len(fixed_var2_bins), 
                           figsize=(15, 12), sharey=True)
    fig.subplots_adjust(hspace=0.15, wspace=0.1)
    
    for i, bin1 in enumerate(fixed_var1_bins):
        for j, bin2 in enumerate(fixed_var2_bins):
            ax = axs[i, j] if len(fixed_var1_bins) > 1 else axs[j]
            
            # Create mask for current bins
            if fixed_var1_name == "Q²":
                mask = (Q2 >= bin1[0]) & (Q2 < bin1[1]) & (t_minus >= bin2[0]) & (t_minus < bin2[1])
            elif fixed_var1_name == "ξ":
                mask = (xi >= bin1[0]) & (xi < bin1[1]) & (Q2 >= bin2[0]) & (Q2 < bin2[1])
            else:  # -t
                mask = (t_minus >= bin1[0]) & (t_minus < bin1[1]) & (xi >= bin2[0]) & (xi < bin2[1])
            
            if np.sum(mask) > 0:
                ax.errorbar(x_var[mask], y_var[mask], yerr=sigma[mask], 
                           fmt='o', color='b', ms=6, capsize=4, alpha=0.8)
            
            # Set titles for first row and first column
            if i == 0:
                ax.set_title(f"{fixed_var2_name} = {bin2[0]}-{bin2[1]}", fontsize=10)
            if j == 0:
                ax.set_ylabel(f"{fixed_var1_name} = {bin1[0]}-{bin1[1]}\n$A_{{LU}}$", fontsize=10)
            
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_xlabel(x_label, fontsize=10)
    
    return fig

# Create PDF for ξ dependence
with PdfPages('output/plots/xi_dependence.pdf') as pdf:
    fig = create_subplots(
        x_var=xi, 
        y_var=Amp,
        fixed_var1_bins=Q2_bins,
        fixed_var2_bins=t_bins,
        x_label="ξ",
        fixed_var1_name="Q² (GeV²)",
        fixed_var2_name="-t (GeV²)"
    )
    fig.suptitle("Beam Spin Asymmetry vs ξ", fontsize=16)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# Create PDF for -t dependence
with PdfPages('output/plots/t_dependence.pdf') as pdf:
    fig = create_subplots(
        x_var=t_minus, 
        y_var=Amp,
        fixed_var1_bins=xi_bins,
        fixed_var2_bins=Q2_bins,
        x_label="$-t$ (GeV$^{2}$)",
        fixed_var1_name="ξ",
        fixed_var2_name="Q² (GeV²)"
    )
    fig.suptitle("Beam Spin Asymmetry vs $-t$", fontsize=16)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# Create PDF for Q² dependence
with PdfPages('output/plots/Q2_dependence.pdf') as pdf:
    fig = create_subplots(
        x_var=Q2, 
        y_var=Amp,
        fixed_var1_bins=xi_bins,
        fixed_var2_bins=t_bins,
        x_label="$Q^{2}$ (GeV$^{2}$)",
        fixed_var1_name="ξ",
        fixed_var2_name="-t (GeV²)"
    )
    fig.suptitle("Beam Spin Asymmetry vs $Q^{2}$", fontsize=16)
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)