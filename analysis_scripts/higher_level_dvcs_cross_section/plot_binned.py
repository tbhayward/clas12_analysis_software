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

# Print data statistics for verification
print(f"Total data points: {len(xi)}")
print(f"xi range: {xi.min():.4f} - {xi.max():.4f}")
print(f"Q² range: {Q2.min():.4f} - {Q2.max():.4f} GeV²")
print(f"-t range: {t_minus.min():.4f} - {t_minus.max():.4f} GeV²")

# Define bin ranges based on data distribution
xi_bins = [(0.04, 0.11), (0.11, 0.22), (0.22, 0.45)]
Q2_bins = [(1.0, 1.8), (1.8, 2.5), (2.5, 8.0)]
t_bins = [(0.14, 0.25), (0.25, 0.45), (0.45, 1.2), (1.2, 3.0)]

# Function to create subplots for each bin combination
def create_subplots(x_var, y_var, fixed_var1_bins, fixed_var2_bins, 
                   x_label, fixed_var1_name, fixed_var2_name, title):
    # Count total data points that will be plotted
    total_points = 0
    
    # Calculate how many subplots we need (only bins with data)
    valid_subplots = []
    for i, bin1 in enumerate(fixed_var1_bins):
        for j, bin2 in enumerate(fixed_var2_bins):
            # Create mask for current bins
            if fixed_var1_name == "Q² (GeV²)":
                mask = (Q2 >= bin1[0]) & (Q2 < bin1[1]) & (t_minus >= bin2[0]) & (t_minus < bin2[1])
            elif fixed_var1_name == "ξ":
                mask = (xi >= bin1[0]) & (xi < bin1[1]) & (Q2 >= bin2[0]) & (Q2 < bin2[1])
            else:  # -t
                mask = (t_minus >= bin1[0]) & (t_minus < bin1[1]) & (xi >= bin2[0]) & (xi < bin2[1])
            
            if np.sum(mask) > 0:
                valid_subplots.append((i, j, bin1, bin2, mask))
                total_points += np.sum(mask)
    
    # Create figure with appropriate grid
    rows = len(fixed_var1_bins)
    cols = len(fixed_var2_bins)
    fig, axs = plt.subplots(rows, cols, figsize=(4*cols, 4*rows), 
                           sharey=True, squeeze=False)
    fig.subplots_adjust(hspace=0.25, wspace=0.15)
    fig.suptitle(title, fontsize=16, y=0.98)
    
    # Plot data in valid subplots
    for i, j, bin1, bin2, mask in valid_subplots:
        ax = axs[i, j]
        x_vals = x_var[mask]
        y_vals = y_var[mask]
        y_err = sigma[mask]
        
        # Plot with error bars
        ax.errorbar(x_vals, y_vals, yerr=y_err, fmt='o', color='b', 
                   ms=6, capsize=4, alpha=0.8)
        
        # Add bin information and count
        bin1_label = f"{fixed_var1_name}: {bin1[0]:.2f}-{bin1[1]:.2f}"
        bin2_label = f"{fixed_var2_name}: {bin2[0]:.2f}-{bin2[1]:.2f}"
        count_label = f"N = {len(x_vals)}"
        ax.text(0.05, 0.95, bin1_label, transform=ax.transAxes, 
               fontsize=10, va='top')
        ax.text(0.05, 0.88, bin2_label, transform=ax.transAxes, 
               fontsize=10, va='top')
        ax.text(0.95, 0.95, count_label, transform=ax.transAxes, 
               fontsize=10, ha='right', va='top')
        
        ax.grid(True, linestyle='--', alpha=0.5)
        ax.set_xlabel(x_label, fontsize=12)
        if j == 0:
            ax.set_ylabel('$A_{LU}$', fontsize=12)
    
    # Hide empty subplots
    for i in range(rows):
        for j in range(cols):
            if (i, j) not in [(idx[0], idx[1]) for idx in valid_subplots]:
                axs[i, j].axis('off')
    
    print(f"Plotted {total_points} points for {title}")
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
        fixed_var2_name="-t (GeV²)",
        title="Beam Spin Asymmetry vs ξ"
    )
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
        fixed_var2_name="Q² (GeV²)",
        title="Beam Spin Asymmetry vs $-t$"
    )
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
        fixed_var2_name="-t (GeV²)",
        title="Beam Spin Asymmetry vs $Q^{2}$"
    )
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)