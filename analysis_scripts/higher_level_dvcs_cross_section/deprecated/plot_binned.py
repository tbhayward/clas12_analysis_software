import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import MaxNLocator

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

# Define 5 bins for each variable based on data distribution
xi_edges = np.percentile(xi, [0, 20, 40, 60, 80, 100])
Q2_edges = np.percentile(Q2, [0, 20, 40, 60, 80, 100])
t_edges = np.percentile(t_minus, [0, 20, 40, 60, 80, 100])

# Create bin ranges from edges
xi_bins = [(xi_edges[i], xi_edges[i+1]) for i in range(5)]
Q2_bins = [(Q2_edges[i], Q2_edges[i+1]) for i in range(5)]
t_bins = [(t_edges[i], t_edges[i+1]) for i in range(5)]

# Print bin definitions for verification
print("\nξ bins:")
for i, bin in enumerate(xi_bins):
    print(f"Bin {i+1}: {bin[0]:.4f} - {bin[1]:.4f}")
    
print("\nQ² bins:")
for i, bin in enumerate(Q2_bins):
    print(f"Bin {i+1}: {bin[0]:.4f} - {bin[1]:.4f} GeV²")
    
print("\n-t bins:")
for i, bin in enumerate(t_bins):
    print(f"Bin {i+1}: {bin[0]:.4f} - {bin[1]:.4f} GeV²")

# Function to create 5x5 subplot grid
def create_5x5_grid(x_var, y_var, row_bins, col_bins, 
                   x_label, row_label, col_label, title):
    fig, axs = plt.subplots(5, 5, figsize=(25, 20), 
                           sharey=True, squeeze=True)
    fig.subplots_adjust(hspace=0.3, wspace=0.2)
    fig.suptitle(title, fontsize=24, y=0.98)
    
    # Plot data in each subplot
    total_points = 0
    for i in range(5):
        for j in range(5):
            ax = axs[i, j]
            
            # Create mask for current bins
            if row_label == "Q² (GeV²)":
                row_mask = (Q2 >= row_bins[i][0]) & (Q2 < row_bins[i][1])
                col_mask = (t_minus >= col_bins[j][0]) & (t_minus < col_bins[j][1])
            elif row_label == "ξ":
                row_mask = (xi >= row_bins[i][0]) & (xi < row_bins[i][1])
                col_mask = (Q2 >= col_bins[j][0]) & (Q2 < col_bins[j][1])
            else:  # -t
                row_mask = (t_minus >= row_bins[i][0]) & (t_minus < row_bins[i][1])
                col_mask = (xi >= col_bins[j][0]) & (xi < col_bins[j][1])
                
            mask = row_mask & col_mask
            
            # Plot data if available
            if np.any(mask):
                x_vals = x_var[mask]
                y_vals = y_var[mask]
                y_err = sigma[mask]
                
                # Plot with error bars
                ax.errorbar(x_vals, y_vals, yerr=y_err, fmt='o', color='blue', 
                           ms=8, capsize=5, alpha=0.8, markeredgecolor='black')
                
                # Add bin information
                bin_info = (f"{row_label}: {row_bins[i][0]:.3f}-{row_bins[i][1]:.3f}\n"
                            f"{col_label}: {col_bins[j][0]:.3f}-{col_bins[j][1]:.3f}\n"
                            f"N = {len(x_vals)}")
                ax.text(0.05, 0.95, bin_info, transform=ax.transAxes, 
                       fontsize=12, va='top', bbox=dict(facecolor='white', alpha=0.8))
                
                total_points += len(x_vals)
            
            # Configure axes
            ax.grid(True, linestyle='--', alpha=0.5)
            ax.set_xlabel(x_label, fontsize=14)
            ax.xaxis.set_major_locator(MaxNLocator(5))
            if j == 0:
                ax.set_ylabel('$A_{LU}$', fontsize=14)
    
    print(f"Plotted {total_points} points for {title}")
    return fig

# Create PDF for ξ dependence
with PdfPages('output/plots/xi_dependence_5x5.pdf') as pdf:
    fig = create_5x5_grid(
        x_var=xi, 
        y_var=Amp,
        row_bins=Q2_bins,
        col_bins=t_bins,
        x_label="ξ",
        row_label="Q² (GeV²)",
        col_label="-t (GeV²)",
        title="Beam Spin Asymmetry vs ξ (Grouped by Q² and -t)"
    )
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# Create PDF for -t dependence
with PdfPages('output/plots/t_dependence_5x5.pdf') as pdf:
    fig = create_5x5_grid(
        x_var=t_minus, 
        y_var=Amp,
        row_bins=xi_bins,
        col_bins=Q2_bins,
        x_label="$-t$ (GeV$^{2}$)",
        row_label="ξ",
        col_label="Q² (GeV²)",
        title="Beam Spin Asymmetry vs $-t$ (Grouped by ξ and Q²)"
    )
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)

# Create PDF for Q² dependence
with PdfPages('output/plots/Q2_dependence_5x5.pdf') as pdf:
    fig = create_5x5_grid(
        x_var=Q2, 
        y_var=Amp,
        row_bins=xi_bins,
        col_bins=t_bins,
        x_label="$Q^{2}$ (GeV$^{2}$)",
        row_label="ξ",
        col_label="-t (GeV²)",
        title="Beam Spin Asymmetry vs $Q^{2}$ (Grouped by ξ and -t)"
    )
    pdf.savefig(fig, bbox_inches='tight')
    plt.close(fig)