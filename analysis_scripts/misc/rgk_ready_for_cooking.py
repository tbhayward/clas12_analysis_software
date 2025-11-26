import uproot
import numpy as np
import matplotlib.pyplot as plt
import os

# Create output directory if it doesn't exist
output_dir = "output/RGK_ready_for_cooking/"
os.makedirs(output_dir, exist_ok=True)

# File paths and labels
files = [
    "/volatile/clas12/thayward/RGK_pass2/rga_fa18_out_epi+X.root",
    "/volatile/clas12/thayward/RGK_pass2/rgk_pass4_epi+X.root",
    "/volatile/clas12/thayward/RGK_pass2/rgk_pass3_epi+X.root"
]
labels = ["RGA 5-pass", "RGK 4-pass", "RGK 3-pass"]
colors = ['red', 'blue', 'green']

# Plot 1: 2D distribution of epsilon vs Q2 for 0.20 < x < 0.22
print("Creating 2D distribution plot...")
fig, ax = plt.subplots(figsize=(10, 8))

for file_path, label, color in zip(files, labels, colors):
    with uproot.open(file_path) as file:
        tree = file["PhysicsEvents"]
        
        # Load branches
        DepA = tree["DepA"].array(library="np")
        DepB = tree["DepB"].array(library="np")
        Q2 = tree["Q2"].array(library="np")
        x = tree["x"].array(library="np")
        
        # Calculate epsilon
        epsilon = DepB / DepA
        
        # Apply x cut
        mask = (x > 0.20) & (x < 0.22)
        epsilon_cut = epsilon[mask]
        Q2_cut = Q2[mask]
        
        # Plot as scatter
        ax.scatter(Q2_cut, epsilon_cut, c=color, label=label, alpha=0.5, s=1)

ax.set_xlabel(r"$Q^{2}$ (GeV$^{2}$)", fontsize=14)
ax.set_ylabel(r"$\epsilon$", fontsize=14)
ax.set_xlim(1, 6)
ax.set_ylim(0.2, 1)
ax.set_title(r"$0.20 < x_{B} < 0.22$", fontsize=16)
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(os.path.join(output_dir, "epsilon_vs_Q2.png"), dpi=300)
plt.close()
print(f"Saved: {os.path.join(output_dir, 'epsilon_vs_Q2.png')}")

# Plot 2: Beam spin asymmetry for all datasets on one plot
print("\nCreating beam spin asymmetry plot...")
from scipy.optimize import curve_fit

# Define the fit function: A*sin(phi)
def sin_fit(phi, A):
    return A * np.sin(phi)

fig, ax = plt.subplots(figsize=(10, 8))

for idx, (file_path, label, color) in enumerate(zip(files, labels, colors)):
    with uproot.open(file_path) as file:
        tree = file["PhysicsEvents"]
        
        # Load branches
        x = tree["x"].array(library="np")
        phi = tree["phi"].array(library="np")
        helicity = tree["helicity"].array(library="np")
        DepA = tree["DepA"].array(library="np")
        DepW = tree["DepW"].array(library="np")
        
        # Apply x cut
        mask = (x > 0.1) & (x < 0.15)
        phi_cut = phi[mask]
        helicity_cut = helicity[mask]
        DepA_cut = DepA[mask]
        DepW_cut = DepW[mask]
        
        # Create 12 bins from 0 to 2pi
        n_bins = 12
        phi_bins = np.linspace(0, 2*np.pi, n_bins + 1)
        phi_centers = (phi_bins[:-1] + phi_bins[1:]) / 2
        
        asymmetry = []
        asymmetry_err = []
        
        for i in range(n_bins):
            # Mask for this phi bin
            bin_mask = (phi_cut >= phi_bins[i]) & (phi_cut < phi_bins[i+1])
            
            # Count events with helicity +1 and -1
            N_plus = np.sum((helicity_cut[bin_mask] == 1))
            N_minus = np.sum((helicity_cut[bin_mask] == -1))
            
            # Calculate asymmetry
            if (N_plus + N_minus) > 0:
                asym = (N_plus - N_minus) / (N_plus + N_minus)
                
                # Calculate uncertainty using error propagation
                # d(A)/dN+ = 2*N- / (N+ + N-)^2
                # d(A)/dN- = -2*N+ / (N+ + N-)^2
                # sigma_A = sqrt((dA/dN+)^2 * N+ + (dA/dN-)^2 * N-)
                if N_plus + N_minus > 0:
                    denom = (N_plus + N_minus)**2
                    dA_dNplus = 2 * N_minus / denom
                    dA_dNminus = -2 * N_plus / denom
                    sigma_A = np.sqrt(dA_dNplus**2 * N_plus + dA_dNminus**2 * N_minus)
                else:
                    sigma_A = 0
                
                # Calculate average DepW/DepA for this bin
                avg_ratio = np.mean(DepW_cut[bin_mask] / DepA_cut[bin_mask])
                
                # Divide by average ratio
                if avg_ratio != 0:
                    asym /= avg_ratio
                    sigma_A /= avg_ratio
                
                asymmetry.append(asym)
                asymmetry_err.append(sigma_A)
            else:
                asymmetry.append(0)
                asymmetry_err.append(0)
        
        asymmetry = np.array(asymmetry)
        asymmetry_err = np.array(asymmetry_err)
        
        # Plot data points
        ax.errorbar(phi_centers, asymmetry, yerr=asymmetry_err, fmt='o', 
                   capsize=5, markersize=6, label=label, color=color)
        
        # Fit to A*sin(phi)
        try:
            # Use weights based on uncertainties (avoid division by zero)
            weights = 1.0 / np.where(asymmetry_err > 0, asymmetry_err, 1.0)
            popt, pcov = curve_fit(sin_fit, phi_centers, asymmetry, sigma=asymmetry_err, 
                                   absolute_sigma=True, p0=[0.1])
            A_fit = popt[0]
            A_err = np.sqrt(pcov[0, 0])
            
            # Plot fit
            phi_fine = np.linspace(0, 2*np.pi, 200)
            ax.plot(phi_fine, sin_fit(phi_fine, A_fit), '-', color=color, 
                   label=f'{label} fit: A={A_fit:.3f}±{A_err:.3f}', linewidth=2)
            
            print(f"{label}: A = {A_fit:.4f} ± {A_err:.4f}")
        except Exception as e:
            print(f"Fit failed for {label}: {e}")

ax.set_xlabel(r"$\phi$", fontsize=14)
ax.set_ylabel(r"$F_{LU}^{sin\phi}/F_{UU}$", fontsize=14)
ax.set_title(r"Beam Spin Asymmetry ($0.1 < x_{B} < 0.15$)", fontsize=16)
ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
ax.grid(True, alpha=0.3)
ax.set_xlim(0, 2*np.pi)
ax.legend(fontsize=10)

# Set x-axis ticks to show multiples of pi
ax.set_xticks([0, np.pi/2, np.pi, 3*np.pi/2, 2*np.pi])
ax.set_xticklabels(['0', r'$\pi/2

print("\nAll plots created successfully!"), r'$\pi

print("\nAll plots created successfully!"), r'$3\pi/2

print("\nAll plots created successfully!"), r'$2\pi

print("\nAll plots created successfully!")])

plt.tight_layout()
plt.savefig(os.path.join(output_dir, "beam_spin_asymmetry.png"), dpi=300)
plt.close()
print(f"Saved: {os.path.join(output_dir, 'beam_spin_asymmetry.png')}")

print("\nAll plots created successfully!")