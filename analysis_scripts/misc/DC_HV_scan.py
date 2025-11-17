import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import colors
from scipy.optimize import curve_fit
import pandas as pd

# Create output directory
os.makedirs('output/DC_HV_scan', exist_ok=True)

# List of root files
file_paths = [
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_10_10.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_10_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_11_10.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_11_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_11_12.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_12_10.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_10_12_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_11_11_11.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_12_13_13.root',
    '/volatile/clas12/thayward/RGK_DC_HV_scan/processed_files/rgk_epi-_9_10_10.root'
]

# Define the fit function: Gaussian + quadratic background
def gaussian_quadratic(x, a, mu, sigma, b, c, d):
    gaussian = a * np.exp(-0.5 * ((x - mu) / sigma)**2)
    quadratic = b + c*x + d*x**2
    return gaussian + quadratic

# =============================================================================
# PART 1: Kinematics Plots (Combined Data)
# =============================================================================

print("Creating kinematics plots...")

# Initialize arrays to store all data for kinematics plots
all_Mx = []
all_Q2 = []
all_x = []
all_e_theta = []
all_e_p = []
all_p_theta = []
all_p_p = []

# Loop over all files and read trees for kinematics
for file_path in file_paths:
    try:
        with uproot.open(file_path) as file:
            tree = file['PhysicsEvents']
            
            # Read the required branches
            Mx2 = tree['Mx2'].array()
            Q2 = tree['Q2'].array()
            x = tree['x'].array()
            e_theta = tree['e_theta'].array()
            e_p = tree['e_p'].array()
            p_theta = tree['p_theta'].array()
            p_p = tree['p_p'].array()
            detector = tree['detector'].array()
            
            # Calculate Mx from Mx2 (taking square root)
            Mx = np.sqrt(Mx2)
            
            # Convert angles from radians to degrees
            e_theta_deg = np.degrees(e_theta)
            p_theta_deg = np.degrees(p_theta)
            
            # Apply universal cuts: Mx (0.94 to 1.02 GeV), detector == 1, p_theta < 40 degrees
            mask = (Mx >= 0.94) & (Mx <= 1.02) & (detector == 1) & (p_theta_deg < 40)
            
            # Append data after applying cut
            all_Mx.extend(Mx[mask])
            all_Q2.extend(Q2[mask])
            all_x.extend(x[mask])
            all_e_theta.extend(e_theta_deg[mask])
            all_e_p.extend(e_p[mask])
            all_p_theta.extend(p_theta_deg[mask])
            all_p_p.extend(p_p[mask])
            
    except Exception as e:
        print(f"Error processing {file_path} for kinematics: {e}")

# Convert to numpy arrays for easier handling
all_Mx = np.array(all_Mx)
all_Q2 = np.array(all_Q2)
all_x = np.array(all_x)
all_e_theta = np.array(all_e_theta)
all_e_p = np.array(all_e_p)
all_p_theta = np.array(all_p_theta)
all_p_p = np.array(all_p_p)

print(f"Total events for kinematics plots: {len(all_Mx)}")

# Create the 2x2 kinematics plot
if len(all_Mx) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Top left: Mx histogram as line plot
    counts, bins, _ = axes[0,0].hist(all_Mx, bins=100, range=(0.94, 1.02), 
                                    color='blue', alpha=0.7, histtype='step', linewidth=2)
    axes[0,0].set_xlabel('M$_x$ (GeV)')
    axes[0,0].set_ylabel('Counts')
    axes[0,0].grid(True, alpha=0.3)
    axes[0,0].set_title('Missing Mass Distribution')
    axes[0,0].set_xlim(0.94, 1.02)

    # Top right: Q2 vs x with updated ranges
    if len(all_x) > 0 and len(all_Q2) > 0:
        hist1 = axes[0,1].hist2d(all_x, all_Q2, bins=100, range=[[0, 0.6], [1, 5]], 
                                cmap='viridis', norm=colors.LogNorm())
        axes[0,1].set_xlabel('x$_B$')
        axes[0,1].set_ylabel('Q$^{2}$ (GeV$^{2}$)')
        axes[0,1].grid(True, alpha=0.3)
        axes[0,1].set_title('Q$^{2}$ vs x$_B$')
        axes[0,1].set_xlim(0, 0.6)
        axes[0,1].set_ylim(1, 5)
        plt.colorbar(hist1[3], ax=axes[0,1])

    # Bottom left: e_theta vs e_p with updated e_p range
    if len(all_e_p) > 0 and len(all_e_theta) > 0:
        hist2 = axes[1,0].hist2d(all_e_p, all_e_theta, bins=100, 
                                range=[[2, 5], [np.min(all_e_theta), np.max(all_e_theta)]],
                                cmap='viridis', norm=colors.LogNorm())
        axes[1,0].set_xlabel('e$_p$ (GeV)')
        axes[1,0].set_ylabel('e$_{\\theta}$ (degrees)')
        axes[1,0].grid(True, alpha=0.3)
        axes[1,0].set_title('Electron $\\theta$ vs Electron Momentum')
        axes[1,0].set_xlim(2, 5)
        plt.colorbar(hist2[3], ax=axes[1,0])

    # Bottom right: p_theta vs p_p with angles in degrees
    if len(all_p_p) > 0 and len(all_p_theta) > 0:
        hist3 = axes[1,1].hist2d(all_p_p, all_p_theta, bins=100, cmap='viridis', norm=colors.LogNorm())
        axes[1,1].set_xlabel(r'$\pi^{+}_{p}$ (GeV)')
        axes[1,1].set_ylabel(r'$\pi^{+}_{\theta}$ (degrees)')
        axes[1,1].grid(True, alpha=0.3)
        axes[1,1].set_title(r'$\pi^{+}$ $\theta$ vs $\pi^{+}$ Momentum')
        plt.colorbar(hist3[3], ax=axes[1,1])

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('output/DC_HV_scan/kinematics.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Kinematics plots saved to output/DC_HV_scan/kinematics.png")

# =============================================================================
# PART 2: Integrated Mx Fits (Individual Files)
# =============================================================================

print("\nCreating integrated Mx fits...")

# Create figure for integrated Mx distributions
plt.figure(figsize=(12, 8))

# Colors for different files
colors_cycle = plt.cm.tab10(np.linspace(0, 1, len(file_paths)))

# Store fit results for CSV
fit_results = []

# Process each file separately for fitting
for i, file_path in enumerate(file_paths):
    try:
        with uproot.open(file_path) as file:
            tree = file['PhysicsEvents']
            
            # Read the required branches
            Mx2 = tree['Mx2'].array()
            detector = tree['detector'].array()
            p_theta = tree['p_theta'].array()
            
            # Calculate Mx from Mx2 (taking square root)
            Mx = np.sqrt(Mx2)
            
            # Convert p_theta from radians to degrees
            p_theta_deg = np.degrees(p_theta)
            
            # Apply universal cuts:
            # Mx: 0.94 to 1.02 GeV
            # detector == 1
            # p_theta < 40 degrees
            mask = (Mx >= 0.94) & (Mx <= 1.02) & (detector == 1) & (p_theta_deg < 40)
            
            Mx_cut = Mx[mask]
            
            if len(Mx_cut) == 0:
                print(f"No events in {file_path} after cuts")
                continue
            
            # Create histogram for fitting
            hist, bin_edges = np.histogram(Mx_cut, bins=80, range=(0.94, 1.02))
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            # Initial parameter guess: [amplitude, mean, sigma, const, linear, quadratic]
            # Guess neutron peak around 0.98 GeV with sigma ~ 0.02
            initial_guess = [
                np.max(hist),  # amplitude
                0.98,          # mean (neutron mass)
                0.02,          # sigma
                np.min(hist),  # constant background
                0,             # linear term
                0              # quadratic term
            ]
            
            # Perform the fit
            try:
                popt, pcov = curve_fit(gaussian_quadratic, bin_centers, hist, p0=initial_guess, maxfev=5000)
                
                # Extract the fitted parameters
                amplitude, mu, sigma, b, c, d = popt
                sigma_err = np.sqrt(pcov[2, 2])  # Error on sigma
                
                # Plot the histogram
                plt.hist(Mx_cut, bins=80, range=(0.94, 1.02), 
                        histtype='step', linewidth=1.5, color=colors_cycle[i], alpha=0.8)
                
                # Plot the fit
                x_fit = np.linspace(0.94, 1.02, 200)
                y_fit = gaussian_quadratic(x_fit, *popt)
                plt.plot(x_fit, y_fit, color=colors_cycle[i], linewidth=2, 
                        label=f"{file_path.split('_')[-3]},{file_path.split('_')[-2]},{file_path.split('_')[-1].split('.')[0]}; σ = {sigma:.4f}±{sigma_err:.4f}")
                
                # Store results for CSV
                hv_settings = file_path.split('_')[-3] + ',' + file_path.split('_')[-2] + ',' + file_path.split('_')[-1].split('.')[0]
                fit_results.append({
                    'HV_Settings': hv_settings,
                    'Sigma': sigma,
                    'Sigma_Error': sigma_err,
                    'Mean': mu,
                    'Mean_Error': np.sqrt(pcov[1, 1]),
                    'Amplitude': amplitude,
                    'Amplitude_Error': np.sqrt(pcov[0, 0]),
                    'Events': len(Mx_cut)
                })
                
                print(f"Fit for {file_path}: μ = {mu:.4f}, σ = {sigma:.4f} ± {sigma_err:.4f}")
                
            except Exception as e:
                print(f"Fit failed for {file_path}: {e}")
            
    except Exception as e:
        print(f"Error processing {file_path} for fitting: {e}")

# Format the integrated plot
plt.xlabel('M$_x$ (GeV)')
plt.ylabel('Counts')
plt.title('Missing Mass Distributions with Gaussian+Quadratic Fits')
plt.legend(loc='upper right', fontsize=9)
plt.grid(True, alpha=0.3)
plt.yscale('log')  # Log scale to better see all distributions

# Save the integrated plot
plt.tight_layout()
plt.savefig('output/DC_HV_scan/integrated.png', dpi=300, bbox_inches='tight')
plt.close()

print("Integrated plot saved to output/DC_HV_scan/integrated.png")

# Save fit results to CSV
if fit_results:
    df = pd.DataFrame(fit_results)
    df.to_csv('output/DC_HV_scan/fit_results.csv', index=False)
    print("Fit results saved to output/DC_HV_scan/fit_results.csv")
    
    # Print summary
    print("\nFit Results Summary:")
    print(df[['HV_Settings', 'Sigma', 'Sigma_Error', 'Events']].to_string(index=False))

print("\nAll tasks completed successfully!")