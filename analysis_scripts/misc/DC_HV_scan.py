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

# Create figure for integrated Mx distributions
plt.figure(figsize=(12, 8))

# Colors for different files
colors_cycle = plt.cm.tab10(np.linspace(0, 1, len(file_paths)))

# Store fit results for CSV
fit_results = []

# Process each file separately
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
        print(f"Error processing {file_path}: {e}")

# Format the plot
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