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
    gaussian = a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    quadratic = b + c * x + d * x ** 2
    return gaussian + quadratic

# Small helpers
def np_branch(tree, name):
    """Read a TTree branch as a 1D NumPy array (float64/int64 as appropriate)."""
    try:
        arr = tree[name].array(library='np')
    except Exception as e:
        raise RuntimeError(f"Missing or unreadable branch '{name}': {e}")
    return np.asarray(arr)
# endif

def parse_hv_label(path):
    """Extract the trailing three underscore-separated tokens before .root."""
    base = os.path.basename(path).replace('.root', '')
    toks = base.split('_')
    if len(toks) >= 3:
        hv = ','.join(toks[-3:])
    else:
        hv = base
    return hv
# endif

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
            if 'PhysicsEvents' not in file:
                print(f"Error: 'PhysicsEvents' not found in {file_path}")
                continue
            tree = file['PhysicsEvents']

            # Read the required branches as NumPy
            Mx2 = np_branch(tree, 'Mx2')
            Q2 = np_branch(tree, 'Q2')
            x = np_branch(tree, 'x')
            e_theta = np_branch(tree, 'e_theta')
            e_p = np_branch(tree, 'e_p')
            p_theta = np_branch(tree, 'p_theta')
            p_p = np_branch(tree, 'p_p')
            detector = np_branch(tree, 'detector')

            # Calculate Mx from Mx2 (taking square root)
            # Guard against negative numerical noise
            Mx2_clipped = np.clip(Mx2, a_min=0.0, a_max=None)
            Mx = np.sqrt(Mx2_clipped)

            # Convert angles from radians to degrees
            e_theta_deg = np.degrees(e_theta)
            p_theta_deg = np.degrees(p_theta)

            # Apply universal cuts: Mx (0.94 to 1.02 GeV), detector == 1, p_theta < 40 degrees
            mask = (Mx >= 0.94) & (Mx <= 1.02) & (detector == 1) & (p_theta_deg < 40.0)

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
# endfor

# Convert to numpy arrays for easier handling
all_Mx = np.asarray(all_Mx, dtype=np.float64)
all_Q2 = np.asarray(all_Q2, dtype=np.float64)
all_x = np.asarray(all_x, dtype=np.float64)
all_e_theta = np.asarray(all_e_theta, dtype=np.float64)
all_e_p = np.asarray(all_e_p, dtype=np.float64)
all_p_theta = np.asarray(all_p_theta, dtype=np.float64)
all_p_p = np.asarray(all_p_p, dtype=np.float64)

print(f"Total events for kinematics plots: {len(all_Mx)}")

# Create the 2x2 kinematics plot
if len(all_Mx) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Top left: Mx histogram as line plot
    counts, bins, _ = axes[0, 0].hist(
        all_Mx, bins=100, range=(0.94, 1.02),
        color='blue', alpha=0.7, histtype='step', linewidth=2
    )
    axes[0, 0].set_xlabel('M_x (GeV)')
    axes[0, 0].set_ylabel('Counts')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title('Missing Mass Distribution')
    axes[0, 0].set_xlim(0.94, 1.02)

    # Top right: Q2 vs x with updated ranges
    if len(all_x) > 0 and len(all_Q2) > 0:
        hist1 = axes[0, 1].hist2d(
            all_x, all_Q2, bins=100, range=[[0.0, 0.6], [1.0, 5.0]],
            cmap='viridis', norm=colors.LogNorm()
        )
        axes[0, 1].set_xlabel('x_B')
        axes[0, 1].set_ylabel('Q^{2} (GeV^{2})')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title('Q^{2} vs x_B')
        axes[0, 1].set_xlim(0.0, 0.6)
        axes[0, 1].set_ylim(1.0, 5.0)
        plt.colorbar(hist1[3], ax=axes[0, 1])
    # endif

    # Bottom left: e_theta vs e_p with updated e_p range
    if len(all_e_p) > 0 and len(all_e_theta) > 0:
        hist2 = axes[1, 0].hist2d(
            all_e_p, all_e_theta, bins=100,
            range=[[2.0, 5.0], [float(np.min(all_e_theta)), float(np.max(all_e_theta))]],
            cmap='viridis', norm=colors.LogNorm()
        )
        axes[1, 0].set_xlabel('e_p (GeV)')
        axes[1, 0].set_ylabel('e_{theta} (degrees)')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].set_title('Electron theta vs Electron Momentum')
        axes[1, 0].set_xlim(2.0, 5.0)
        plt.colorbar(hist2[3], ax=axes[1, 0])
    # endif

    # Bottom right: p_theta vs p_p with angles in degrees
    if len(all_p_p) > 0 and len(all_p_theta) > 0:
        hist3 = axes[1, 1].hist2d(
            all_p_p, all_p_theta, bins=100, cmap='viridis', norm=colors.LogNorm()
        )
        axes[1, 1].set_xlabel('pi^{+}_{p} (GeV)')
        axes[1, 1].set_ylabel('pi^{+}_{theta} (degrees)')
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_title('pi^{+} theta vs pi^{+} Momentum')
        plt.colorbar(hist3[3], ax=axes[1, 1])
    # endif

    # Adjust layout and save
    plt.tight_layout()
    plt.savefig('output/DC_HV_scan/kinematics.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("Kinematics plots saved to output/DC_HV_scan/kinematics.png")
# endif

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
any_plotted = False
line_handles = []
line_labels = []

# Process each file separately for fitting
for i, file_path in enumerate(file_paths):
    try:
        with uproot.open(file_path) as file:
            if 'PhysicsEvents' not in file:
                print(f"Error: 'PhysicsEvents' not found in {file_path}")
                continue

            tree = file['PhysicsEvents']

            # Read the required branches as NumPy
            Mx2 = np_branch(tree, 'Mx2')
            detector = np_branch(tree, 'detector')
            p_theta = np_branch(tree, 'p_theta')

            # Calculate Mx from Mx2 (taking square root, guarded)
            Mx = np.sqrt(np.clip(Mx2, a_min=0.0, a_max=None))

            # Convert p_theta from radians to degrees
            p_theta_deg = np.degrees(p_theta)

            # Apply universal cuts:
            # Mx: 0.94 to 1.02 GeV
            # detector == 1
            # p_theta < 40 degrees
            mask = (Mx >= 0.94) & (Mx <= 1.02) & (detector == 1) & (p_theta_deg < 40.0)

            Mx_cut = Mx[mask]

            if Mx_cut.size == 0:
                print(f"No events in {file_path} after cuts")
                continue

            # Create histogram for fitting
            hist, bin_edges = np.histogram(Mx_cut, bins=80, range=(0.94, 1.02))
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0

            # Initial parameter guess: [amplitude, mean, sigma, const, linear, quadratic]
            # Guess neutron peak around 0.98 GeV with sigma ~ 0.02
            initial_guess = [
                float(np.max(hist)),  # amplitude
                0.98,                 # mean (neutron mass region)
                0.02,                 # sigma
                float(np.min(hist)),  # constant background
                0.0,                  # linear term
                0.0                   # quadratic term
            ]

            # Reasonable bounds for stability
            lower_bounds = [0.0, 0.94, 0.005, -np.inf, -np.inf, -np.inf]
            upper_bounds = [np.inf, 1.02, 0.06,  np.inf,  np.inf,  np.inf]

            # Perform the fit
            popt, pcov = curve_fit(
                gaussian_quadratic,
                bin_centers.astype(np.float64),
                hist.astype(np.float64),
                p0=initial_guess,
                bounds=(lower_bounds, upper_bounds),
                maxfev=20000
            )

            # Extract the fitted parameters
            amplitude, mu, sigma, b, c, d = popt
            # Protect against non-positive definite covariance
            sigma_err = float(np.sqrt(abs(pcov[2, 2]))) if np.isfinite(pcov[2, 2]) else np.nan
            mu_err = float(np.sqrt(abs(pcov[1, 1]))) if np.isfinite(pcov[1, 1]) else np.nan
            amp_err = float(np.sqrt(abs(pcov[0, 0]))) if np.isfinite(pcov[0, 0]) else np.nan

            # Plot the histogram
            plt.hist(Mx_cut, bins=80, range=(0.94, 1.02),
                     histtype='step', linewidth=1.5, color=colors_cycle[i], alpha=0.8)

            # Plot the fit
            x_fit = np.linspace(0.94, 1.02, 200)
            y_fit = gaussian_quadratic(x_fit, *popt)
            (line_obj,) = plt.plot(
                x_fit, y_fit, color=colors_cycle[i], linewidth=2
            )
            hv_label = parse_hv_label(file_path)
            label_text = f"{hv_label}; sigma = {sigma:.4f} +/- {sigma_err:.4f}"
            line_handles.append(line_obj)
            line_labels.append(label_text)
            any_plotted = True

            # Store results for CSV
            fit_results.append({
                'HV_Settings': hv_label,
                'Sigma': float(sigma),
                'Sigma_Error': sigma_err,
                'Mean': float(mu),
                'Mean_Error': mu_err,
                'Amplitude': float(amplitude),
                'Amplitude_Error': amp_err,
                'Events': int(Mx_cut.size)
            })

            print(f"Fit for {file_path}: mu = {mu:.4f}, sigma = {sigma:.4f} +/- {sigma_err:.4f}")

    except Exception as e:
        print(f"Error processing {file_path} for fitting: {e}")
# endfor

# Format the integrated plot
plt.xlabel('M_x (GeV)')
plt.ylabel('Counts')
plt.title('Missing Mass Distributions with Gaussian+Quadratic Fits')
plt.grid(True, alpha=0.3)

if any_plotted:
    plt.yscale('log')  # Log scale to better see all distributions
    plt.legend(line_handles, line_labels, loc='upper right', fontsize=9)
else:
    # Avoid empty-legend warning; annotate for debugging
    plt.text(0.945, 5.0, 'No successful fits (check cuts and input branches)',
             fontsize=10, va='top')
# endif

# Save the integrated plot
plt.tight_layout()
plt.savefig('output/DC_HV_scan/integrated.png', dpi=300, bbox_inches='tight')
plt.close()

print("Integrated plot saved to output/DC_HV_scan/integrated.png")

# Save fit results to CSV
if fit_results:
    df = pd.DataFrame(fit_results)
    df.sort_values(by='HV_Settings', inplace=True)
    df.to_csv('output/DC_HV_scan/fit_results.csv', index=False)
    print("Fit results saved to output/DC_HV_scan/fit_results.csv")

    # Print summary
    print("\nFit Results Summary:")
    print(df[['HV_Settings', 'Sigma', 'Sigma_Error', 'Events']].to_string(index=False))
# endif

print("\nAll tasks completed successfully!")