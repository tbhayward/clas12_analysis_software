import uproot
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib import colors
from scipy.optimize import curve_fit
import pandas as pd

# -----------------------------
# Config / Output
# -----------------------------
OUT_DIR = 'output/DC_HV_scan'
os.makedirs(OUT_DIR, exist_ok=True)

# List of root files (we'll reorder to put 9,10,10 first)
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

# -----------------------------
# Helpers
# -----------------------------
def np_branch(tree, name):
    """Read a TTree branch as a 1D NumPy array."""
    try:
        arr = tree[name].array(library='np')
    except Exception as e:
        raise RuntimeError(f"Missing or unreadable branch '{name}': {e}")
    return np.asarray(arr)
# endif

def parse_hv_label(path):
    """Extract the trailing three underscore-separated tokens before .root, joined by commas."""
    base = os.path.basename(path).replace('.root', '')
    toks = base.split('_')
    if len(toks) >= 3:
        hv = ','.join(toks[-3:])
    else:
        hv = base
    return hv
# endif

def hv_tuple(hv_str):
    """Turn 'a,b,c' into (a,b,c) of ints when possible; fallback keeps sort stable."""
    try:
        a, b, c = hv_str.split(',')
        return (int(a), int(b), int(c))
    except Exception:
        return (999, 999, 999)
# endif

def hv_sort_key(path):
    hv = parse_hv_label(path)
    # Put 9,10,10 first; then sort numerically
    return (0 if hv == '9,10,10' else 1, hv_tuple(hv), hv)
# endif

# Reorder files so 9,10,10 is first
file_paths = sorted(file_paths, key=hv_sort_key)

# Consistent color assignment per HV across all plots
HV_LABELS = [parse_hv_label(p) for p in file_paths]
color_map = dict(zip(HV_LABELS, plt.cm.tab10(np.linspace(0, 1, len(HV_LABELS)))))

# Fit function: Gaussian + linear background
def gaussian_linear(x, a, mu, sigma, b, c):
    gaussian = a * np.exp(-0.5 * ((x - mu) / sigma) ** 2)
    linear   = b + c * x
    return gaussian + linear
# endif

# Common histogram settings for Mx
MX_MIN, MX_MAX = 0.95, 1.00
MX_NBINS = 80

# Universal selection toggles
REQ_DETECTOR_1 = True
REQ_PTHETA_LT_40 = True

# -----------------------------
# PART 1: Kinematics Plots (Combined Data)
# -----------------------------
print("Creating kinematics plots...")

all_Mx, all_Q2, all_x = [], [], []
all_e_theta, all_e_p = [], []
all_p_theta, all_p_p = [], []

for file_path in file_paths:
    try:
        with uproot.open(file_path) as f:
            if 'PhysicsEvents' not in f:
                print(f"Error: 'PhysicsEvents' not found in {file_path}")
                continue
            tree = f['PhysicsEvents']

            # Read required branches as NumPy
            Mx2      = np_branch(tree, 'Mx2')
            Q2       = np_branch(tree, 'Q2')
            x        = np_branch(tree, 'x')
            e_theta  = np_branch(tree, 'e_theta')   # radians
            e_p      = np_branch(tree, 'e_p')       # GeV
            p_theta  = np_branch(tree, 'p_theta')   # radians
            p_p      = np_branch(tree, 'p_p')       # GeV
            detector = np_branch(tree, 'detector')  # int

            # Derived
            Mx = np.sqrt(np.clip(Mx2, a_min=0.0, a_max=None))
            e_theta_deg = np.degrees(e_theta)
            p_theta_deg = np.degrees(p_theta)

            # Universal cuts for the combined kinematics display
            mask = (Mx > MX_MIN) & (Mx < MX_MAX)
            if REQ_DETECTOR_1:
                mask &= (detector == 1)
            # endif
            if REQ_PTHETA_LT_40:
                mask &= (p_theta_deg < 40.0)
            # endif

            all_Mx.extend(Mx[mask])
            all_Q2.extend(Q2[mask])
            all_x.extend(x[mask])
            all_e_theta.extend(e_theta_deg[mask])
            all_e_p.extend(e_p[mask])
            all_p_theta.extend(p_theta_deg[mask])
            all_p_p.extend(p_p[mask])
        # endwith
    except Exception as e:
        print(f"Error processing {file_path} for kinematics: {e}")
# endfor

all_Mx     = np.asarray(all_Mx, dtype=np.float64)
all_Q2     = np.asarray(all_Q2, dtype=np.float64)
all_x      = np.asarray(all_x, dtype=np.float64)
all_e_theta= np.asarray(all_e_theta, dtype=np.float64)
all_e_p    = np.asarray(all_e_p, dtype=np.float64)
all_p_theta= np.asarray(all_p_theta, dtype=np.float64)
all_p_p    = np.asarray(all_p_p, dtype=np.float64)

print(f"Total events for kinematics plots: {len(all_Mx)}")

if len(all_Mx) > 0:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # Top-left: Mx distribution
    axes[0, 0].hist(all_Mx, bins=100, range=(MX_MIN, MX_MAX),
                    histtype='step', linewidth=2)
    axes[0, 0].set_xlabel(r'$M_{x}$ (GeV)')
    axes[0, 0].set_ylabel('Counts')
    axes[0, 0].grid(True, alpha=0.3)
    axes[0, 0].set_title(r'Missing Mass Distribution')
    axes[0, 0].set_xlim(MX_MIN, MX_MAX)

    # Top-right: Q2 vs x_B
    if (len(all_x) > 0) and (len(all_Q2) > 0):
        h = axes[0, 1].hist2d(all_x, all_Q2, bins=100,
                              range=[[0.0, 0.6], [1.0, 5.0]],
                              cmap='viridis', norm=colors.LogNorm())
        axes[0, 1].set_xlabel(r'$x_{B}$')
        axes[0, 1].set_ylabel(r'$Q^{2}$ (GeV$^{2}$)')
        axes[0, 1].grid(True, alpha=0.3)
        axes[0, 1].set_title(r'$Q^{2}$ vs $x_{B}$')
        axes[0, 1].set_xlim(0.0, 0.6)
        axes[0, 1].set_ylim(1.0, 5.0)
        plt.colorbar(h[3], ax=axes[0, 1])
    # endif

    # Bottom-left: e_theta vs e_p
    if (len(all_e_p) > 0) and (len(all_e_theta) > 0):
        h2 = axes[1, 0].hist2d(all_e_p, all_e_theta, bins=100,
                               range=[[2.0, 5.0],
                                      [float(np.min(all_e_theta)),
                                       float(np.max(all_e_theta))]],
                               cmap='viridis', norm=colors.LogNorm())
        axes[1, 0].set_xlabel(r'$e_{p}$ (GeV)')
        axes[1, 0].set_ylabel(r'$e_{\theta}$ (^\circ)')
        axes[1, 0].grid(True, alpha=0.3)
        axes[1, 0].set_title(r'Electron $\theta$ vs Electron Momentum')
        axes[1, 0].set_xlim(2.0, 5.0)
        plt.colorbar(h2[3], ax=axes[1, 0])
    # endif

    # Bottom-right: pi+ theta vs pi+ momentum
    if (len(all_p_p) > 0) and (len(all_p_theta) > 0):
        h3 = axes[1, 1].hist2d(all_p_p, all_p_theta, bins=100,
                               cmap='viridis', norm=colors.LogNorm())
        axes[1, 1].set_xlabel(r'$\pi^{+}_{p}$ (GeV)')
        axes[1, 1].set_ylabel(r'$\pi^{+}_{\theta}$ (^\circ)')
        axes[1, 1].grid(True, alpha=0.3)
        axes[1, 1].set_title(r'$\pi^{+}$ $\theta$ vs $\pi^{+}$ Momentum')
        plt.colorbar(h3[3], ax=axes[1, 1])
    # endif

    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, 'kinematics.png'),
                dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Kinematics plots saved to {OUT_DIR}/kinematics.png")
# endif

# -----------------------------
# PART 2: Integrated Mx Fits (Individual Files)
# -----------------------------
print("\nCreating integrated Mx fits...")

plt.figure(figsize=(12, 8))
fit_results = []
any_plotted = False
handles, labels = [], []

for path in file_paths:
    try:
        with uproot.open(path) as f:
            if 'PhysicsEvents' not in f:
                print(f"Error: 'PhysicsEvents' not found in {path}")
                continue
            t = f['PhysicsEvents']

            Mx2      = np_branch(t, 'Mx2')
            detector = np_branch(t, 'detector')
            p_theta  = np_branch(t, 'p_theta')   # radians

            Mx = np.sqrt(np.clip(Mx2, a_min=0.0, a_max=None))
            p_theta_deg = np.degrees(p_theta)

            mask = (Mx > MX_MIN) & (Mx < MX_MAX)
            if REQ_DETECTOR_1:
                mask &= (detector == 1)
            # endif
            if REQ_PTHETA_LT_40:
                mask &= (p_theta_deg < 40.0)
            # endif

            data = Mx[mask]
            if data.size == 0:
                print(f"No events in {path} after cuts")
                continue

            hist, edges = np.histogram(data, bins=MX_NBINS, range=(MX_MIN, MX_MAX))
            centers = 0.5 * (edges[:-1] + edges[1:])

            # Gaussian + linear initial guesses/bounds
            p0 = [float(np.max(hist)), 0.975, 0.015, float(np.min(hist)), 0.0]
            lb = [0.0, MX_MIN, 0.003, -np.inf, -np.inf]
            ub = [np.inf, MX_MAX, 0.040,  np.inf,  np.inf]

            popt, pcov = curve_fit(
                gaussian_linear,
                centers.astype(np.float64),
                hist.astype(np.float64),
                p0=p0, bounds=(lb, ub), maxfev=20000
            )

            a, mu, sigma, b, c = popt
            sigma_err = float(np.sqrt(abs(pcov[2, 2]))) if np.isfinite(pcov[2, 2]) else np.nan
            mu_err    = float(np.sqrt(abs(pcov[1, 1]))) if np.isfinite(pcov[1, 1]) else np.nan
            a_err     = float(np.sqrt(abs(pcov[0, 0]))) if np.isfinite(pcov[0, 0]) else np.nan

            # Draw histogram (outline) and fit line
            hv = parse_hv_label(path)
            color = color_map[hv]

            plt.hist(data, bins=MX_NBINS, range=(MX_MIN, MX_MAX),
                     histtype='step', linewidth=1.2, color=color, alpha=0.9)

            x_fit = np.linspace(MX_MIN, MX_MAX, 250)
            y_fit = gaussian_linear(x_fit, *popt)
            (line,) = plt.plot(x_fit, y_fit, color=color, linewidth=2.0)

            # Legend line label with #sigma
            label_text = rf"{hv}; $\#\sigma = {sigma:.4f} \pm {sigma_err:.4f}$"
            handles.append(line)
            labels.append(label_text)
            any_plotted = True

            fit_results.append({
                'HV_Settings': hv,
                'Sigma': float(sigma),
                'Sigma_Error': sigma_err,
                'Mean': float(mu),
                'Mean_Error': mu_err,
                'Amplitude': float(a),
                'Amplitude_Error': a_err,
                'Events': int(data.size)
            })

            print(f"Fit for {path}: mu = {mu:.4f}, sigma = {sigma:.4f} +/- {sigma_err:.4f}")

    except Exception as e:
        print(f"Error processing {path} for fitting: {e}")
# endfor

plt.xlabel(r'$M_{x}$ (GeV)')
plt.ylabel('Counts')
plt.title(r'Missing Mass Distributions with Gaussian + Linear Fits')
plt.grid(True, alpha=0.3)
plt.xlim(MX_MIN, MX_MAX)

if any_plotted:
    plt.legend(handles, labels, loc='upper right', fontsize=9)
# endif

plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, 'integrated.png'), dpi=300, bbox_inches='tight')
plt.close()
print(f"Integrated plot saved to {OUT_DIR}/integrated.png")

# Write CSV (sorted with 9,10,10 first; rounded to 4 decimals)
if fit_results:
    df = pd.DataFrame(fit_results)

    # Round numeric columns to 4 decimals
    num_cols = ['Sigma', 'Sigma_Error', 'Mean', 'Mean_Error', 'Amplitude', 'Amplitude_Error']
    df[num_cols] = df[num_cols].round(4)

    # Sort: 9,10,10 first, then numeric hv ordering
    def df_hv_key(hv):
        return (0 if hv == '9,10,10' else 1, hv_tuple(hv), hv)
    # endif
    df = df.sort_values(by='HV_Settings', key=lambda col: col.map(df_hv_key))

    csv_path = os.path.join(OUT_DIR, 'integrated_fit_results.csv')
    df.to_csv(csv_path, index=False)
    print(f"Fit results saved to {csv_path}")

    print("\nFit Results Summary:")
    print(df[['HV_Settings', 'Sigma', 'Sigma_Error', 'Events']].to_string(index=False))
# endif

# -----------------------------
# PART 3: 2x3 canvas with e_theta (top) and p_theta (bottom), with fits
# -----------------------------
print("\nCreating 2x3 binned overlays with fits...")

e_bins = [(12.0, 15.0), (15.0, 20.0), (20.0, 35.0)]  # degrees
p_bins = [(9.0, 16.0), (16.0, 24.0), (24.0, 40.0)]   # degrees

fig, axes = plt.subplots(2, 3, figsize=(18, 9), sharey='row')

binned_rows = [
    ('e_theta', e_bins, 0),
    ('p_theta', p_bins, 1),
]

binned_sigma_rows = []  # will accumulate rows for the binned CSV

for cat, bins, row in binned_rows:
    for col, (lo, hi) in enumerate(bins):
        ax = axes[row, col]
        handles_loc, labels_loc = [], []
        any_here = False

        for path in file_paths:
            try:
                with uproot.open(path) as f:
                    if 'PhysicsEvents' not in f:
                        continue
                    t = f['PhysicsEvents']

                    Mx2      = np_branch(t, 'Mx2')
                    e_theta  = np_branch(t, 'e_theta')   # radians
                    p_theta  = np_branch(t, 'p_theta')   # radians
                    detector = np_branch(t, 'detector')

                    Mx   = np.sqrt(np.clip(Mx2, a_min=0.0, a_max=None))
                    e_deg = np.degrees(e_theta)
                    p_deg = np.degrees(p_theta)

                    mask = (Mx > MX_MIN) & (Mx < MX_MAX)
                    if REQ_DETECTOR_1:
                        mask &= (detector == 1)
                    # endif
                    if REQ_PTHETA_LT_40:
                        mask &= (p_deg < 40.0)
                    # endif

                    if cat == 'e_theta':
                        mask &= (e_deg >= lo) & (e_deg < hi)
                    else:
                        mask &= (p_deg >= lo) & (p_deg < hi)
                    # endif

                    data = Mx[mask]
                    if data.size == 0:
                        continue

                    hv = parse_hv_label(path)
                    color = color_map[hv]

                    # Draw histogram
                    ax.hist(data, bins=MX_NBINS, range=(MX_MIN, MX_MAX),
                            histtype='step', linewidth=1.1, color=color, alpha=0.95)

                    # Fit (Gaussian + linear)
                    hist, edges = np.histogram(data, bins=MX_NBINS, range=(MX_MIN, MX_MAX))
                    centers = 0.5 * (edges[:-1] + edges[1:])

                    p0 = [float(np.max(hist)), 0.975, 0.015, float(np.min(hist)), 0.0]
                    lb = [0.0, MX_MIN, 0.003, -np.inf, -np.inf]
                    ub = [np.inf, MX_MAX, 0.040,  np.inf,  np.inf]

                    try:
                        popt, pcov = curve_fit(
                            gaussian_linear,
                            centers.astype(np.float64),
                            hist.astype(np.float64),
                            p0=p0, bounds=(lb, ub), maxfev=20000
                        )
                        a, mu, sigma, b, c = popt
                        sigma_err = float(np.sqrt(abs(pcov[2, 2]))) if np.isfinite(pcov[2, 2]) else np.nan

                        x_fit = np.linspace(MX_MIN, MX_MAX, 250)
                        y_fit = gaussian_linear(x_fit, *popt)
                        (line,) = ax.plot(x_fit, y_fit, color=color, linewidth=2.0)

                        # Legend label with #sigma
                        lbl = rf"{hv}; $\#\sigma = {sigma:.4f} \pm {sigma_err:.4f}$"
                        handles_loc.append(line)
                        labels_loc.append(lbl)
                        any_here = True

                        # Save to binned CSV rows (only sigma)
                        bin_label = f"[{int(lo)},{int(hi)})"
                        binned_sigma_rows.append({
                            'Category': 'e_theta' if cat == 'e_theta' else 'p_theta',
                            'Bin': bin_label,
                            'HV_Settings': hv,
                            'Sigma': round(float(sigma), 4)
                        })
                    except Exception:
                        # If fit fails, skip adding curve/legend for this HV/bin
                        continue
                    # endif
                # endwith
            except Exception:
                continue
        # endfor

        ax.set_xlabel(r'$M_{x}$ (GeV)')
        if col == 0:
            ax.set_ylabel('Counts')
        # endif
        ax.grid(True, alpha=0.3)
        ax.set_xlim(MX_MIN, MX_MAX)

        if cat == 'e_theta':
            ax.set_title(rf'$e_{{\theta}} \in [{lo:.0f}, {hi:.0f})^\circ$')
        else:
            ax.set_title(rf'$p_{{\theta}} \in [{lo:.0f}, {hi:.0f})^\circ$')
        # endif

        if any_here:
            ax.legend(handles_loc, labels_loc, loc='upper right', fontsize=8)
        # endif
    # endfor
# endfor

plt.tight_layout()
binned_png = os.path.join(OUT_DIR, 'mx_binned_theta_2x3.png')
plt.savefig(binned_png, dpi=300, bbox_inches='tight')
plt.close()
print(f"Saved {binned_png}")

# Write binned CSV with only sigma
if binned_sigma_rows:
    bdf = pd.DataFrame(binned_sigma_rows)
    # Sort: by Category (e_theta first), then bin order, then HV with 9,10,10 first
    cat_order = {'e_theta': 0, 'p_theta': 1}
    def bin_key(s):
        # s like "[12,15)"; extract numbers
        s = s.strip('[]')
        lo, hi = s.replace(')', '').split(',')
        return (int(lo), int(hi))
    # endif
    def hv_key(hv):
        return (0 if hv == '9,10,10' else 1, hv_tuple(hv), hv)

    bdf = bdf.sort_values(
        by=['Category', 'Bin', 'HV_Settings'],
        key=lambda col: (
            col.map(cat_order) if col.name == 'Category'
            else col.map(lambda x: bin_key(x)) if col.name == 'Bin'
            else col.map(hv_key)
        )
    )
    binned_csv = os.path.join(OUT_DIR, 'binned_fit_sigmas.csv')
    bdf.to_csv(binned_csv, index=False)
    print(f"Binned sigmas saved to {binned_csv}")
# endif

print("\nAll tasks completed successfully!")