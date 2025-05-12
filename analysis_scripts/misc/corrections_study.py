import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.optimize import curve_fit

# Rho0 mass and squared mass in GeV²
m_rho0 = 0.77526
m_rho0_sq = m_rho0**2

def generate_phase_space_plots(channel, correction, plot_type, parent_dir, output_dir):
    """
    Generate a 1x3 canvas of 2D histograms for proton phase space analysis.
    """
    # Configuration based on channel and plot type
    if channel == 'ep':
        file_prefix = 'resIncl'
        branch_config = {
            'theta_phi': ('p_theta', 'p_phi'),
            'theta_p': ('p_theta', 'p_p')
        }
    elif channel == 'eppi+pi-':
        file_prefix = 'nSidis'
        branch_config = {
            'theta_phi': ('p1_theta', 'p1_phi'),
            'theta_p': ('p1_theta', 'p1_p')
        }
    else:
        raise ValueError(f"Invalid channel: {channel}")

    try:
        y_branch, x_branch = branch_config[plot_type]
    except KeyError:
        raise ValueError(f"Invalid plot_type: {plot_type}. Must be 'theta_phi' or 'theta_p'")

    # Plot parameters
    plot_config = {
        'theta_phi': {
            'x_range': [0, 360],
            'y_range': [0, 140],
            'x_label': r'$\phi$ [deg]',
            'y_label': r'$\theta$ [deg]'
        },
        'theta_p': {
            'x_range': [0, 7],
            'y_range': [0, 140],
            'x_label': 'p [GeV]',
            'y_label': r'$\theta$ [deg]'
        }
    }

    run_periods = {
        'fa18_inb': 'RGA Fa18 Inb',
        'fa18_out': 'RGA Fa18 Out',
        'sp19_inb': 'RGA Sp19 Inb'
    }

    # Create figure
    fig, axs = plt.subplots(1, 3, figsize=(18, 6))
    plt.subplots_adjust(wspace=0.3, left=0.05, right=0.95)

    for idx, (run_key, run_label) in enumerate(run_periods.items()):
        filename = f"{file_prefix}_{run_key}_{channel}_{correction}.root"
        filepath = os.path.join(parent_dir, filename)
        
        if not os.path.exists(filepath):
            print(f"Warning: File not found - {filepath}")
            continue

        try:
            with uproot.open(filepath) as f:
                tree = f['PhysicsEvents']
                theta = np.degrees(tree[y_branch].array(library='np'))
                x_data = tree[x_branch].array(library='np')
                
                if plot_type == 'theta_phi':
                    x_data = np.degrees(x_data) % 360

        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            continue

        ax = axs[idx]
        h = ax.hist2d(x_data, theta, bins=100,
                     range=[plot_config[plot_type]['x_range'], 
                            plot_config[plot_type]['y_range']],
                     cmap='viridis', norm=LogNorm(), cmin=1)

        ax.set_xlabel(plot_config[plot_type]['x_label'], fontsize=12)
        ax.set_ylabel(plot_config[plot_type]['y_label'], fontsize=12)
        ax.set_title(run_label, fontsize=14)
        fig.colorbar(h[3], ax=ax).set_label('Counts', fontsize=10)

    output_file = os.path.join(output_dir, f"{channel}_{correction}_{plot_type}_phase_space.pdf")
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()
    print(f"Saved: {output_file}")

def plot_mx2_comparison(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) spectra for different corrections
    """
    detectors = {
        1: {
            'name': 'Forward',
            'theta_bins': [0, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 80]
        },
        2: {
            'name': 'Central',
            'theta_bins': [0, 36, 39, 42, 45, 48, 51, 54, 57, 180]
        }
    }

    # build theta_labels for each detector, skipping the first (0–edge) bin
    for det_config in detectors.values():
        bins = det_config['theta_bins']
        det_config['theta_labels'] = ['All θ'] + [
            f"{bins[i]}-{bins[i+1]}" for i in range(1, len(bins)-1)
        ]
    #endfor

    corrections = ['noCorrections', 'timothy', 'krishna', 'mariana']
    corr_labels = {
        'noCorrections': 'No Corrections',
        'timothy': "Timothy's",
        'krishna': "Krishna's",
        'mariana': "Mariana's"
    }
    colors = ['black', 'red', 'green', 'blue']
    line_styles = ['-', '--', ':', '-.']
    run_periods = ['fa18_inb', 'fa18_out', 'sp19_inb']

    mx2_bins = np.linspace(-0.2, 0.2, 100)

    for run in run_periods:
        for det_num, det_config in detectors.items():
            n_total_plots = len(det_config['theta_labels'])
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots - 1) / n_cols)) + 1

            fig = plt.figure(figsize=(20, 5 * n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols, figure=fig)

            all_data = {}
            for corr in corrections:
                filename = f"resIncl_{run}_ep_{corr}.root"
                filepath = os.path.join(parent_dir, filename)

                try:
                    with uproot.open(filepath) as f:
                        tree = f['PhysicsEvents']
                        data = tree.arrays(['Mx2', 'p_theta', 'detector'], library='np')
                        mask = (data['detector'] == det_num)

                        if mask.sum() > 0:
                            all_data[corr] = {
                                'Mx2': data['Mx2'][mask],
                                'theta': np.degrees(data['p_theta'][mask])
                            }
                        #endif
                except Exception as e:
                    print(f"Error loading {filepath}: {e}")
                #endtry
            #endfor

            # Integrated spectrum plot
            ax_int = fig.add_subplot(gs[0, :])
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr in all_data:
                    ax_int.hist(
                        all_data[corr]['Mx2'], bins=mx2_bins,
                        histtype='step', color=color, linestyle=ls,
                        label=corr_labels[corr]
                    )
                #endif
            #endfor

            ax_int.set(
                xlabel=r'$M_{x}^{2}$ (GeV²)', ylabel='Counts',
                xlim=(-0.2, 0.2),
                title=f"{det_config['name']} Detector - {run}"
            )
            ax_int.legend()
            ax_int.grid(True, alpha=0.3)

            # Theta-binned spectra
            for idx in range(1, n_total_plots):
                row = (idx - 1) // n_cols + 1
                col = (idx - 1) % n_cols
                ax = fig.add_subplot(gs[row, col])

                theta_min = det_config['theta_bins'][idx]
                theta_max = det_config['theta_bins'][idx + 1]

                artists = []
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr in all_data:
                        mask = (
                            (all_data[corr]['theta'] >= theta_min) &
                            (all_data[corr]['theta'] < theta_max)
                        )
                        mx2_data = all_data[corr]['Mx2'][mask]
                        if len(mx2_data) > 0:
                            _, _, patches = ax.hist(
                                mx2_data, bins=mx2_bins,
                                histtype='step', color=color, linestyle=ls,
                                label=corr_labels[corr]
                            )
                            artists.append(patches[0])
                        #endif
                    #endif
                #endfor

                ax.set(
                    xlabel=r'$M_{x}^{2}$ (GeV²)',
                    ylabel='Counts',
                    xlim=(-0.2, 0.2),
                    title=f'θ: {det_config["theta_labels"][idx]}°'
                )
                if artists:
                    ax.legend(handles=artists, fontsize=8)
                #endif
            #endfor

            output_file = os.path.join(
                output_dir,
                f'Mx2_comparison_{run}_{det_config["name"]}.pdf'
            )
            plt.savefig(output_file, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved: {output_file}")
        #endfor
    #endfor

def gauss_poly(x, A, mu, sigma, p2, p1, p0):
    """Gaussian plus quadratic polynomial."""
    gauss = A * np.exp(-(x - mu)**2 / (2 * sigma**2))
    poly = p2*x**2 + p1*x + p0
    return gauss + poly

def plot_three_particles(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) for eppi+pi- channel
    with three particle final state, fits a Gaussian+quadratic to each distribution,
    overlays the fit, and adds μ and σ to the legend.
    """
    detectors = {
        1: {
            'name': 'Forward',
            'theta_bins': [0, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 80],
            'theta_labels': ['All θ'] + [f'{b}-{b+3}' for b in [8,11,14,17,20,23,26,29,32,35,38]] + ['38-80']
        },
        2: {
            'name': 'Central',
            'theta_bins': [0, 36, 39, 42, 45, 48, 51, 54, 57, 180],
            'theta_labels': ['All θ'] + [f'{b}-{b+3}' for b in [36,39,42,45,48,51,54]] + ['54-180']
        }
    }

    corrections = ['noCorrections', 'timothy', 'krishna', 'mariana']
    corr_labels = {
        'noCorrections': 'No Corrections',
        'timothy': "Timothy's",
        'krishna': "Krishna's",
        'mariana': "Mariana's"
    }
    colors = ['black', 'red', 'green', 'blue']
    line_styles = ['-', '--', ':', '-.']
    run_periods = ['fa18_inb', 'fa18_out', 'sp19_inb']
    
    # Halved number of bins, range [0.4, 0.9]
    mx2_bins = np.linspace(0.4, 0.9, 50)

    for run in run_periods:
        print(f"\n{'#'*80}")
        print(f"Processing run period: {run}")
        print(f"{'#'*80}")
        
        for det_num, det_config in detectors.items():
            print(f"\n{'='*60}")
            print(f"Analyzing {det_config['name']} Detector (ID: {det_num})")
            print(f"Theta bins: {det_config['theta_bins']}")
            print(f"Theta labels: {det_config['theta_labels']}")
            print(f"{'='*60}\n")

            bins = det_config['theta_bins']
            labels = det_config['theta_labels']
            if len(bins) != len(labels) + 1:
                print("! Critical Error: Bin/Label mismatch!")
                continue

            n_total_plots = len(labels)
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots-1)/n_cols)) + 1
            
            all_data = {}
            
            # Load data
            for corr in corrections:
                filename = f"nSidis_{run}_{corr}.root"
                filepath = os.path.join(parent_dir, filename)
                if not os.path.exists(filepath):
                    print(f"! File not found: {filepath}")
                    continue
                try:
                    with uproot.open(filepath) as f:
                        if 'PhysicsEvents' not in f:
                            print("! Missing PhysicsEvents tree")
                            continue
                        tree = f['PhysicsEvents']
                        required = ['Mx2_1', 'p1_theta', 'detector1']
                        if any(b not in tree for b in required):
                            print(f"! Missing branches in {filename}")
                            continue
                        data = tree.arrays(required, library='np')
                        mask = (data['detector1'] == det_num)
                        if np.sum(mask) == 0:
                            continue
                        all_data[corr] = {
                            'Mx2_1': data['Mx2_1'][mask],
                            'theta': np.degrees(data['p1_theta'][mask])
                        }
                except Exception as e:
                    print(f"! Error loading {filename}: {e}")
                    continue

            if not all_data:
                print("! No data for any corrections, skipping")
                continue

            # Plotting
            fig = plt.figure(figsize=(20, 5*n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols, figure=fig)

            # Integrated
            ax_int = fig.add_subplot(gs[0, :])
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr not in all_data: 
                    continue
                vals = all_data[corr]['Mx2_1']
                counts, edges = np.histogram(vals, bins=mx2_bins)
                centers = (edges[:-1] + edges[1:]) / 2
                try:
                    popt, _ = curve_fit(
                        gauss_poly, centers, counts,
                        p0=[counts.max(), m_rho0_sq, 0.1, 0, 0, 0]
                    )
                    mu, sigma = popt[1], abs(popt[2])
                except Exception as e:
                    print(f"! Fit failed for {corr}: {e}")
                    mu, sigma = np.nan, np.nan

                label = f"{corr_labels[corr]} (μ={mu:.3f}, σ={sigma:.3f})"
                _, _, ph = ax_int.hist(
                    vals, bins=mx2_bins,
                    histtype='step',
                    color=color,
                    linestyle=ls,
                    label=label
                )
                # overlay full fit
                x_fit = np.linspace(0.4, 0.9, 200)
                y_fit = gauss_poly(x_fit, *popt)
                ax_int.plot(x_fit, y_fit, color=color, linestyle=ls)

            ax_int.set(
                xlim=(0.4, 0.9),
                xlabel=r'$M_{x (ep)}^{2}$ (GeV²)',
                ylabel='Counts',
                title=f"{det_config['name']} Detector - {run}"
            )
            ax_int.legend()
            ax_int.grid(True, alpha=0.3)

            # Theta-binned
            for idx in range(1, n_total_plots):
                ax = fig.add_subplot(gs[(idx-1)//n_cols + 1, (idx-1)%n_cols])
                tmin, tmax = bins[idx], bins[idx+1]
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr not in all_data:
                        continue
                    mask = ((all_data[corr]['theta'] >= tmin) &
                            (all_data[corr]['theta'] < tmax))
                    slice_vals = all_data[corr]['Mx2_1'][mask]
                    if len(slice_vals) == 0:
                        continue
                    counts, edges = np.histogram(slice_vals, bins=mx2_bins)
                    centers = (edges[:-1] + edges[1:]) / 2
                    try:
                        popt, _ = curve_fit(
                            gauss_poly, centers, counts,
                            p0=[counts.max(), m_rho0_sq, 0.1, 0, 0, 0]
                        )
                        mu, sigma = popt[1], abs(popt[2])
                    except Exception as e:
                        print(f"! Fit failed for {corr} in θ bin {labels[idx]}: {e}")
                        mu, sigma = np.nan, np.nan

                    label = f"{corr_labels[corr]} (μ={mu:.3f}, σ={sigma:.3f})"
                    _, _, ph = ax.hist(
                        slice_vals, bins=mx2_bins,
                        histtype='step',
                        color=color,
                        linestyle=ls,
                        label=label
                    )
                    x_fit = np.linspace(0.4, 0.9, 200)
                    ax.plot(x_fit, gauss_poly(x_fit, *popt), color=color, linestyle=ls)

                ax.set(
                    xlim=(0.4, 0.9),
                    xlabel=r'$M_{x (ep)}^{2}$ (GeV²)',
                    ylabel='Counts',
                    title=f'θ: {labels[idx]}°'
                )
                ax.legend(fontsize=8) if any(ax.get_legend_handles_labels()[0]) else ax.text(
                    0.5, 0.5, 'No Data', ha='center', va='center'
                )

            plt.tight_layout()
            out = os.path.join(output_dir, f'Mx2_3particle_{run}_{det_config["name"]}.pdf')
            plt.savefig(out, bbox_inches='tight')
            plt.close()
            print(f"Saved {out}")
        #endfor
    #endfor
#endfor

def gaussian(x, A, mu, sigma):
    """Gaussian function for curve fitting."""
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def plot_dvcs(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) for eppi+pi- channel
    with three particle final state (debug version),
    applies the same kinematic cuts as the C++ code,
    fits a Gaussian to each distribution, overlays the fit,
    and adds μ and σ to the legend.
    """
    detectors = {
        1: {
            'name': 'Forward',
            'theta_bins': [0, 8, 11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 80],
            'theta_labels': ['All θ'] + [f'{b}-{b+3}' for b in [8,11,14,17,20,23,26,29,32,35,38]] + ['38-80']
        },
        2: {
            'name': 'Central',
            'theta_bins': [0, 36, 39, 42, 45, 48, 51, 54, 57, 180],
            'theta_labels': ['All θ'] + [f'{b}-{b+3}' for b in [36,39,42,45,48,51,54]] + ['54-180']
        }
    }

    corrections = ['noCorrections', 'timothy', 'krishna', 'mariana']
    corr_labels = {
        'noCorrections': 'No Corrections',
        'timothy': "Timothy's",
        'krishna': "Krishna's",
        'mariana': "Mariana's"
    }
    colors = ['black', 'red', 'green', 'blue']
    line_styles = ['-', '--', ':', '-.']
    run_periods = ['fa18_inb', 'fa18_out', 'sp19_inb']

    # same Mx² binning as before
    mx2_bins = np.linspace(-0.2, 0.2, 20)

    for run in run_periods:
        print(f"\n{'#'*80}")
        print(f"Processing run period: {run}")
        print(f"{'#'*80}")

        for det_num, det_config in detectors.items():
            print(f"\n{'='*60}")
            print(f"Analyzing {det_config['name']} Detector (ID: {det_num})")
            print(f"Theta bins: {det_config['theta_bins']}")
            print(f"Theta labels: {det_config['theta_labels']}")
            print(f"{'='*60}\n")

            # verify bins/labels
            bins = det_config['theta_bins']
            labels = det_config['theta_labels']
            if len(bins) != len(labels) + 1:
                print("! Critical Error: Bin/Label mismatch!")
                continue

            n_total_plots = len(labels)
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots - 1) / n_cols)) + 1

            all_data = {}

            # load & cut data for each correction
            for corr in corrections:
                fname = f"DVCSWagon_{run}_{corr}.root"
                path = os.path.join(parent_dir, fname)
                print(f"\n{'='*40}")
                print(f"Processing: {fname}")
                if not os.path.exists(path):
                    print(f"! File not found: {path}")
                    continue

                try:
                    with uproot.open(path) as f:
                        if 'PhysicsEvents' not in f:
                            print("! Missing PhysicsEvents tree")
                            continue
                        tree = f['PhysicsEvents']

                        required = [
                            'Mx2_1','p1_theta','detector1',
                            'eta2','t1','theta_gamma_gamma',
                            'Emiss2','pTmiss'
                        ]
                        missing = [b for b in required if b not in tree]
                        if missing:
                            print(f"! Missing branches: {missing}")
                            continue

                        data = tree.arrays(required, library='np')
                        print(f"Total events: {len(data['Mx2_1']):,}")

                        # convert & grab arrays
                        theta = np.degrees(data['p1_theta'])
                        mx2   = data['Mx2_1']
                        det   = data['detector1']
                        eta2  = data['eta2']
                        t1    = data['t1']
                        thgg  = data['theta_gamma_gamma']
                        Emiss = data['Emiss2']
                        pTmis = data['pTmiss']

                        # detector selection
                        mask_det = (det == det_num)
                        print(f"Events in det {det_num}: {np.sum(mask_det):,}")

                        if np.sum(mask_det) == 0:
                            continue

                        # kinematic cuts (as in C++):
                        mask_kin = (
                            (theta >= 5) & (theta < 65) &
                            (eta2 < 0) &
                            (t1 > -2) &
                            (thgg < 0.6) &
                            (Emiss < 0.5) &
                            (pTmis < 0.125)
                        )

                        final_mask = mask_det & mask_kin
                        print(f"After kinematic cuts: {np.sum(final_mask):,}")

                        if np.sum(final_mask) == 0:
                            continue

                        all_data[corr] = {
                            'Mx2_1': mx2[final_mask],
                            'theta': theta[final_mask]
                        }
                        print("Kept ranges:",
                              f"Mx²: {np.min(mx2[final_mask]):.3f}–{np.max(mx2[final_mask]):.3f}",
                              f"θ: {np.min(theta[final_mask]):.1f}°–{np.max(theta[final_mask]):.1f}°")

                except Exception as e:
                    print(f"! Error loading {fname}: {e}")
                    continue

            if not all_data:
                print("! No data for any corrections, skipping plots")
                continue

            # summary printout
            print(f"\n{'='*40}")
            print(f"Data Summary for {det_config['name']} Detector")
            print(f"{'Corr':<15} | {'Events':>10} | {'Mx2 Range':<25} | {'θ Range':<15}")
            print(f"{'-'*70}")
            for corr in corrections:
                if corr in all_data:
                    arr = all_data[corr]['Mx2_1']
                    th  = all_data[corr]['theta']
                    r1, r2 = np.min(arr), np.max(arr)
                    t1_, t2_ = np.min(th), np.max(th)
                    print(f"{corr:<15} | {len(arr):>10,} | "
                          f"{r1:.3f}–{r2:.3f:<17} | {t1_:.1f}°–{t2_:.1f}°")
                else:
                    print(f"{corr:<15} | {'N/A':>10} | {'N/A':<25} | {'N/A':<15}")
            print(f"{'='*40}\n")

            # plotting
            fig = plt.figure(figsize=(20, 5 * n_rows))
            gs  = gridspec.GridSpec(n_rows, n_cols, figure=fig)

            # integrated
            ax0 = fig.add_subplot(gs[0, :])
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr not in all_data:
                    continue
                arr = all_data[corr]['Mx2_1']
                counts, edges = np.histogram(arr, bins=mx2_bins)
                centers = 0.5*(edges[:-1] + edges[1:])
                try:
                    popt, _ = curve_fit(
                        gaussian, centers, counts,
                        p0=[counts.max(), np.mean(arr), np.std(arr)]
                    )
                    mu, sigma = popt[1], abs(popt[2])
                except Exception as e:
                    print(f"! Fit failed for {corr}: {e}")
                    mu, sigma = np.nan, np.nan

                label = f"{corr_labels[corr]} (μ={mu:.3f}, σ={sigma:.3f})"
                _, _, patches = ax0.hist(
                    arr, bins=mx2_bins,
                    histtype='step',
                    color=color, linestyle=ls,
                    label=label
                )
                # overlay fit
                xfit = np.linspace(mx2_bins.min(), mx2_bins.max(), 200)
                ax0.plot(xfit, gaussian(xfit, *popt), color=color, linestyle=ls)

            ax0.set(
                xlabel=r'$M_{x (ep)}^{2}$ (GeV²)',
                ylabel='Counts',
                xlim=(mx2_bins.min(), mx2_bins.max()),
                title=f"{det_config['name']} Detector — {run}"
            )
            ax0.legend()
            ax0.grid(True, alpha=0.3)

            # θ-binned subplots
            for idx in range(1, n_total_plots):
                ax = fig.add_subplot(gs[(idx-1)//n_cols+1, (idx-1)%n_cols])
                tmin, tmax = bins[idx], bins[idx+1]
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr not in all_data:
                        continue
                    th_arr = all_data[corr]['theta']
                    arr    = all_data[corr]['Mx2_1']
                    mask   = (th_arr >= tmin) & (th_arr < tmax)
                    slice_vals = arr[mask]
                    if len(slice_vals) == 0:
                        continue

                    counts, edges = np.histogram(slice_vals, bins=mx2_bins)
                    centers = 0.5*(edges[:-1] + edges[1:])
                    try:
                        popt, _ = curve_fit(
                            gaussian, centers, counts,
                            p0=[counts.max(), np.mean(slice_vals), np.std(slice_vals)]
                        )
                        mu, sigma = popt[1], abs(popt[2])
                    except Exception as e:
                        print(f"! Fit failed for {corr} in θ bin {labels[idx]}: {e}")
                        mu, sigma = np.nan, np.nan

                    label = f"{corr_labels[corr]} (μ={mu:.3f}, σ={sigma:.3f})"
                    _, _, patches = ax.hist(
                        slice_vals, bins=mx2_bins,
                        histtype='step',
                        color=color, linestyle=ls,
                        label=label
                    )
                    xfit = np.linspace(mx2_bins.min(), mx2_bins.max(), 200)
                    ax.plot(xfit, gaussian(xfit, *popt), color=color, linestyle=ls)

                ax.set(
                    xlabel=r'$M_{x (ep)}^{2}$ (GeV²)',
                    ylabel='Counts',
                    xlim=(mx2_bins.min(), mx2_bins.max()),
                    title=f'θ: {labels[idx]}°'
                )
                handles, labs = ax.get_legend_handles_labels()
                if handles:
                    ax.legend(handles=handles, fontsize=8)
                else:
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center')

            # save and close
            out = os.path.join(output_dir, f'Mx2_dvcs_{run}_{det_config["name"]}.pdf')
            plt.savefig(out, bbox_inches='tight')
            plt.close()
            print(f"Saved: {out}")
        #endfor detectors
    #endfor run_periods

if __name__ == "__main__":
    PARENT_DIR = "/volatile/clas12/thayward/corrections_study/results/proton_energy_loss/"
    OUTPUT_DIR = "output/correction_study"
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Generate phase space plots (uncomment to run)
    # for channel in ['ep', 'eppi+pi-']:
    #     for plot_type in ['theta_phi', 'theta_p']:
    #         generate_phase_space_plots(
    #             channel=channel,
    #             correction='noCorrections',
    #             plot_type=plot_type,
    #             parent_dir=PARENT_DIR,
    #             output_dir=OUTPUT_DIR
    #         )

    # # Generate Mx² comparison plots
    # try:
    #     # Generate plots
    #     plot_mx2_comparison(PARENT_DIR, OUTPUT_DIR)
        
    #     # Uncomment to run phase space plots
    #     # generate_phase_space_plots(...)
        
    # finally:
    #     # Clean up matplotlib resources
    #     plt.close('all')

    # Generate Mx² comparison plots
    try:
        # Generate plots
        # plot_three_particles(PARENT_DIR, OUTPUT_DIR)

        plot_dvcs(PARENT_DIR, OUTPUT_DIR)
        
        # Uncomment to run phase space plots
        # generate_phase_space_plots(...)
        
    finally:
        # Clean up matplotlib resources
        plt.close('all')


