import os
import numpy as np
import uproot
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

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

def plot_three_particles(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) for eppi+pi- channel
    with three particle final state (debug version)
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
    
    # mx2_bins = np.linspace(-0.2, 0.2, 100)
    mx2_bins = np.linspace(-0.4, 0.9, 100)

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

            # Verify bin/label alignment
            bins = det_config['theta_bins']
            labels = det_config['theta_labels']
            if len(bins) != len(labels) + 1:
                print(f"! Critical Error: Bin/Label mismatch!")
                print(f"Bins: {len(bins)} items, Labels: {len(labels)} items")
                print("Skipping this detector configuration")
                continue

            n_total_plots = len(labels)
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots-1)/n_cols)) + 1
            
            all_data = {}
            
            # Load data for all corrections
            for corr in corrections:
                filename = f"nSidis_{run}_{corr}.root"
                filepath = os.path.join(parent_dir, filename)
                print(f"\n{'='*40}")
                print(f"Processing: {filename}")
                
                if not os.path.exists(filepath):
                    print(f"! File not found: {filepath}")
                    continue
                
                try:
                    with uproot.open(filepath) as f:
                        # Verify PhysicsEvents tree exists
                        if 'PhysicsEvents' not in f:
                            print("! Missing PhysicsEvents tree")
                            continue
                        tree = f['PhysicsEvents']
                        
                        # Check required branches
                        # required_branches = ['Mx2', 'p1_theta', 'detector1']
                        required_branches = ['Mx2_1', 'p1_theta', 'detector1']
                        missing = [b for b in required_branches if b not in tree]
                        if missing:
                            print(f"! Missing branches: {missing}")
                            continue
                        
                        # Inside the data loading try-block:
                        # Load data
                        data = tree.arrays(required_branches, library='np')
                        # print(f"Total events in file: {len(data['Mx2']):,}")
                        print(f"Total events in file: {len(data['Mx2_1']):,}")

                        # Debug: Check actual detector numbers
                        unique_detectors = np.unique(data['detector1'])
                        print(f"Unique detector1 values: {unique_detectors}")

                        # Apply detector selection
                        mask = (data['detector1'] == det_num)
                        print(f"Events in detector {det_num}: {np.sum(mask):,}")

                        if np.sum(mask) == 0:
                            print("Detector1 value ranges:")
                            print(f"Min: {np.min(data['detector1'])}, Max: {np.max(data['detector1'])}")
                            print("Possible detector numbers:", np.unique(data['detector1']))
                            continue
                            
                        # Store data with debug info
                        all_data[corr] = {
                            # 'Mx2': data['Mx2'][mask],
                            'Mx2_1': data['Mx2_1'][mask],
                            'theta': np.degrees(data['p1_theta'][mask])
                        }
                        print("Data ranges:")
                        # print(f"  Mx²: {np.min(all_data[corr]['Mx2']):.3f} to {np.max(all_data[corr]['Mx2']):.3f} GeV²")
                        print(f"  Mx²: {np.min(all_data[corr]['Mx2_1']):.3f} to {np.max(all_data[corr]['Mx2_1']):.3f} GeV²")
                        print(f"  Theta: {np.min(all_data[corr]['theta']):.1f}° to {np.max(all_data[corr]['theta']):.1f}°")
                        
                except Exception as e:
                    print(f"! Error loading file: {str(e)}")
                    continue

            # Data summary after loading
            print(f"\n{'='*40}")
            print(f"Data Summary for {det_config['name']} Detector")
            # print(f"{'Correction':<15} | {'Events':>10} | {'Mx2 Range':<25} | {'Theta Range':<15}")
            print(f"{'Correction':<15} | {'Events':>10} | {'Mx2_1 Range':<25} | {'Theta Range':<15}")
            print(f"{'-'*70}")
            for corr in corrections:
                if corr in all_data:
                    # mx2_min = np.min(all_data[corr]['Mx2'])
                    # mx2_max = np.max(all_data[corr]['Mx2'])
                    mx2_min = np.min(all_data[corr]['Mx2_1'])
                    mx2_max = np.max(all_data[corr]['Mx2_1'])
                    theta_min = np.min(all_data[corr]['theta'])
                    theta_max = np.max(all_data[corr]['theta'])
                    # print(f"{corr:<15} | {len(all_data[corr]['Mx2']):>10,} | {f'{mx2_min:.3f}-{mx2_max:.3f}':<25} | {f'{theta_min:.1f}°-{theta_max:.1f}°':<15}")
                    print(f"{corr:<15} | {len(all_data[corr]['Mx2_1']):>10,} | {f'{mx2_min:.3f}-{mx2_max:.3f}':<25} | {f'{theta_min:.1f}°-{theta_max:.1f}°':<15}")
                else:
                    print(f"{corr:<15} | {'N/A':>10} | {'N/A':<25} | {'N/A':<15}")
            print(f"{'='*40}\n")

            # Skip plotting if no data
            if not all_data:
                print("! No data available for any corrections - skipping plots")
                continue
                
            # Create figure and plot
            fig = plt.figure(figsize=(20, 5*n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols, figure=fig)
            
            # Integrated spectrum plot
            ax_int = fig.add_subplot(gs[0, :])
            artists = []
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr in all_data:
                    hist = ax_int.hist(
                        # all_data[corr]['Mx2'], bins=mx2_bins,
                        all_data[corr]['Mx2_1'], bins=mx2_bins,
                        histtype='step', color=color, linestyle=ls,
                        label=corr_labels[corr]
                    )
                    artists.extend(hist[2])
            
            if artists:
                ax_int.legend()
                # ax_int.set(xlabel=r'$M_{x}^{2}$ (GeV²)', ylabel='Counts',
                ax_int.set(xlabel=r'$M_{x (ep)}^{2}$ (GeV²)', ylabel='Counts',
                          # xlim=(-0.2, 0.2), title=f"{det_config['name']} Detector - {run}")
                          xlim=(0.4, 0.9), title=f"{det_config['name']} Detector - {run}")
                ax_int.grid(True, alpha=0.3)
            else:
                print("! No data for integrated plot")

            # Theta-binned spectra
            for idx in range(1, n_total_plots):
                row = (idx-1) // n_cols + 1
                col = (idx-1) % n_cols
                ax = fig.add_subplot(gs[row, col])
                
                theta_min = bins[idx]
                theta_max = bins[idx+1]
                current_label = labels[idx]
                
                sub_artists = []
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr in all_data:
                        mask = (all_data[corr]['theta'] >= theta_min) & (all_data[corr]['theta'] < theta_max)
                        # mx2_data = all_data[corr]['Mx2'][mask]
                        mx2_data = all_data[corr]['Mx2_1'][mask]
                        if len(mx2_data) > 0:
                            hist = ax.hist(mx2_data, bins=mx2_bins,
                                         histtype='step', color=color, linestyle=ls,
                                         label=corr_labels[corr])
                            sub_artists.extend(hist[2])
                
                if sub_artists:
                    # ax.set(xlabel=r'$M_{x}^{2}$ (GeV²)', ylabel='Counts',
                    ax.set(xlabel=r'$M_{x (ep)}^{2}$ (GeV²)', ylabel='Counts',
                          # xlim=(-0.2, 0.2), title=f'θ: {current_label}°')
                          xlim=(0.4, 0.9), title=f'θ: {current_label}°')
                    ax.legend(handles=sub_artists, fontsize=8)
                else:
                    ax.set_title(f'θ: {current_label}° (No Data)')
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center')

            output_file = os.path.join(output_dir, f'Mx2_3particle_{run}_{det_config["name"]}.pdf')
            plt.savefig(output_file, bbox_inches='tight')
            plt.close()
            print(f"\nSaved: {output_file}")

def gaussian(x, A, mu, sigma):
    """Gaussian function for curve fitting."""
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

def plot_dvcs(parent_dir, output_dir):
    """
    Analyzes missing mass squared (Mx²) for eppi+pi- channel
    with three particle final state (debug version),
    fits a Gaussian to each distribution, and adds μ and σ to the legend.
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

    # Define histogram bins for Mx²
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

            # Verify bin/label alignment
            bins = det_config['theta_bins']
            labels = det_config['theta_labels']
            if len(bins) != len(labels) + 1:
                print(f"! Critical Error: Bin/Label mismatch!")
                print(f"Bins: {len(bins)} items, Labels: {len(labels)} items")
                print("Skipping this detector configuration")
                continue
            #endif

            n_total_plots = len(labels)
            n_cols = 4
            n_rows = int(np.ceil((n_total_plots - 1) / n_cols)) + 1

            all_data = {}

            # Load data for all corrections
            for corr in corrections:
                filename = f"DVCSWagon_{run}_{corr}.root"
                filepath = os.path.join(parent_dir, filename)
                print(f"\n{'='*40}")
                print(f"Processing: {filename}")

                if not os.path.exists(filepath):
                    print(f"! File not found: {filepath}")
                    continue
                #endif

                try:
                    with uproot.open(filepath) as f:
                        if 'PhysicsEvents' not in f:
                            print("! Missing PhysicsEvents tree")
                            continue
                        #endif

                        tree = f['PhysicsEvents']
                        required_branches = ['Mx2_1', 'p1_theta', 'detector1']
                        missing = [b for b in required_branches if b not in tree]
                        if missing:
                            print(f"! Missing branches: {missing}")
                            continue
                        #endif

                        data = tree.arrays(required_branches, library='np')
                        print(f"Total events in file: {len(data['Mx2_1']):,}")

                        mask = (data['detector1'] == det_num)
                        print(f"Events in detector {det_num}: {np.sum(mask):,}")

                        if np.sum(mask) == 0:
                            print("Detector1 value ranges:")
                            print(f"Min: {np.min(data['detector1'])}, Max: {np.max(data['detector1'])}")
                            print("Possible detector numbers:", np.unique(data['detector1']))
                            continue
                        #endif

                        all_data[corr] = {
                            'Mx2_1': data['Mx2_1'][mask],
                            'theta': np.degrees(data['p1_theta'][mask])
                        }
                        print("Data ranges:")
                        print(f"  Mx²: {np.min(all_data[corr]['Mx2_1']):.3f} to {np.max(all_data[corr]['Mx2_1']):.3f} GeV²")
                        print(f"  Theta: {np.min(all_data[corr]['theta']):.1f}° to {np.max(all_data[corr]['theta']):.1f}°")

                except Exception as e:
                    print(f"! Error loading file: {e}")
                    continue
                #endtry
            #endfor corrections

            # Summary of loaded data
            print(f"\n{'='*40}")
            print(f"Data Summary for {det_config['name']} Detector")
            print(f"{'Correction':<15} | {'Events':>10} | {'Mx2_1 Range':<25} | {'Theta Range':<15}")
            print(f"{'-'*70}")
            for corr in corrections:
                if corr in all_data:
                    mx2_min = np.min(all_data[corr]['Mx2_1'])
                    mx2_max = np.max(all_data[corr]['Mx2_1'])
                    theta_min = np.min(all_data[corr]['theta'])
                    theta_max = np.max(all_data[corr]['theta'])
                    print(f"{corr:<15} | {len(all_data[corr]['Mx2_1']):>10,} | "
                          f"{f'{mx2_min:.3f}-{mx2_max:.3f}':<25} | "
                          f"{f'{theta_min:.1f}°-{theta_max:.1f}°':<15}")
                else:
                    print(f"{corr:<15} | {'N/A':>10} | {'N/A':<25} | {'N/A':<15}")
                #endif
            #endfor corrections
            print(f"{'='*40}\n")

            if not all_data:
                print("! No data available for any corrections - skipping plots")
                continue
            #endif

            # Create figure and gridspec
            fig = plt.figure(figsize=(20, 5 * n_rows))
            gs = gridspec.GridSpec(n_rows, n_cols, figure=fig)

            # Integrated spectrum plot
            ax_int = fig.add_subplot(gs[0, :])
            artists = []
            for corr, color, ls in zip(corrections, colors, line_styles):
                if corr in all_data:
                    vals = all_data[corr]['Mx2_1']
                    counts, hist_edges = np.histogram(vals, bins=mx2_bins)
                    centers = (hist_edges[:-1] + hist_edges[1:]) / 2
                    try:
                        popt, _ = curve_fit(gaussian, centers, counts,
                                            p0=[counts.max(), np.mean(vals), np.std(vals)])
                        mu, sigma = popt[1], abs(popt[2])
                    except Exception as e:
                        print(f"! Gaussian fit failed for {corr}: {e}")
                        mu, sigma = np.nan, np.nan

                    label = f"{corr_labels[corr]} (μ={mu:.3f}, σ={sigma:.3f})"
                    _, _, patches = ax_int.hist(vals, bins=mx2_bins,
                                                histtype='step',
                                                color=color,
                                                linestyle=ls,
                                                label=label)
                    artists.extend(patches)
                #endif
            #endfor integrated

            if artists:
                ax_int.legend()
                ax_int.set(xlabel=r'$M_{x (ep)}^{2}$ (GeV²)',
                           ylabel='Counts',
                           xlim=(-0.2, 0.2),
                           title=f"{det_config['name']} Detector - {run}")
                ax_int.grid(True, alpha=0.3)
            else:
                print("! No data for integrated plot")

            # Theta-binned spectra
            for idx in range(1, n_total_plots):
                row = (idx - 1) // n_cols + 1
                col = (idx - 1) % n_cols
                ax = fig.add_subplot(gs[row, col])

                theta_min = bins[idx]
                theta_max = bins[idx + 1]
                current_label = labels[idx]

                sub_artists = []
                for corr, color, ls in zip(corrections, colors, line_styles):
                    if corr in all_data:
                        mask = ((all_data[corr]['theta'] >= theta_min) &
                                (all_data[corr]['theta'] < theta_max))
                        mx2_data = all_data[corr]['Mx2_1'][mask]
                        if len(mx2_data) > 0:
                            counts, hist_edges = np.histogram(mx2_data, bins=mx2_bins)
                            centers = (hist_edges[:-1] + hist_edges[1:]) / 2
                            try:
                                popt, _ = curve_fit(gaussian, centers, counts,
                                                    p0=[counts.max(),
                                                        np.mean(mx2_data),
                                                        np.std(mx2_data)])
                                mu, sigma = popt[1], abs(popt[2])
                            except Exception as e:
                                print(f"! Gaussian fit failed for {corr} in θ bin {current_label}: {e}")
                                mu, sigma = np.nan, np.nan

                            label = f"{corr_labels[corr]} (μ={mu:.3f}, σ={sigma:.3f})"
                            _, _, patches = ax.hist(mx2_data,
                                                     bins=mx2_bins,
                                                     histtype='step',
                                                     color=color,
                                                     linestyle=ls,
                                                     label=label)
                            sub_artists.extend(patches)
                        #endif
                    #endif
                #endfor theta subplot

                if sub_artists:
                    ax.set(xlabel=r'$M_{x (ep)}^{2}$ (GeV²)',
                           ylabel='Counts',
                           xlim=(-0.2, 0.2),
                           title=f'θ: {current_label}°')
                    ax.legend(handles=sub_artists, fontsize=8)
                else:
                    ax.set_title(f'θ: {current_label}° (No Data)')
                    ax.text(0.5, 0.5, 'No Data', ha='center', va='center')
            #endfor theta-binned

            output_file = os.path.join(output_dir,
                                       f'Mx2_dvcs_{run}_{det_config["name"]}.pdf')
            plt.savefig(output_file, bbox_inches='tight')
            plt.close()
            print(f"\nSaved: {output_file}")
        #endfor detectors
    #endfor runs

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


