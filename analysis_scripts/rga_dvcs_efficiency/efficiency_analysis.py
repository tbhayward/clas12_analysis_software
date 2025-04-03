def plot_normalized_efficiencies(output_dir):
    """
    Calculates and plots normalized efficiencies as a function of beam current
    for each run period. For each current x nA, the efficiency is calculated as:

        efficiency = (n_reco(x nA)/n_gen(x nA)) / (n_reco(0 nA)/n_gen(0 nA))

    The results are then fitted to a linear polynomial forced through (0,1), i.e.,
    y = m*x + 1. Three separate canvases (one per run period) are saved as PDFs.
    """
    import os
    import glob
    import numpy as np
    import matplotlib.pyplot as plt

    # Define base directory where the ROOT files are stored.
    base_dir = "/work/clas12/thayward/CLAS12_exclusive/dvcs/data/pass2/efficiency_study"
    
    # Define the run periods with their file prefix mapping.
    run_periods = {
        "rga_fa18_inb": "DVCS_Fa18_inb",
        "rga_fa18_out": "DVCS_Fa18_out",
        "rga_sp19_inb": "DVCS_Sp19_inb"
    }
    
    # Map run_prefix to a more stylized title.
    title_map = {
        "rga_fa18_inb": "RGA Fa18 Inb",
        "rga_fa18_out": "RGA Fa18 Out",
        "rga_sp19_inb": "RGA Sp19 Inb"
    }
    
    # Define the selected topology and analysis type.
    selected_topology = "(FD,FD)"
    analysis_type = "dvcs"
    
    # Loop over each run period.
    for run_prefix, period_code in run_periods.items():
        # Find all generated and reconstructed ROOT files for this run period.
        # Use a less restrictive pattern to include both "nobkg" and "XXnA" files.
        pattern_gen = os.path.join(base_dir, f"gen_{run_prefix}_dvcs_*.root")
        pattern_reco = os.path.join(base_dir, f"{run_prefix}_dvcs_*.root")
        gen_files = glob.glob(pattern_gen)
        reco_files = glob.glob(pattern_reco)
        
        # Create dictionaries mapping current (in nA) to file path for both gen and reco.
        gen_dict = {}
        reco_dict = {}
        
        def extract_current(filename):
            """
            Extracts the beam current from the filename.
            Expected pattern: ..._{current}nA.root, where current is a number
            or 'nobkg' (which is mapped to 0).
            """
            import os
            base = os.path.basename(filename)
            parts = base.split("_")
            last_part = parts[-1].replace(".root", "")
            if "nobkg" in last_part:
                return 0
            else:
                try:
                    current_str = last_part.replace("nA", "")
                    return int(current_str)
                except:
                    return None
                #endif
            #endif
        #enddef

        for fpath in gen_files:
            cur = extract_current(fpath)
            if cur is not None:
                gen_dict[cur] = fpath
            #endif
        #endfor

        for fpath in reco_files:
            cur = extract_current(fpath)
            if cur is not None:
                reco_dict[cur] = fpath
            #endif
        #endfor

        # Check that we have a 0 nA (nobkg) file for normalization.
        if 0 not in gen_dict or 0 not in reco_dict:
            print(f"Missing nobkg files for run period {run_prefix}. Skipping.")
            continue
        #endif
        
        # Load the exclusivity cuts for this period.
        cuts = load_cuts(period_code, selected_topology)
        
        # Calculate efficiency for each current.
        currents = []
        efficiencies = []
        efficiency_errors = []
        
        # Calculate normalization factor from the nobkg files (0 nA).
        n_gen_0 = get_total_tree_entries(gen_dict[0])
        n_reco_0 = get_tree_entries_with_cuts(reco_dict[0], cuts, selected_topology, analysis_type, period_code)
        if n_gen_0 == 0:
            print(f"Warning: 0 generated entries for normalization in run period {run_prefix}")
            norm_eff = 1.0
        else:
            norm_eff = (n_reco_0 / n_gen_0)
        #endif

        # Loop over all currents (ensuring 0 is included for normalization).
        for cur in sorted(gen_dict.keys()):
            if cur == 0:
                currents.append(0)
                efficiencies.append(1.0)  # by definition, normalized to itself
                efficiency_errors.append(0.0)
                continue
            #endif
            n_gen = get_total_tree_entries(gen_dict[cur])
            n_reco = get_tree_entries_with_cuts(reco_dict[cur], cuts, selected_topology, analysis_type, period_code)
            
            if n_gen == 0:
                eff = 0
                err = 0
            else:
                eff = (n_reco / n_gen) / norm_eff
                # Estimate error from Poisson statistics (assuming error in n_reco dominates).
                if n_reco > 0:
                    err = (np.sqrt(n_reco) / n_gen) / norm_eff
                else:
                    err = 0
                #endif
            #endif
            currents.append(cur)
            efficiencies.append(eff)
            efficiency_errors.append(err)
        #endfor
        
        # Convert lists to numpy arrays.
        currents_arr = np.array(currents)
        efficiencies_arr = np.array(efficiencies)
        
        # For weighted least squares, use weights = 1/err^2.
        weights2 = np.array([1/(err**2) if err > 0 else 1 for err in efficiency_errors])
        
        # Force the fit through (0,1): model: efficiency = m*x + 1.
        # Calculate best-fit slope and its uncertainty.
        # m = sum(w * x*(y-1)) / sum(w * x^2)
        numerator = np.sum(weights2 * currents_arr * (efficiencies_arr - 1))
        denominator = np.sum(weights2 * currents_arr**2)
        m = numerator / denominator if denominator != 0 else 0
        sigma_m = np.sqrt(1 / denominator) if denominator != 0 else 0
        
        # Calculate chi2:
        model = 1 + m * currents_arr
        chi2 = np.sum(weights2 * (efficiencies_arr - model)**2)
        ndf = len(currents_arr) - 1  # one parameter fit
        chi2_ndf = chi2 / ndf if ndf > 0 else 0
        
        # Generate fitted line data.
        fit_x = np.linspace(min(currents_arr), max(currents_arr), 100)
        fit_y = 1 + m * fit_x
        
        # Create the plot.
        plt.figure()
        plt.errorbar(currents_arr, efficiencies_arr, yerr=efficiency_errors, fmt='ko', label='Data')
        plt.plot(fit_x, fit_y, 'r--', label=f"m = {m:.3e} ± {sigma_m:.3e}\n" + r"$\chi^{2}/ndf$ = " + f"{chi2_ndf:.2f}")
        plt.xlabel("current (nA)")
        plt.ylabel("normalized efficiency")
        plt.xlim(-5, 60)
        plt.ylim(0.7, 1.1)
        plt.legend(loc='upper right')
        plt.title(f"Normalized Efficiency for {title_map.get(run_prefix, run_prefix)}")
        
        # Save the plot as a PDF.
        output_path = os.path.join(output_dir, f"{run_prefix}_integrated.pdf")
        plt.savefig(output_path, format='pdf')
        plt.close()
        print(f"Saved plot for {run_prefix} to {output_path}")
    #endfor
#enddef