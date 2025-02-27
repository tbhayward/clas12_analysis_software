#!/usr/bin/env python3
"""
CLAS12 DVBSA Full Analysis Program
Processes all runs/topologies for DVCS and eppi0 analyses automatically
"""

import os
import json
import logging
import numpy as np
import uproot
import matplotlib.pyplot as plt
from pathlib import Path
from datetime import datetime
from matplotlib.gridspec import GridSpec

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler('dvbsa_analysis.log'),
        logging.StreamHandler()
    ]
)

# Analysis configuration
CONFIG = {
    'base_path': "/work/clas12/thayward/CLAS12_exclusive/",
    'output_dir': "./dvbsa_results/",
    'run_periods': ['fa18_inb', 'fa18_out', 'sp19_inb'],
    'analyses': ['dvcs', 'eppi0'],
    'topologies': {
        '(FD,FD)': (1, 1),
        '(CD,FD)': (2, 1),
        '(CD,FT)': (2, 0)
    },
    'variables': {
        'common': ['open_angle_ep2', 'pTmiss', 'xF', 'Emiss2', 'Mx2', 'Mx2_1', 'Mx2_2'],
        'dvcs': ['theta_gamma_gamma'],
        'eppi0': ['theta_pi0_pi0']
    },
    'hist_config': {
        'theta_gamma_gamma': {'bins': 100, 'range': (0, 180)},
        'theta_pi0_pi0': {'bins': 100, 'range': (0, 180)},
        'open_angle_ep2': {'bins': 100, 'range': (0, 180)},
        'pTmiss': {'bins': 100, 'range': (0, 1.5)},
        'xF': {'bins': 100, 'range': (-1, 1)},
        'Emiss2': {'bins': 100, 'range': (-1, 5)},
        'Mx2': {'bins': 100, 'range': (-1, 5)},
        'Mx2_1': {'bins': 100, 'range': (-1, 5)},
        'Mx2_2': {'bins': 100, 'range': (-1, 5)}
    }
}

def load_files(analysis_type):
    """Load all files for a specific analysis type"""
    files = {}
    base = Path(CONFIG['base_path'])
    
    for run in CONFIG['run_periods']:
        files[run] = {}
        
        # Data files
        data_path = base / f"{analysis_type}/data/pass2/data/{analysis_type}/rga_{run}_epgamma.root"
        if data_path.exists():
            files[run]['data'] = uproot.open(data_path)['PhysicsEvents']
        else:
            logging.warning(f"Missing data file for {run} {analysis_type}")
            continue

        # MC files
        if analysis_type == 'dvcs':
            mc_path = base / f"{analysis_type}/data/pass2/mc/dvcsgen/rec_dvcsgen_rga_{run}_*_epgamma.root"
        else:
            mc_path = base / f"eppi0/data/pass2/mc/hipo_files/rec_aaogen_norad_rga_{run}_*_eppi0.root"
            
        mc_files = list(mc_path.parent.glob(mc_path.name))
        if mc_files:
            files[run]['mc'] = uproot.open(mc_files[0])['PhysicsEvents']
        else:
            logging.warning(f"Missing MC file for {run} {analysis_type}")

    return files

def process_analysis(analysis_type):
    """Process complete analysis for one type (DVCS/eppi0)"""
    logging.info(f"\n{'='*40}\nStarting {analysis_type.upper()} analysis\n{'='*40}")
    
    all_cuts = {}
    files = load_files(analysis_type)
    
    for run in CONFIG['run_periods']:
        if run not in files or 'data' not in files[run] or 'mc' not in files[run]:
            continue
            
        logging.info(f"\nProcessing {run} ({analysis_type.upper()})")
        all_cuts[run] = {}
        
        data_tree = files[run]['data']
        mc_tree = files[run]['mc']
        
        for topology, (d1, d2) in CONFIG['topologies'].items():
            logging.info(f"  Analyzing {topology} topology")
            
            # Create topology masks
            data_mask = (data_tree['detector1'].array() == d1) & (data_tree['detector2'].array() == d2)
            mc_mask = (mc_tree['detector1'].array() == d1) & (mc_tree['detector2'].array() == d2)
            
            # Get variables for this analysis
            variables = CONFIG['variables']['common'] + CONFIG['variables'][analysis_type]
            
            # Process data and MC
            data = {var: data_tree[var].array()[data_mask] for var in variables}
            mc = {var: mc_tree[var].array()[mc_mask] for var in variables}
            
            # Calculate and store cuts
            all_cuts[run][topology] = {
                'data': calculate_cuts(data),
                'mc': calculate_cuts(mc)
            }
            
            # Generate plots
            create_comparison_plots(data, mc, run, topology, analysis_type)
    
    # Save cuts
    output_dir = Path(CONFIG['output_dir']) / analysis_type
    output_dir.mkdir(parents=True, exist_ok=True)
    
    with open(output_dir / 'exclusivity_cuts.json', 'w') as f:
        json.dump(all_cuts, f, indent=2)
    
    logging.info(f"Completed {analysis_type.upper()} analysis")

def calculate_cuts(data):
    """Calculate 3-sigma cuts for all variables"""
    cuts = {}
    for var, values in data.items():
        if var not in CONFIG['hist_config']:
            continue
        mean = np.mean(values)
        std = np.std(values)
        cuts[var] = {
            'mean': float(mean),
            'std': float(std),
            'low': float(mean - 3*std),
            'high': float(mean + 3*std)
        }
    return cuts

def create_comparison_plots(data, mc, run, topology, analysis_type):
    """Create normalized comparison plots"""
    fig = plt.figure(figsize=(20, 10))
    fig.suptitle(f"{run} {topology} ({analysis_type.upper()})", fontsize=14)
    gs = GridSpec(2, 4, figure=fig)
    
    variables = CONFIG['variables']['common'] + CONFIG['variables'][analysis_type]
    
    for i, var in enumerate(variables):
        ax = fig.add_subplot(gs[i//4, i%4])
        config = CONFIG['hist_config'].get(var, {'bins': 50})
        
        # Plot histograms
        counts, bins, _ = ax.hist(data[var], **config, alpha=0.5, label='Data', density=True)
        ax.hist(mc[var], bins=bins, alpha=0.5, label='MC', density=True)
        
        ax.set_xlabel(var.replace('_', ' ').title())
        ax.set_ylabel('Normalized Counts')
        ax.legend()
    
    # Save plot
    output_dir = Path(CONFIG['output_dir']) / analysis_type / run
    output_dir.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_dir / f"{topology}_comparison.png")
    plt.close()

def main():
    """Main entry point - processes everything automatically"""
    start_time = datetime.now()
    logging.info("Starting full DVBSA analysis")
    
    # Process both analysis types
    for analysis in CONFIG['analyses']:
        process_analysis(analysis)
    
    # Generate final report
    logging.info(f"\nAnalysis complete! Total duration: {datetime.now() - start_time}")
    logging.info(f"Results saved to: {Path(CONFIG['output_dir']).resolve()}")

if __name__ == "__main__":
    main()