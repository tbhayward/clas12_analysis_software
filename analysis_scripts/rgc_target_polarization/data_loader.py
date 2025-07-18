#!/usr/bin/env python3
"""
data_loader.py

Module to load PhysicsEvents trees from ROOT files for each run period
and target type using uproot, and report unique Sp23 run numbers.
"""

import uproot
import numpy as np


def load_root_trees():
    """
    Load PhysicsEvents trees from ROOT files for each run period and target type.

    Returns:
        dict: Nested dict of uproot TTree objects:
              trees[period][target] = uproot TTree
    """
    # Define ROOT file paths for each run period and target
    file_paths = {
        "RGC_Su22": {
            "NH3": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_NH3_eX.root",
            "C":   "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_C_eX.root",
            "CH2": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_CH2_eX.root",
            "He":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_He_eX.root",
            "ET":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_su22_inb_ET_eX.root",
        },
        "RGC_Fa22": {
            "NH3": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_NH3_eX.root",
            "C":   "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_C_eX.root",
            "CH2": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_CH2_eX.root",
            "He":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_He_eX.root",
            "ET":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_fa22_inb_ET_eX.root",
        },
        "RGC_Sp23": {
            "NH3": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_NH3_eX.root",
            "C":   "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_C_eX.root",
            "CH2": "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_CH2_eX.root",
            "He":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_He_eX.root",
            "ET":  "/work/clas12/thayward/CLAS12_SIDIS/processed_data/pass2/data/eX/rgc_sp23_inb_ET_eX.root",
        }
    }

    trees = {}
    for period, targets in file_paths.items():
        trees[period] = {}
        for target, path in targets.items():
            # Open the ROOT file and grab the PhysicsEvents tree
            root_file = uproot.open(path)
            tree = root_file["PhysicsEvents"]
            trees[period][target] = tree

            # Print confirmation of successful import
            print(f"Successfully loaded PhysicsEvents tree for period '{period}', target '{target}'.")

    # After loading, report unique run numbers for Sp23
    if "RGC_Sp23" in trees:
        # gather runnum arrays from each Sp23 target
        runnums = []
        for target, tree in trees["RGC_Sp23"].items():
            try:
                arr = tree["runnum"].array(library="np")
                unique_vals = np.unique(arr)
                print(f"Sp23 unique runnum values for target '{target}': {unique_vals}")
                runnums.append(unique_vals)
            except KeyError:
                print(f"Warning: 'runnum' branch not found in Sp23 tree for target '{target}'.")
        # combined unique values across all targets
        if runnums:
            combined = np.unique(np.concatenate(runnums))
            print(f"Combined unique 'runnum' values across all Sp23 targets: {combined}")

    return trees
