#!/usr/bin/env python3
"""
data_loader.py

Module to load PhysicsEvents trees from ROOT files for each run period
and target type using uproot.
"""

import uproot

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
        }
    }

    trees = {}
    for period, targets in file_paths.items():
        """
        Loop over run periods
        period : str, e.g. "RGC_Su22"
        targets: dict of {target_label: filepath}
        """
        trees[period] = {}
        for target, path in targets.items():
            # Open the ROOT file and grab the PhysicsEvents tree
            root_file = uproot.open(path)
            tree = root_file["PhysicsEvents"]
            trees[period][target] = tree

            # Print confirmation of successful import
            print(f"Successfully loaded PhysicsEvents tree for period '{period}', target '{target}'.")
        #endfor inner loop
    #endfor outer loop

    return trees