#ifndef RADIATIVE_CORRECTIONS_H
#define RADIATIVE_CORRECTIONS_H

#include <map>
#include <string>

class TTree;

/**
 * Compute radiative correction factors Frad for the DVCS bins defined
 * in dvcs_pass2_analysis.csv, and write them into the columns
 *
 *   "Frad, 10.6 GeV"
 *   "Frad, 10.2 GeV"
 *
 * as three-tuples "(value, stat, sys)" with sys = 0.0 for now.
 *
 * The 10.6 GeV factor is filled only for rows that have a non-empty
 * "xBavg, 10.6 GeV" entry.
 *
 * The 10.2 GeV factor is filled only for rows that have a non-empty
 * "xBavg, Sp19 Inb" entry.
 *
 * genMcTrees:    generator MC (no-rad), keys like "DVCS_Sp18_inb_gen"
 * radGenMcTrees: generator MC (radiative), keys like "DVCS_Sp18_inb_gen_rad"
 *
 * out_root_dir:  base output directory (e.g. "output"), used to place
 *                diagnostic radiative_correction_plots.
 *
 * Returns true on success, false on any fatal problem (missing columns,
 * missing MC trees, etc.).
 */
bool update_radiative_corrections_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& radGenMcTrees,
    const std::string& out_root_dir);

#endif