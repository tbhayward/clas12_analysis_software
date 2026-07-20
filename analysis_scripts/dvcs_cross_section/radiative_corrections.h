#ifndef RADIATIVE_CORRECTIONS_H
#define RADIATIVE_CORRECTIONS_H

#include <map>
#include <string>

class TTree;

// Compute Frad and write it into the CSV columns
//   "Frad, 10.6 GeV" and "Frad, 10.2 GeV"
// as "(value, stat, sys)" triples.
//
// Uses generated Born and radiative MC grouped by beam energy:
//
//   10.6 GeV: Sp18 Inb, Sp18 Out, Fa18 Inb, Fa18 Out
//   10.2 GeV: Sp19 Inb
//
// For each energy group, Frad in a CSV row i is defined as
//   Frad_i = (a_i / N_born) / (b_i / N_rad)
//          = a_i * N_rad / (b_i * N_born)
// where a_i and b_i are the Born and Rad counts in that row, and
// N_born, N_rad are the global totals of generated events in that
// group (summing all periods).
//
// Bins come from dvcs_pass2_analysis.csv (Lee-style).
//
// Generated Born and radiative MC use generator-applicable global cuts and
// analysis-bin phase space only. Reconstructed exclusivity cuts are not
// applied at generator level.
//
// Returns true on success, false on any fatal failure.
bool update_radiative_corrections_csv(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& genMcTrees,
    const std::map<std::string, TTree*>& radGenMcTrees,
    const std::string& out_root_dir);

#endif // RADIATIVE_CORRECTIONS_H