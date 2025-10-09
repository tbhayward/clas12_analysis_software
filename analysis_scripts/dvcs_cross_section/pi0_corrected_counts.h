#ifndef PI0_CORRECTED_COUNTS_H
#define PI0_CORRECTED_COUNTS_H

#include <string>
#include <vector>

// Use the project's existing Binning struct
#include "load_binning_scheme.h"

/// Build helicity-resolved π0-corrected counts per φ-bin.
/// Inputs:
///  - periods: DVCS period names (e.g. "DVCS_Fa18_inb", "DVCS_Sp18_out", ...).
///  - binning_scheme: for establishing (ix,iQ,it) grid and φ centers.
///  - total_counts_json: path to output/jsons/total_counts.json (from compute_total_counts).
///  - contamination_json_dir: directory containing contamination_<period>.json files
///       (from compute_pi0_contamination_helicity), typically "output/jsons/contamination".
///  - out_root_dir: root "output" (will write JSONs under output/jsons and plots under
///       output/pi0_corrected_counts/<runTag>/...)
///
/// Output:
///  - Per-period JSON:  output/jsons/pi0_corrected_counts_<period>.json
///  - Per-period canvases: one per xB slice, showing RAW vs CORR counts (± helicities)
///    in each (Q2, t) subpad:  output/pi0_corrected_counts/<runTag>/corr_counts_<period>_xB_<ix>.png
void compute_pi0_corrected_counts(
    const std::vector<std::string>& periods,
    const std::vector<Binning>&     binning_scheme,
    const std::string&              total_counts_json,
    const std::string&              contamination_json_dir,
    const std::string&              out_root_dir);

#endif // PI0_CORRECTED_COUNTS_H