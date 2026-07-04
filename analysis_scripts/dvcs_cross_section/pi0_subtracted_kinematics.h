#ifndef PI0_SUBTRACTED_KINEMATICS_H
#define PI0_SUBTRACTED_KINEMATICS_H

#include <map>
#include <string>

class TTree;

/**
 * Make DATA/MC comparison plots for selected reconstructed DVCS kinematic
 * branches after the DATA-side pi0 subtraction has been applied, plus companion
 * raw post-exclusivity diagnostics.
 *
 * Placement:
 *   Call this after update_pi0_corrected_counts_csv(...), because the
 *   pi0-subtracted mode reads the CSV signal-yield columns:
 *
 *       signal yield, ep->epg, exp, <period>, unpol
 *
 * DATA weighting, pi0-subtracted mode:
 *   For every accepted epgamma DATA event, the event is assigned to its
 *   xB-Q2-|t|-phi CSV row. The event weight is
 *
 *       w_data(row, period) = S_unpol(row, period) / N_raw_unpol(row, period)
 *
 *   where S_unpol is the pi0-corrected signal yield from the CSV and
 *   N_raw_unpol is the sum of raw epgamma DATA counts over the three topologies
 *   in that same row and period. This makes the total DATA weight in each 4D
 *   bin equal to the pi0-subtracted signal yield while preserving the observed
 *   within-bin shape of the plotted branch.
 *
 * DATA weighting, raw mode:
 *   Accepted epgamma DATA events are filled with unit weight after the same
 *   global and topology-dependent 3sigma exclusivity cuts.
 *
 * MC weighting:
 *   Reconstructed DVCS MC events are filled with unit weight after the same
 *   global and topology-dependent 3sigma exclusivity cuts. Each DATA and MC
 *   histogram is normalized to its own integral before drawing.
 *
 * Output subdirectories below out_dir:
 *
 *   inclusive/pi0_subtracted/
 *   inclusive/raw_post_exclusivity/
 *   topology/pi0_subtracted/FD_FD/, CD_FD/, CD_FT/
 *   topology/raw_post_exclusivity/FD_FD/, CD_FD/, CD_FT/
 *   electron_sector_sp18_out/pi0_subtracted/
 *   electron_sector_sp18_out/raw_post_exclusivity/
 *   diagnostics/
 *
 * Each plot set contains one 2x3 canvas for each branch:
 *
 *   e_p, e_theta, p1_p, p1_theta, p2_p, p2_theta
 *
 * Angles are converted from radians to degrees.
 */
bool plot_pi0_subtracted_dvcs_kinematics(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& out_dir,
    int max_workers = 5);

#endif // PI0_SUBTRACTED_KINEMATICS_H
