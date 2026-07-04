#ifndef PI0_SUBTRACTED_KINEMATICS_H
#define PI0_SUBTRACTED_KINEMATICS_H

#include <map>
#include <string>

class TTree;

/**
 * Make DATA/MC comparison plots for selected reconstructed DVCS kinematic
 * branches after the DATA-side pi0 subtraction has been applied.
 *
 * Placement:
 *   Call this after update_pi0_corrected_counts_csv(...), because this routine
 *   reads the CSV signal-yield columns:
 *
 *       signal yield, ep->epg, exp, <period>, unpol
 *
 * DATA weighting:
 *   For every accepted epgamma DATA event, the event is assigned to its
 *   xB-Q2-|t|-phi CSV row.  The event weight is
 *
 *       w_data(row, period) = S_unpol(row, period) / N_raw_unpol(row, period)
 *
 *   where S_unpol is the pi0-corrected signal yield from the CSV and
 *   N_raw_unpol is the sum of the raw epgamma DATA counts over the three
 *   topologies in that same row and period.  This makes the total DATA weight
 *   in each 4D bin equal to the pi0-subtracted signal yield while preserving
 *   the observed within-bin shape of the plotted branch.
 *
 * MC weighting:
 *   Reconstructed DVCS MC events are filled with unit weight after the same
 *   global and topology-dependent 3sigma exclusivity cuts.  Each DATA and MC
 *   histogram is normalized to its own integral before drawing.
 *
 * Output:
 *   One 2x3 canvas per plotted branch is written to out_dir:
 *
 *       e_p_pi0_subtracted_data_vs_rec_mc.png
 *       e_theta_pi0_subtracted_data_vs_rec_mc.png
 *       p1_p_pi0_subtracted_data_vs_rec_mc.png
 *       p1_theta_pi0_subtracted_data_vs_rec_mc.png
 *       p2_p_pi0_subtracted_data_vs_rec_mc.png
 *       p2_theta_pi0_subtracted_data_vs_rec_mc.png
 */
bool plot_pi0_subtracted_dvcs_kinematics(
    const std::string& csv_path,
    const std::map<std::string, TTree*>& dvcsDataTrees,
    const std::map<std::string, TTree*>& dvcsRecMcTrees,
    const std::string& combined_cuts_json,
    const std::string& out_dir,
    int max_workers = 5);

#endif // PI0_SUBTRACTED_KINEMATICS_H
