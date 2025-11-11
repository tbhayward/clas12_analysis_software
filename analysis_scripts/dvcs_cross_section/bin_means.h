#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include <TTree.h>

#include <map>
#include <string>

/**
 * Update grouped bin-averaged kinematics directly in the analysis CSV.
 *
 * For each valid CSV row (i.e., bin), we:
 *   - Infer the bin grid from columns: xBmin/max, Q2min/max, t_abs_min/max, phimin/phimax
 *   - Determine which topology applies for that row (by inspecting raw-yield columns);
 *     if ambiguous or empty, we accept all topologies.
 *   - Loop events in each DVCS period (in parallel, up to max_workers) and accumulate
 *     sums for xB, Q2, |t|, and phi (with circular mean for phi).
 *   - Write eight grouped-average columns for: Fa18 Inb, Fa18 Out, Sp19 Inb,
 *     Sp18 Inb, Sp18 Out, Fa18 (combined), Sp18 (combined), 10.6 GeV (combined).
 *
 * Notes:
 *   - CSV is updated in place (callers should back up first).
 *   - Data trees are expected in dataTrees keyed by canonical period keys
 *     (e.g. "DVCS_Fa18_inb", "DVCS_Fa18_out", "DVCS_Sp19_inb", "DVCS_Sp18_inb", "DVCS_Sp18_out").
 *   - We auto-detect CSV phi units (degrees vs radians) from phimin/phimax values.
 *   - We call an optional 3-sigma exclusivity cut function if available; otherwise
 *     we fall back to global simple cuts.
 */
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers = 5);

#endif // BIN_MEANS_H