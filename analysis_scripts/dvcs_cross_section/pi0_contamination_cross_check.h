#ifndef PI0_CONTAMINATION_CROSS_CHECK_H
#define PI0_CONTAMINATION_CROSS_CHECK_H

#include <string>
#include <vector>

// Forward: uses LeeRow from load_csv.h
struct LeeRow;

/**
 * Compare our helicity-averaged pi0 contamination (from JSON)
 * to colleague values (from all_bin_v3.csv).
 *
 * - pi0_combined_json must be the combined JSON produced by your
 *   regular analysis (pi0_contamination_combined.json).
 * - output_base_dir is where plots are written, e.g.
 *   output/cross_check/lee/pi0_contamination
 */
void plot_pi0_contam_cross_checks(const std::vector<LeeRow>& rows,
                                  const std::string& pi0_combined_json,
                                  const std::string& output_base_dir);

#endif // PI0_CONTAMINATION_CROSS_CHECK_H