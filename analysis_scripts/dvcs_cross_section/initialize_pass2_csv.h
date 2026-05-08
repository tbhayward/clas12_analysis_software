#ifndef INITIALIZE_PASS2_CSV_H
#define INITIALIZE_PASS2_CSV_H

#include <string>

/**
 * initialize_pass2_csv
 *
 * Reads Lee's all_bin_v3.csv, filters to rows where "valid bin" == 1,
 * and writes a new pass-2 analysis CSV with the extended schema.
 *
 * The initializer copies only stable bin-definition and theory-prefactor
 * information from Lee's CSV. All analysis-produced quantities are created as
 * blank cells and must be filled by later modules.
 *
 * The generated schema includes columns for:
 *
 * - grouped average kinematics,
 * - raw DVCS and eppi0 yields,
 * - channel-aware current-efficiency factors for data and MC,
 * - eppi0 p1_theta data/MC cross-section normalization factors for six FD sectors plus CD,
 * - normalized raw DVCS and eppi0 yields,
 * - generated and reconstructed MC yields for DVCS and eppi0,
 * - pi0 contamination,
 * - signal yields,
 * - acceptance,
 * - acceptance-corrected yields,
 * - radiative corrections,
 * - bin-centering corrections,
 * - bin volumes,
 * - luminosities,
 * - cross sections,
 * - BH-normalization factors,
 * - normed cross sections,
 * - BSA quantities,
 * - copied 10.6 GeV theory/prefactor columns.
 *
 * Current-efficiency factor cells are intended to store:
 *
 *   (value,stat)
 *
 * eppi0 cross-section normalization cubic cells are intended to store:
 *
 *   (p0,p1,p2,p3)
 *
 * where:
 *
 *   R_pi0(theta_p) = p0 + p1*theta_p + p2*theta_p^2 + p3*theta_p^3
 *
 * and downstream data yield corrections divide event weights by R_pi0(theta_p).
 *
 * Return:
 *   true on success, false on any fatal error.
 */
bool initialize_pass2_csv(const std::string& lee_csv_path,
                          const std::string& out_csv_path);

#endif // INITIALIZE_PASS2_CSV_H