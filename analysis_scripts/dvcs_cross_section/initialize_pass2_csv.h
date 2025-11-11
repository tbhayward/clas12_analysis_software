#ifndef INITIALIZE_PASS2_CSV_H
#define INITIALIZE_PASS2_CSV_H

#include <string>

/**
 * initialize_pass2_csv
 *
 * Reads Lee's all_bin_v3.csv, filters to rows where "valid bin" == 1,
 * and writes a NEW CSV "dvcs_pass2_analysis.csv" with the extended schema:
 * - exact bin definitions cloned from Lee where requested
 * - helicity-split columns for yields and cross sections (unpol, pos, neg)
 * - contamination, acceptance, BSA as triplets ("value,stat,syst")
 * - prefactors and Fourier coefficients for 10.6 GeV copied from Lee
 * - analogous 10.2 GeV columns left blank
 *
 * All measurement-like columns that must eventually store triplets are created
 * as blank cells (empty string) in this initializer. We only copy simple bin
 * definitions and 10.6 GeV prefactors/coefficients from Lee at this stage.
 *
 * Return: true on success, false on any fatal error.
 */
bool initialize_pass2_csv(const std::string& lee_csv_path,
                          const std::string& out_csv_path);

#endif // INITIALIZE_PASS2_CSV_H