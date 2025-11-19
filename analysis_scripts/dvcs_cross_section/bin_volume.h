#ifndef BIN_VOLUME_H
#define BIN_VOLUME_H

#include <string>

// Compute kinematic bin volumes for each phi-binned row in the CSV and
// write them into the columns:
//   "bin_volume, 10.6 GeV"
//   "bin_volume, 10.2 GeV"
//
// Binning comes from dvcs_pass2_analysis.csv (Lee-style).
// No event input is needed; volumes are computed by deterministic
// integration over (xB, Q2, |t|, phi) with y, W, and t_min masks.
//
// Also produces bin-volume vs phi canvases per beam energy and xB bin under
//   output/bin_volume/10.60
//   output/bin_volume/10.2
//
// Returns true on success, false on failure.
bool update_bin_volume_csv(const std::string& csv_path,
                           const std::string& out_root_dir);

#endif // BIN_VOLUME_H