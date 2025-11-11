#ifndef BIN_MEANS_H
#define BIN_MEANS_H

#include <TTree.h>

#include <map>
#include <string>

/**
 * Update grouped bin-averaged kinematics directly in the analysis CSV.
 *
 * Behavior:
 *   - Reads the bin grid (xB, Q2, |t|, phi) from the CSV; no hardcoded phi bins.
 *   - ALWAYS accepts all three topologies (FD,FD), (CD,FD), (CD,FT) for the averages.
 *   - ALWAYS applies BOTH:
 *        (a) global “simple” cuts and
 *        (b) 3-sigma exclusivity cuts (if exclusivity_cuts.h is available).
 *   - Computes per-row means for xB, Q2, |t|, phi for:
 *        Fa18 Inb, Fa18 Out, Sp19 Inb, Sp18 Inb, Sp18 Out,
 *        and combined Fa18, combined Sp18, combined 10.6 GeV.
 *   - Writes the eight grouped-average columns back into the CSV in place.
 *
 * CSV phi units are auto-detected (degrees if any phimax > 7, else radians).
 */
bool update_bin_means_csv(const std::string& csv_path,
                          const std::map<std::string, TTree*>& dataTrees,
                          int max_workers = 5);

#endif // BIN_MEANS_H