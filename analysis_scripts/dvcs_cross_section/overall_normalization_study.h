#ifndef OVERALL_NORMALIZATION_STUDY_H
#define OVERALL_NORMALIZATION_STUDY_H

#include <string>

// Print a per-(xB,Q2,t) bin study of how close-to-BH the measured cross section is
// at the phi point closest to 0 or 360 degrees.
//
// Required CSV columns (fail-fast if missing):
//   xBmin, xBmax, Q2min, Q2max, t_abs_min, t_abs_max
//   xBavg, <label>
//   Q2avg, <label>
//   t_abs_avg, <label>
//   phiavg, <label>
//   cross sections, ep->epg, exp, <label>, <helicity>
//
// Example label: "10.6 GeV"
// Example helicity: "unpol"
//
// Returns true on success, false on failure (with messages to stderr).
bool print_bh_normalization_study(const std::string &csv_path,
                                  const std::string &label,
                                  const std::string &helicity);

#endif