#ifndef NORM_CROSS_SECTIONS_H
#define NORM_CROSS_SECTIONS_H

#include <string>
#include <vector>

// Update CSV cells:
//   "normed cross sections, ep->epg, exp, <Label>, <hel>"
// from:
//   "cross sections, ep->epg, exp, <Label>, <hel>"
// and:
//   "norm, <Label>"
//
// Wrapper overload uses an internal default label list.
bool update_normed_cross_sections_csv(const std::string &csv_main);

bool update_normed_cross_sections_csv(const std::string &csv_main,
                                      const std::vector<std::string> &labels_to_update);

// Plot normed cross sections for a label, overlaying theory curves
// from theory_json_root/<EnergyDir>/xs_phi_all.json keyed by CSV row index.
bool plot_normed_cross_sections_for_label(const std::string &csv_main,
                                          const std::string &label,
                                          const std::string &theory_json_root,
                                          const std::string &out_root_dir);

#endif