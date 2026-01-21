#ifndef NORM_CROSS_SECTIONS_H
#define NORM_CROSS_SECTIONS_H

#include <string>

// Apply overall-normalization factors to the already-computed cross sections:
//
//   input:  "cross sections, ep->epg, exp, <Label>, <hel>" (tuple3)
//   norm:   "overall normalization, ep->epg, exp, <Label>, <hel>" (double OR tuple3)
//   output: "normed cross sections, ep->epg, exp, <Label>, <hel>" (tuple3)
//
// Returns true on success.
bool update_normed_cross_sections_csv(const std::string &csv_main);

// Plot the normalized cross sections exactly like cross_sections.cpp plotting,
// but reading from:
//
//   "normed cross sections, ep->epg, exp, <Label>, <hel>"
//
// Theory curves are unchanged and loaded from the same xs_phi_all.json files.
bool plot_normed_cross_sections_for_label(const std::string &csv_main,
                                          const std::string &label,
                                          const std::string &theory_json_root,
                                          const std::string &out_root_dir);

#endif