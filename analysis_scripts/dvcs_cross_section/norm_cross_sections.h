#ifndef NORM_CROSS_SECTIONS_H
#define NORM_CROSS_SECTIONS_H

#include <string>
#include <vector>

bool update_normed_cross_sections_csv(const std::string &csv_main);

bool update_normed_cross_sections_csv(const std::string &csv_main,
                                      const std::vector<std::string> &labels_to_update);

bool plot_normed_cross_sections_for_label(const std::string &csv_main,
                                          const std::string &label,
                                          const std::string &theory_json_root,
                                          const std::string &out_root_dir);

#endif