// load_run_info_from_csv.h

#ifndef LOAD_RUN_INFO_FROM_CSV_H
#define LOAD_RUN_INFO_FROM_CSV_H

#include <string>
#include <vector>

struct RunInfo {
  int runnum;
  float total_charge;
  float positive_charge;
  float negative_charge;
  float target_polarization;
  float target_polarization_uncertainty;
};

extern std::vector<RunInfo> run_info_list;

void load_run_info_from_csv(const std::string& filename);

#endif // LOAD_RUN_INFO_FROM_CSV_H
