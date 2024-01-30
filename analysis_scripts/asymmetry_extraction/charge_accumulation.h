#ifndef CHARGE_ACCUMULATION_H
#define CHARGE_ACCUMULATION_H

#include <vector>
#include "load_run_info_from_csv.h" // Include this if RunInfo is defined here
#include <TTreeReader.h>

// Updated function declaration
void charge_accumulation(TTreeReader& dataReader, const std::vector<RunInfo>& run_info_list);

#endif // CHARGE_ACCUMULATION_H