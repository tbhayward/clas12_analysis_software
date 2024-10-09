#ifndef WRITE_CSV_H
#define WRITE_CSV_H

#include <string>
#include <vector>
#include "plot_unfolding.h"  // Include this to access the `UnfoldingData` struct

// Function prototype for write_csv
void write_csv(const std::string& filename, const std::vector<UnfoldingData>& unfolding_data);

#endif // WRITE_CSV_H