#ifndef WRITE_CSV_H
#define WRITE_CSV_H

#include <string>
#include <vector>
#include "plot_unfolding.h"  // Assuming the new structure YieldData is defined here

// Function to write yield data to a CSV file
void write_csv(std::vector<YieldData> &yields, const std::string &filename);

#endif // WRITE_CSV_H