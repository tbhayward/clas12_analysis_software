#ifndef WRITE_CSV_H
#define WRITE_CSV_H

#include <string>
#include <vector>
#include <map>
#include "unfolding_data.h"

// Function to write unfolding data to a CSV file
void write_csv(const std::string& filename, const std::map<std::string, std::vector<UnfoldingData>>& all_unfolding_data);

#endif // WRITE_CSV_H