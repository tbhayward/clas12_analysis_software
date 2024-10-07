#ifndef BIN_BOUNDARIES_H
#define BIN_BOUNDARIES_H

#include <map>
#include <string>
#include <vector>

// Function to initialize the bin boundaries
std::map<std::string, std::vector<double>> initialize_bin_boundaries();

// Function to map values to bin indices
std::tuple<int, int, int> map_values_to_bin(double xB, double Q2, double t);

#endif