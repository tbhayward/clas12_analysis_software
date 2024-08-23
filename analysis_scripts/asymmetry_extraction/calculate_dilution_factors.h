#ifndef CALCULATE_DILUTION_FACTORS_H
#define CALCULATE_DILUTION_FACTORS_H

#include <vector>
#include <utility>

// Declaration of the calculate_dilution_factors function
std::vector<std::pair<double, double>> calculate_dilution_factors();

// Declaration of the calculate_dilution_error function
double calculate_dilution_error(double nA, double nC, double nCH, double nMT, double nf);

#endif