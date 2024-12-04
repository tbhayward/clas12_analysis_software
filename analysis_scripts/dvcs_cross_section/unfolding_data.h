// UnfoldingData.h

#ifndef UNFOLDING_DATA_H
#define UNFOLDING_DATA_H

#include <vector>
#include <string>
#include <map>

struct UnfoldingData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;
    std::vector<double> phi_min;
    std::vector<double> phi_max;
    std::vector<double> phi_avg;

    // Raw yields per topology and period:
    std::map<std::string, std::vector<std::vector<int>>> raw_yields; // [topology][period][phi_idx]
    std::vector<std::vector<int>> combined_raw_yields;               // [period][phi_idx]

    // Acceptances, their uncertainties, and unfolded yields per period:
    std::vector<std::vector<double>> acceptance;                 // [period][phi_idx]
    std::vector<std::vector<double>> acceptance_uncertainty;     // [period][phi_idx]
    std::vector<std::vector<double>> unfolded_yields;            // [period][phi_idx]
    std::vector<std::vector<double>> unfolded_yield_uncertainty; // [period][phi_idx]

    // Contamination fractions and signal yields per period:
    std::vector<std::vector<double>> contamination_fraction;    // [period][phi_idx]
    std::vector<std::vector<double>> contamination_uncertainty; // [period][phi_idx]
    std::vector<std::vector<double>> signal_yield;              // [period][phi_idx]
};

#endif // UNFOLDING_DATA_H