#ifndef UNFOLDING_DATA_H
#define UNFOLDING_DATA_H

#include <vector>
#include <string>

struct UnfoldingData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;
    std::vector<double> phi_min;
    std::vector<double> phi_max;
    // Additional phi_avg vector
    std::vector<double> phi_avg;

    // Raw yields per topology and period:
    // raw_yields[topology][period][phi_idx]
    std::map<std::string, std::vector<std::vector<int>>> raw_yields;

    // Combined raw yields per period:
    // combined_raw_yields[period][phi_idx]
    std::vector<std::vector<int>> combined_raw_yields;

    // Acceptances and unfolded yields per period:
    // acceptance[period][phi_idx]
    std::vector<std::vector<double>> acceptance;
    std::vector<std::vector<double>> unfolded_yields;
};

#endif // UNFOLDING_DATA_H