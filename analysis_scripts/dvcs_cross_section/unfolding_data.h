// unfolding_data.h
#ifndef UNFOLDING_DATA_H
#define UNFOLDING_DATA_H

#include <vector>

// Struct to hold the unfolding data for each phi bin within an xB-Q2-t bin
struct UnfoldingData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;

    // Phi bin data (min and max for each bin)
    std::vector<double> phi_min;  // Lower boundary for each phi bin
    std::vector<double> phi_max;  // Upper boundary for each phi bin

    // 2D Vectors to store raw yields, acceptance, and unfolded yields for each period and topology
    std::vector<std::vector<int>> raw_yields;        // [period][topology]
    std::vector<std::vector<double>> acceptance;     // [period]
    std::vector<std::vector<double>> unfolded_yields; // [period]

    std::vector<std::vector<double>> contamination_ratio; // [period][phi_bin]
    std::vector<std::vector<double>> contamination_error; // [period][phi_bin]
};

#endif // UNFOLDING_DATA_H