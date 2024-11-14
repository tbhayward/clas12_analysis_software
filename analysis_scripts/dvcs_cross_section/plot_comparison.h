#ifndef PLOT_COMPARISON_H
#define PLOT_COMPARISON_H

#include <string>
#include <vector>

// Struct to hold bin data from the CSV
struct BinData {
    int global_bin_number;
    int bin_number;
    double xBmin, xBmax, xBavg;
    double Q2min, Q2max, Q2avg;
    double tmin, tmax, tavg;
    double phimin, phimax, phiavg;
    double unfolded_yield_inbending;
    double unfolded_yield_outbending;
};

// Main function to manage reading, processing, and printing bin data
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second);

// Helper functions to read each CSV format and print bin data for debugging
std::vector<BinData> read_csv_first(const std::string &file_path);
std::vector<BinData> read_csv_second(const std::string &file_path, const std::vector<BinData> &first_csv_data);
void print_bin_data(const std::vector<BinData> &bins);

#endif // PLOT_COMPARISON_H