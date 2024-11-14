#ifndef PLOT_COMPARISON_H
#define PLOT_COMPARISON_H

#include <string>
#include <vector>

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

// Main function to manage reading, processing, and plotting bin data
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second);

// Helper functions
std::vector<BinData> read_csv_first(const std::string &file_path);
std::vector<BinData> read_csv_second(const std::string &file_path, const std::vector<BinData> &first_csv_data);
void print_bin_data(const std::vector<BinData> &bins);
std::vector<std::pair<double, double>> find_unique_xB_bins(const std::vector<BinData> &data);
std::vector<BinData> filter_data_by_xB(const std::vector<BinData> &data, const std::pair<double, double> &xB_range);
void plot_for_xB_bin(const std::vector<BinData> &data, int xB_index);

// New helper function to ensure directories exist
void ensure_directory_exists(const std::string &path);

#endif // PLOT_COMPARISON_H