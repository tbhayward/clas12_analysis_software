#ifndef PLOT_COMPARISON_H
#define PLOT_COMPARISON_H

#include <string>
#include <vector>
#include <tuple>
#include <map>
#include <TCanvas.h>
#include <TGraphErrors.h>

// Define the BinData struct to hold bin information from the CSV
struct BinData {
    int global_bin_number;
    int bin_number;

    double xBmin;
    double xBmax;
    double xBavg;

    double Q2min;
    double Q2max;
    double Q2avg;

    double tmin;
    double tmax;
    double tavg;

    double phimin;
    double phimax;
    double phiavg;

    std::vector<double> signal_yields;              // May have multiple yields
    std::vector<double> signal_yield_uncertainties; // Corresponding uncertainties (Added)
};

// Primary function to manage reading, processing, and plotting of bin data
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second);

// Helper functions for reading CSV files
std::vector<BinData> read_csv_first(const std::string &file_path);
std::vector<BinData> read_csv_second(const std::string &file_path, const std::vector<BinData> &first_csv_data);

// Directory handling function to ensure required output directories exist
void ensure_directory_exists(const std::string &path);

// Function to identify unique xB bins
std::vector<std::pair<double, double>> find_unique_xB_bins(const std::vector<BinData> &data);

// Function to filter data for a specific xB bin range
std::vector<BinData> filter_data_by_xB(const std::vector<BinData> &data, const std::pair<double, double> &xB_range);

// Function to plot data for a specific xB bin with aesthetic enhancements and comparison between two datasets
void plot_for_xB_bin(const std::vector<BinData> &data_first, const std::vector<BinData> &data_second, int xB_index);

#endif // PLOT_COMPARISON_H