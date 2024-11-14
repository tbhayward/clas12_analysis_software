#include "plot_comparison.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

// Helper function to read data from the first CSV file format
std::vector<BinData> read_csv_first(const std::string &file_path) {
    std::vector<BinData> bins;
    std::ifstream file(file_path);
    std::string line;

    // Skip the header line
    std::getline(file, line);

    // Read each line of the CSV
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        BinData bin;

        // Read values based on column order for the first CSV format
        std::getline(ss, value, ','); // First unnamed column
        bin.global_bin_number = std::stoi(value);

        std::getline(ss, value, ','); // Bin name column
        bin.bin_number = std::stoi(value);

        std::getline(ss, value, ','); // xBmin
        bin.xBmin = std::stod(value);
        std::getline(ss, value, ','); // xBmax
        bin.xBmax = std::stod(value);
        std::getline(ss, value, ','); // xBavg
        bin.xBavg = std::stod(value);

        std::getline(ss, value, ','); // Q2min
        bin.Q2min = std::stod(value);
        std::getline(ss, value, ','); // Q2max
        bin.Q2max = std::stod(value);
        std::getline(ss, value, ','); // Q2avg
        bin.Q2avg = std::stod(value);

        std::getline(ss, value, ','); // t_abs_min (tmin)
        bin.tmin = std::stod(value);
        std::getline(ss, value, ','); // t_abs_max (tmax)
        bin.tmax = std::stod(value);
        std::getline(ss, value, ','); // t_abs_avg (tavg)
        bin.tavg = std::stod(value);

        std::getline(ss, value, ','); // phimin
        bin.phimin = std::stod(value);
        std::getline(ss, value, ','); // phimax
        bin.phimax = std::stod(value);
        std::getline(ss, value, ','); // phiavg
        bin.phiavg = std::stod(value);

        // Skip intermediate columns until reaching "signal yield, ep->epg, exp, inbending"
        for (int i = 0; i < 21; ++i) std::getline(ss, value, ',');

        std::getline(ss, value, ','); // signal yield, ep->epg, exp, inbending
        bin.unfolded_yield_inbending = std::stod(value);

        std::getline(ss, value, ','); // signal yield, ep->epg, exp, outbending
        bin.unfolded_yield_outbending = std::stod(value);

        bins.push_back(bin);
    }

    return bins;
}

// Helper function to read data from the second CSV file format
std::vector<BinData> read_csv_second(const std::string &file_path) {
    std::vector<BinData> bins;
    std::ifstream file(file_path);
    std::string line;

    // Skip the header line
    std::getline(file, line);

    // Read each line of the CSV
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string value;

        BinData bin;

        // Read values based on column order for the second CSV format
        std::getline(ss, value, ','); // First unnamed column
        bin.global_bin_number = std::stoi(value);

        std::getline(ss, value, ','); // Bin name column
        bin.bin_number = std::stoi(value);

        std::getline(ss, value, ','); // xBmin
        bin.xBmin = std::stod(value);
        std::getline(ss, value, ','); // xBmax
        bin.xBmax = std::stod(value);
        std::getline(ss, value, ','); // xBavg
        bin.xBavg = std::stod(value);

        std::getline(ss, value, ','); // Q2min
        bin.Q2min = std::stod(value);
        std::getline(ss, value, ','); // Q2max
        bin.Q2max = std::stod(value);
        std::getline(ss, value, ','); // Q2avg
        bin.Q2avg = std::stod(value);

        std::getline(ss, value, ','); // t_abs_min (tmin)
        bin.tmin = std::stod(value);
        std::getline(ss, value, ','); // t_abs_max (tmax)
        bin.tmax = std::stod(value);
        std::getline(ss, value, ','); // t_abs_avg (tavg)
        bin.tavg = std::stod(value);

        std::getline(ss, value, ','); // phimin
        bin.phimin = std::stod(value);
        std::getline(ss, value, ','); // phimax
        bin.phimax = std::stod(value);
        std::getline(ss, value, ','); // phiavg
        bin.phiavg = std::stod(value);

        // Skip to columns "AJ" and "AK" for unfolded yields inbending and outbending
        for (int i = 0; i < 35; ++i) std::getline(ss, value, ',');

        std::getline(ss, value, ','); // "ep->e'pgamma unfolded_yield_Fa18Inb"
        bin.unfolded_yield_inbending = std::stod(value);

        std::getline(ss, value, ','); // "ep->e'pgamma unfolded_yield_Fa18Out"
        bin.unfolded_yield_outbending = std::stod(value);

        bins.push_back(bin);
    }

    return bins;
}

// Helper function to print out values in the struct for debugging
void print_bin_data(const std::vector<BinData> &bins) {
    for (const auto &bin : bins) {
        std::cout << "Global Bin Number: " << bin.global_bin_number << '\n';
        std::cout << "Bin Number: " << bin.bin_number << '\n';
        std::cout << "xBmin: " << bin.xBmin << ", xBmax: " << bin.xBmax << ", xBavg: " << bin.xBavg << '\n';
        std::cout << "Q2min: " << bin.Q2min << ", Q2max: " << bin.Q2max << ", Q2avg: " << bin.Q2avg << '\n';
        std::cout << "tmin: " << bin.tmin << ", tmax: " << bin.tmax << ", tavg: " << bin.tavg << '\n';
        std::cout << "phimin: " << bin.phimin << ", phimax: " << bin.phimax << ", phiavg: " << bin.phiavg << '\n';
        std::cout << "Unfolded Yield Inbending: " << bin.unfolded_yield_inbending << '\n';
        std::cout << "Unfolded Yield Outbending: " << bin.unfolded_yield_outbending << '\n';
        std::cout << "----------------------\n";
    }
}

// Main function to control the import, processing, and debug output
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second) {
    // Step 1: Read the first CSV data into a vector of BinData
    std::vector<BinData> bin_data_first = read_csv_first(csv_file_path_first);

    // Optional: Print out the first CSV data for verification (commented out for now)
    // print_bin_data(bin_data_first);

    // Step 2: Read the second CSV data into a vector of BinData
    std::vector<BinData> bin_data_second = read_csv_second(csv_file_path_second);

    // Step 3: Print out the second CSV data to verify correctness
    print_bin_data(bin_data_second);
}