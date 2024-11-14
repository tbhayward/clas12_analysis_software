#include "plot_comparison.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>     
#include <set>    
#include <sys/stat.h> 
#include <unistd.h>   

// ROOT includes for plotting
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TAxis.h>
#include <TColor.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TString.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLegend.h>

// Helper function to ensure a directory exists
void ensure_directory_exists(const std::string &path) {
    struct stat info;
    if (stat(path.c_str(), &info) != 0) {
        // Directory does not exist, so create it
        mkdir(path.c_str(), 0777);
    }
}

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

// Helper function to read data from the second CSV file format, using bin information from the first CSV
std::vector<BinData> read_csv_second(const std::string &file_path, const std::vector<BinData> &first_csv_data) {
    std::vector<BinData> bins;
    std::ifstream file(file_path);
    std::string line;

    // Skip the header line
    std::getline(file, line);

    // Read each line of the CSV
    int index = 0;
    while (std::getline(file, line) && index < first_csv_data.size()) {
        std::stringstream ss(line);
        std::string value;

        BinData bin;

        // Assign global_bin_number and bin_number based on the first CSV
        bin.global_bin_number = first_csv_data[index].global_bin_number;
        bin.bin_number = first_csv_data[index].bin_number;

        // Read values based on column order for the second CSV format
        std::getline(ss, value, ','); // Skip Bin number column
        // The bin number is already assigned from first_csv_data

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

        std::getline(ss, value, ','); // t_min (tmin)
        bin.tmin = std::stod(value);
        std::getline(ss, value, ','); // t_max (tmax)
        bin.tmax = std::stod(value);
        std::getline(ss, value, ','); // t_avg (tavg)
        bin.tavg = std::stod(value);

        std::getline(ss, value, ','); // phimin
        bin.phimin = std::stod(value);
        std::getline(ss, value, ','); // phimax
        bin.phimax = std::stod(value);
        std::getline(ss, value, ','); // phiavg
        bin.phiavg = std::stod(value);

        // Skip columns until reaching unfolded_yield_inbending (column S)
        for (int i = 0; i < 5; ++i) std::getline(ss, value, ',');

        std::getline(ss, value, ','); // ep->e'pgamma unfolded_yield_Fa18Inb
        bin.unfolded_yield_inbending = std::stod(value);


        // Skip to unfolded_yield_outbending (column Y)
        for (int i = 0; i < 5; ++i) std::getline(ss, value, ',');

        std::getline(ss, value, ','); // ep->e'pgamma unfolded_yield_Fa18Out
        bin.unfolded_yield_outbending = std::stod(value);

        bins.push_back(bin);
        index++;
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

// Function to find unique xB bins in the data
std::vector<std::pair<double, double>> find_unique_xB_bins(const std::vector<BinData> &data) {
    std::set<std::pair<double, double>> unique_bins;
    for (const auto &bin : data) {
        unique_bins.emplace(bin.xBmin, bin.xBmax);
    }
    return std::vector<std::pair<double, double>>(unique_bins.begin(), unique_bins.end());
}

// Function to filter data for a specific xB bin
std::vector<BinData> filter_data_by_xB(const std::vector<BinData> &data, const std::pair<double, double> &xB_range) {
    std::vector<BinData> filtered_data;
    for (const auto &bin : data) {
        if (bin.xBmin == xB_range.first && bin.xBmax == xB_range.second) {
            filtered_data.push_back(bin);
        }
    }
    return filtered_data;
}

// Plotting function for each xB bin
void plot_for_xB_bin(const std::vector<BinData> &data, int xB_index, const std::string &type) {
    // Step 1: Identify unique (Q2min, Q2max, tmin, tmax) bins
    std::map<std::tuple<double, double, double, double>, std::vector<BinData>> qt_bins;
    for (const auto &bin : data) {
        // Use a tuple of (Q2min, Q2max, tmin, tmax) to identify unique bins
        auto key = std::make_tuple(bin.Q2min, bin.Q2max, bin.tmin, bin.tmax);
        qt_bins[key].push_back(bin);
    }

    // Step 2: Determine the grid size for subplots
    int num_plots = qt_bins.size();
    int grid_size = std::ceil(std::sqrt(num_plots));

    // Create a canvas and divide it into a grid
    TCanvas canvas("canvas", "Cross Check", 1200, 1200);
    canvas.Divide(grid_size, grid_size);

    int plot_index = 1;
    for (const auto &[qt_key, bins] : qt_bins) {
        // Move to the next pad
        canvas.cd(plot_index);

        // Prepare vectors for phiavg and unfolded yield values
        std::vector<double> phi_values;
        std::vector<double> yields;
        std::vector<double> yield_errors;

        for (const auto &bin : bins) {
            phi_values.push_back(bin.phiavg);
            yields.push_back(bin.unfolded_yield_outbending);
            yield_errors.push_back(std::sqrt(bin.unfolded_yield_outbending)); // Simple sqrt error
        }

        // Create TGraphErrors for each (Q2, t) bin
        int n_points = phi_values.size();
        TGraphErrors *graph = new TGraphErrors(n_points, &phi_values[0], &yields[0], nullptr, &yield_errors[0]);
        
        // Set the title based on xB, Q2, and t bin
        double xB_avg = bins[0].xBavg;
        double Q2avg = (std::get<0>(qt_key) + std::get<1>(qt_key)) / 2.0;
        double tavg = (std::get<2>(qt_key) + std::get<3>(qt_key)) / 2.0;
        graph->SetTitle(Form("%s: x_{B} = %.2f, Q^{2} = %.2f, -t = %.2f", type.c_str(), xB_avg, Q2avg, tavg));
        
        graph->SetMarkerStyle(20);
        graph->SetMarkerColor(kBlack);

        // Style the graph
        graph->GetXaxis()->SetTitle("#phi");
        graph->GetYaxis()->SetTitle("Unfolded Yields");
        graph->Draw("AP");

        plot_index++;
    }

    // Save the canvas
    canvas.SaveAs(Form("output/cross_check/RGAFa18Out/rga_fa18_out_cross_check_xB_%d.pdf", xB_index));
}

// Main function to control the import, processing, and debug output
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second) {
    
    // Ensure output directories exist
    ensure_directory_exists("output");
    ensure_directory_exists("output/cross_check");
    ensure_directory_exists("output/cross_check/RGAFa18Out");
    ensure_directory_exists("output/cross_check/RGAFa18Inb");


    // Step 1: Read the first CSV data into a vector of BinData
    std::vector<BinData> bin_data_first = read_csv_first(csv_file_path_first);
    // // Step 2: Print out the first CSV data to verify correctness
    // print_bin_data(bin_data_first);
    // Step 3: Read the second CSV data, matching bin numbers with the first CSV
    std::vector<BinData> bin_data_second = read_csv_second(csv_file_path_second, bin_data_first);
    // // Step 4: Print out the second CSV data to verify correctness
    // print_bin_data(bin_data_second);  

    auto unique_xB_bins = find_unique_xB_bins(bin_data_second);
    std::cout << "HELLO WORLD" << std::endl;

    int xB_index = 0;
    for (const auto &xB_range : unique_xB_bins) {
        auto filtered_data = filter_data_by_xB(bin_data_second, xB_range);
        plot_for_xB_bin(filtered_data, xB_index, "Out");
        xB_index++;
    }

}