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

    int index = 0;
    while (std::getline(file, line) && index < first_csv_data.size()) {
        std::stringstream ss(line);
        std::string value;

        BinData bin;
        bin.global_bin_number = first_csv_data[index].global_bin_number;
        bin.bin_number = first_csv_data[index].bin_number;

        std::getline(ss, value, ','); // Skip Bin number column
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

        for (int i = 0; i < 5; ++i) std::getline(ss, value, ',');
        std::getline(ss, value, ','); // ep->e'pgamma unfolded_yield_Fa18Inb
        bin.unfolded_yield_inbending = std::stod(value);

        for (int i = 0; i < 5; ++i) std::getline(ss, value, ',');
        std::getline(ss, value, ','); // ep->e'pgamma unfolded_yield_Fa18Out
        bin.unfolded_yield_outbending = std::stod(value);

        bins.push_back(bin);
        index++;
    }

    return bins;
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

// Plotting function for each xB bin with comparison between two datasets
void plot_for_xB_bin(const std::vector<BinData> &data_first, const std::vector<BinData> &data_second, int xB_index) {
    // Step 1: Identify unique (Q2min, Q2max, tmin, tmax) bins in the first dataset
    std::map<std::tuple<double, double, double, double>, std::vector<BinData>> qt_bins_first, qt_bins_second;

    for (const auto &bin : data_first) {
        auto key = std::make_tuple(bin.Q2min, bin.Q2max, bin.tmin, bin.tmax);
        qt_bins_first[key].push_back(bin);
    }

    for (const auto &bin : data_second) {
        auto key = std::make_tuple(bin.Q2min, bin.Q2max, bin.tmin, bin.tmax);
        qt_bins_second[key].push_back(bin);
    }

    // Step 2: Determine the grid size for subplots based on the first dataset's bins
    int num_plots = qt_bins_first.size();
    int grid_size = std::ceil(std::sqrt(num_plots));

    // For canvases _3 and _4, adjust grid size to have one less row if possible
    if ((xB_index == 3 || xB_index == 4) && grid_size * (grid_size - 1) >= num_plots) {
        grid_size -= 1;
    }

    TCanvas canvas("canvas", "Cross Check", 1200, 1200);
    canvas.Divide(grid_size, grid_size, 0.02, 0.02); // Small padding between pads

    int plot_index = 1;
    for (const auto &[qt_key, bins_first] : qt_bins_first) {
        canvas.cd(plot_index);
        gPad->SetLeftMargin(0.15);  // Adds padding to the left
        gPad->SetBottomMargin(0.15); // Adds padding to the bottom

        // Prepare vectors for the first dataset
        std::vector<double> phi_values_first, yields_first, yield_errors_first;

        for (const auto &bin : bins_first) {
            phi_values_first.push_back(bin.phiavg);
            yields_first.push_back(bin.unfolded_yield_outbending);
            yield_errors_first.push_back(std::sqrt(bin.unfolded_yield_outbending)); // Simple sqrt error
        }

        // Create TGraphErrors for the first dataset (blue)
        TGraphErrors *graph_first = new TGraphErrors(phi_values_first.size(), &phi_values_first[0], &yields_first[0], nullptr, &yield_errors_first[0]);
        graph_first->SetMarkerStyle(20);
        graph_first->SetMarkerColor(kBlue);

        // Prepare vectors for the second dataset (if matching (Q2, t) bin exists)
        auto it_second = qt_bins_second.find(qt_key);
        if (it_second != qt_bins_second.end()) {
            const auto &bins_second = it_second->second;
            std::vector<double> phi_values_second, yields_second, yield_errors_second;

            for (const auto &bin : bins_second) {
                phi_values_second.push_back(bin.phiavg);
                yields_second.push_back(bin.unfolded_yield_outbending);
                yield_errors_second.push_back(std::sqrt(bin.unfolded_yield_outbending)); // Simple sqrt error
            }

            // Create TGraphErrors for the second dataset (red)
            TGraphErrors *graph_second = new TGraphErrors(phi_values_second.size(), &phi_values_second[0], &yields_second[0], nullptr, &yield_errors_second[0]);
            graph_second->SetMarkerStyle(21);
            graph_second->SetMarkerColor(kRed);

            // Draw both graphs with legends for comparison
            graph_first->Draw("AP");
            graph_second->Draw("P SAME");
        } else {
            // If no matching (Q2, t) bin in the second dataset, draw only the first dataset graph
            graph_first->Draw("AP");
        }

        // Customize title and axis labels
        double xB_avg = bins_first[0].xBavg;
        double Q2avg = bins_first[0].Q2avg;
        double tavg = bins_first[0].tavg;
        graph_first->SetTitle(Form("Out: x_{B} = %.2f, Q^{2} = %.2f, -t = %.2f", xB_avg, Q2avg, tavg));

        graph_first->GetXaxis()->SetTitle("#phi");
        graph_first->GetYaxis()->SetTitle("Unfolded Yield");
        graph_first->GetXaxis()->SetLabelSize(0.04); // Increased font size
        graph_first->GetYaxis()->SetLabelSize(0.04); // Increased font size

        plot_index++;
    }

    // Save the canvas
    std::string save_path = Form("output/cross_check/RGAFa18Out/rga_fa18_out_cross_check_xB_%d.pdf", xB_index);
    canvas.SaveAs(save_path.c_str());
}

// Main function to control the import, processing, and debug output
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second) {
    ensure_directory_exists("output");
    ensure_directory_exists("output/cross_check");
    ensure_directory_exists("output/cross_check/RGAFa18Out");

    std::vector<BinData> bin_data_first = read_csv_first(csv_file_path_first);
    std::vector<BinData> bin_data_second = read_csv_second(csv_file_path_second, bin_data_first);

    auto unique_xB_bins = find_unique_xB_bins(bin_data_second);

    int xB_index = 0;
    for (const auto &xB_range : unique_xB_bins) {
        // Filter data for the current xB bin range in both datasets
        auto filtered_data_first = filter_data_by_xB(bin_data_first, xB_range);
        auto filtered_data_second = filter_data_by_xB(bin_data_second, xB_range);

        // Pass both datasets to the plotting function
        plot_for_xB_bin(filtered_data_second, filtered_data_second, xB_index);
        xB_index++;
    }
}