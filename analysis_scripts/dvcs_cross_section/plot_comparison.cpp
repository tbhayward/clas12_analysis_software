// plot_comparison.cpp

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
#include "utilities.h"


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

        // Skip intermediate columns until reaching "acceptance corrected yield, ep->epg, exp"
        for (int i = 0; i < 32; ++i) std::getline(ss, value, ',');

        std::getline(ss, value, ','); // Acceptance corrected yield
        double unfolded_yield = std::stod(value); // Store as a single total yield

        bin.signal_yields.push_back(unfolded_yield); // Only one yield in first CSV
        bin.signal_yield_uncertainties.push_back(0.0); // No uncertainties provided in first CSV

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

    size_t index = 0;
    while (std::getline(file, line) && index < first_csv_data.size()) {
        std::stringstream ss(line);
        std::string value;

        BinData bin;
        bin.global_bin_number = first_csv_data[index].global_bin_number;
        bin.bin_number = first_csv_data[index].bin_number;

        // Read initial binning information (13 columns)
        std::getline(ss, value, ','); // Bin number column (already read above)
        std::getline(ss, value, ','); // xB_min
        bin.xBmin = std::stod(value);
        std::getline(ss, value, ','); // xB_max
        bin.xBmax = std::stod(value);
        std::getline(ss, value, ','); // xB_avg
        bin.xBavg = std::stod(value);

        std::getline(ss, value, ','); // Q2_min
        bin.Q2min = std::stod(value);
        std::getline(ss, value, ','); // Q2_max
        bin.Q2max = std::stod(value);
        std::getline(ss, value, ','); // Q2_avg
        bin.Q2avg = std::stod(value);

        std::getline(ss, value, ','); // t_min
        bin.tmin = std::stod(value);
        std::getline(ss, value, ','); // t_max
        bin.tmax = std::stod(value);
        std::getline(ss, value, ','); // t_avg
        bin.tavg = std::stod(value);

        std::getline(ss, value, ','); // phi_min
        bin.phimin = std::stod(value);
        std::getline(ss, value, ','); // phi_max
        bin.phimax = std::stod(value);
        std::getline(ss, value, ','); // phi_avg
        bin.phiavg = std::stod(value);

        // Now, we need to skip columns to reach the signal yield for DVCS periods
        // We've read 13 columns so far

        // For each DVCS period (3 periods), we have:
        // 3 topology raw yields (FD_FD, CD_FD, CD_FT): 3 columns
        // 1 combined raw yield: 1 column
        // 1 acceptance: 1 column
        // 1 unfolded yield: 1 column
        // 1 unfolded yield uncertainty: 1 column
        // 1 contamination fraction: 1 column
        // 1 contamination uncertainty: 1 column
        // 1 signal yield: 1 column
        // 1 signal yield uncertainty: 1 column (newly added)
        // Total per period: 11 columns

        bin.signal_yields.clear();              // Ensure vector is empty
        bin.signal_yield_uncertainties.clear(); // Ensure vector is empty

        double total_signal_yield = 0.0;
        double total_signal_yield_uncertainty_squared = 0.0;

        const int n_periods = 2; // Only periods 0 and 1

        for (int period = 0; period < n_periods; ++period) {
            // Skip columns for this period up to signal yield
            // From current position, skip:
            // 3 (topology raw yields)
            // + 1 (combined raw yield)
            // + 1 (acceptance)
            // + 1 (unfolded yield)
            // + 1 (unfolded yield uncertainty)
            // + 1 (contamination fraction)
            // + 1 (contamination uncertainty)
            // So total of 9 columns to skip to reach signal_yield
            for (int i = 0; i < 9; ++i) std::getline(ss, value, ',');

            // Now read signal yield for this period
            std::getline(ss, value, ','); // Signal yield
            double signal_yield = std::stod(value);

            // Read signal yield uncertainty for this period
            std::getline(ss, value, ','); // Signal yield uncertainty
            double signal_yield_uncertainty = std::stod(value);

            // Sum the signal yields and uncertainties
            total_signal_yield += signal_yield;
            total_signal_yield_uncertainty_squared += signal_yield_uncertainty * signal_yield_uncertainty; // Sum in quadrature
        }

        // Skip columns for the third period (period 2)
        // Since each period has 11 columns, skip 11 columns
        for (int i = 0; i < 11; ++i) std::getline(ss, value, ',');

        // Store the total signal yield and uncertainty
        bin.signal_yields.push_back(total_signal_yield);
        bin.signal_yield_uncertainties.push_back(std::sqrt(total_signal_yield_uncertainty_squared));

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

// Plotting function for each xB bin with comparison between datasets
void plot_for_xB_bin(const std::vector<BinData> &data_first, const std::vector<BinData> &data_second, int xB_index) {
    // Step 1: Identify unique (Q2min, Q2max, tmin, tmax) bins in the first dataset
    std::map<std::tuple<double, double, double, double>, std::vector<BinData>> qt_bins_first, qt_bins_second;

    // Populate qt_bins_first and qt_bins_second
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

    // // Adjust grid size for canvas layout (for specified canvases _3 and _4)
    // if ((xB_index == 3 || xB_index == 4) && grid_size * (grid_size - 1) >= num_plots) {
    //     grid_size -= 1;
    // }

    // Step 3: Calculate the maximum yield value across both datasets for y-axis scaling
    double max_yield = 0.0;
    for (const auto &bin : data_first) {
        if (!bin.signal_yields.empty() && bin.signal_yields[0] > max_yield) {
            max_yield = bin.signal_yields[0];
        }
    }
    for (const auto &bin : data_second) {
        if (!bin.signal_yields.empty() && bin.signal_yields[0] > max_yield) {
            max_yield = bin.signal_yields[0];
        }
    }

    // Set the y-axis range
    double y_min = 0.1;
    double y_max = 1.1 * max_yield;

    TCanvas canvas("canvas", "Cross Check", 1200, 1200);
    canvas.Divide(grid_size, grid_size, 0.02, 0.02); // Small padding between pads

    int plot_index = 1;
    for (const auto &[qt_key, bins_first] : qt_bins_first) {
        canvas.cd(plot_index);
        gPad->SetLeftMargin(0.15);  // Adds padding to the left
        gPad->SetBottomMargin(0.15); // Adds padding to the bottom
        gPad->SetLogy();  // Set log scale for the y-axis

        // Prepare vectors for the first dataset
        std::vector<double> phi_values_first, yields_first, phi_errors_first, yield_uncertainties_first;
        for (const auto &bin : bins_first) {
            phi_values_first.push_back(bin.phiavg);
            if (!bin.signal_yields.empty()) {
                yields_first.push_back(bin.signal_yields[0]); // Use the only yield in first CSV
                yield_uncertainties_first.push_back(0.0);     // No uncertainties provided in first CSV
            } else {
                yields_first.push_back(0);
                yield_uncertainties_first.push_back(0.0);
            }
            double phi_bin_width = (bin.phimax - bin.phimin) / 2.0; // Calculate half-width of the phi bin
            phi_errors_first.push_back(phi_bin_width);
        }

        // Create TGraphErrors for the first dataset (blue)
        TGraphErrors *graph_first = new TGraphErrors(phi_values_first.size(),
                                                     &phi_values_first[0], &yields_first[0],
                                                     &phi_errors_first[0], &yield_uncertainties_first[0]);
        graph_first->SetMarkerStyle(20);
        graph_first->SetMarkerColor(kBlue);

        // Prepare and plot the second dataset if there’s a matching key
        auto it_second = qt_bins_second.find(qt_key);
        if (it_second != qt_bins_second.end()) {
            const auto &bins_second = it_second->second;

            std::vector<double> phi_values_second, yields_second, phi_errors_second, yield_uncertainties_second;

            for (const auto &bin : bins_second) {
                phi_values_second.push_back(bin.phiavg);
                if (!bin.signal_yields.empty()) {
                    yields_second.push_back(bin.signal_yields[0]); // Use the summed yield in second CSV
                    yield_uncertainties_second.push_back(bin.signal_yield_uncertainties[0]);
                } else {
                    yields_second.push_back(0);
                    yield_uncertainties_second.push_back(0);
                }
                double phi_bin_width = (bin.phimax - bin.phimin) / 2.0;
                phi_errors_second.push_back(phi_bin_width);
            }

            // Create TGraphErrors for the second dataset (red) with error bars
            TGraphErrors *graph_second = new TGraphErrors(phi_values_second.size(),
                                                          &phi_values_second[0], &yields_second[0],
                                                          &phi_errors_second[0], &yield_uncertainties_second[0]);
            graph_second->SetMarkerStyle(21);
            graph_second->SetMarkerColor(kRed);

            // Draw both graphs with legends for comparison
            graph_first->Draw("AP");
            graph_second->Draw("P SAME");

            // Create and draw legend at bottom right
            TLegend* legend = new TLegend(0.6, 0.15, 0.9, 0.35); // Adjusted legend position to bottom right
            legend->SetTextSize(0.03);
            legend->AddEntry(graph_first, "Lee, pass-1", "p");
            legend->AddEntry(graph_second, "Hayward, pass-2", "p");
            legend->Draw();
        } else {
            // If no matching (Q2, t) bin in the second dataset, draw only the first dataset graph
            graph_first->Draw("AP");
        }

        // Customize title and axis labels directly from data (no averaging needed)
        double xB_avg = bins_first[0].xBavg;
        double Q2avg = bins_first[0].Q2avg;
        double tavg = bins_first[0].tavg;
        graph_first->SetTitle(Form("x_{B} = %.2f, Q^{2} = %.2f, -t = %.2f", xB_avg, Q2avg, tavg));

        // Adjust axis labels and range
        graph_first->GetXaxis()->SetTitle("#phi [deg]");
        graph_first->GetYaxis()->SetTitle("Signal Yield");
        graph_first->GetXaxis()->SetLabelSize(0.04); // Increased font size
        graph_first->GetYaxis()->SetLabelSize(0.04); // Increased font size
        graph_first->GetYaxis()->SetRangeUser(y_min, y_max);  // Set y-axis range

        plot_index++;
    }

    // Save the canvas
    std::string save_path = Form("output/cross_check/RGAFa18/rga_fa18_cross_check_xB_%d.pdf", xB_index);
    ensure_directory_exists("output/cross_check/RGAFa18");
    canvas.SaveAs(save_path.c_str());
}

// Main function to control the import, processing, and debug output
void plot_comparison(const std::string &csv_file_path_first, const std::string &csv_file_path_second) {
    ensure_directory_exists("output");
    ensure_directory_exists("output/cross_check");

    std::vector<BinData> bin_data_first = read_csv_first(csv_file_path_first);
    std::vector<BinData> bin_data_second = read_csv_second(csv_file_path_second, bin_data_first);

    auto unique_xB_bins = find_unique_xB_bins(bin_data_second);

    int xB_index = 0;
    for (const auto &xB_range : unique_xB_bins) {
        // Filter data for the current xB bin range in both datasets
        auto filtered_data_first = filter_data_by_xB(bin_data_first, xB_range);
        auto filtered_data_second = filter_data_by_xB(bin_data_second, xB_range);

        // Pass both datasets to the plotting function
        plot_for_xB_bin(filtered_data_first, filtered_data_second, xB_index);
        xB_index++;
    }
}