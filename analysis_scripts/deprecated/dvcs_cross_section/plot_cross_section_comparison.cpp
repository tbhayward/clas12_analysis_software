// plot_cross_section_comparison.cpp

#include "plot_cross_section_comparison.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <map>
#include <tuple>
#include <algorithm>
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
#include <TLatex.h>
#include <TGraphAsymmErrors.h>
#include <TH1F.h>  // Added header for TH1F

#include "utilities.h" // Include the utilities header

// Structure to hold bin data
struct CrossSectionBinData {
    int bin_number;
    double xB_min, xB_max, xB_avg;
    double Q2_min, Q2_max, Q2_avg;
    double t_min, t_max, t_avg;
    double phi_min, phi_max, phi_avg;
    double cross_section;
    double stat_uncertainty;
    double sys_uncertainty;
};

// Type alias for a vector of CrossSectionBinData
typedef std::vector<CrossSectionBinData> CrossSectionData;

// Helper function to read data from the first CSV file format
CrossSectionData read_first_csv(const std::string &filename) {
    CrossSectionData data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;
    // Skip header
    std::getline(file, line);

    // Read data lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string value;

        CrossSectionBinData entry;

        // Read columns by index
        // Column 0: Empty
        std::getline(ss, value, ','); // Empty column

        // Column 1: Bin Name
        std::getline(ss, value, ',');
        entry.bin_number = std::stoi(value);

        // Columns 2-4: xB_min, xB_max, xB_avg
        std::getline(ss, value, ',');
        entry.xB_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.xB_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.xB_avg = std::stod(value);

        // Columns 5-7: Q2_min, Q2_max, Q2_avg
        std::getline(ss, value, ',');
        entry.Q2_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.Q2_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.Q2_avg = std::stod(value);

        // Columns 8-10: t_min, t_max, t_avg
        std::getline(ss, value, ',');
        entry.t_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.t_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.t_avg = std::stod(value);

        // Columns 11-13: phi_min, phi_max, phi_avg
        std::getline(ss, value, ',');
        entry.phi_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.phi_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.phi_avg = std::stod(value);

        // Skip intermediate columns until reaching "cross sections, ep->epg, exp"
        // Currently at column index 13, need to skip 41 columns to reach column 55
        for (int i = 0; i < 41; ++i) std::getline(ss, value, ',');

        // Column 55: cross sections, ep->epg, exp
        std::getline(ss, value, ',');
        entry.cross_section = std::stod(value);

        // Skip 3 columns to reach "cross sections, ep->epg, exp, stat. unc."
        for (int i = 0; i < 3; ++i) std::getline(ss, value, ',');

        // Column 59: cross sections, ep->epg, exp, stat. unc.
        std::getline(ss, value, ',');
        entry.stat_uncertainty = std::stod(value);

        // Column 60: cross sections, ep->epg, exp, syst. unc. (up)
        std::getline(ss, value, ',');
        entry.sys_uncertainty = std::stod(value);

        // Skip one column to reach "valid bin"
        std::getline(ss, value, ',');

        // Column: valid bin
        std::getline(ss, value, ',');
        int valid_bin = std::stoi(value);

        // Check if valid_bin is 0, skip the bin if so
        if (valid_bin == 0) {
            continue; // Skip this bin
        }

        data.push_back(entry);
    }

    file.close();
    return data;
}

// Function to read the second CSV file
CrossSectionData read_second_csv(const std::string &filename) {
    CrossSectionData data;
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;
    // Skip header
    std::getline(file, line);

    // Read data lines
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::istringstream ss(line);
        std::string value;

        CrossSectionBinData entry;

        // Read columns by index
        // Column 0: Bin
        std::getline(ss, value, ',');
        entry.bin_number = std::stoi(value);

        // Columns 1-3: xB_min, xB_max, xB_avg
        std::getline(ss, value, ',');
        entry.xB_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.xB_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.xB_avg = std::stod(value);

        // Columns 4-6: Q2_min, Q2_max, Q2_avg
        std::getline(ss, value, ',');
        entry.Q2_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.Q2_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.Q2_avg = std::stod(value);

        // Columns 7-9: t_min, t_max, t_avg
        std::getline(ss, value, ',');
        entry.t_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.t_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.t_avg = std::stod(value);

        // Columns 10-12: phi_min, phi_max, phi_avg
        std::getline(ss, value, ',');
        entry.phi_min = std::stod(value);
        std::getline(ss, value, ',');
        entry.phi_max = std::stod(value);
        std::getline(ss, value, ',');
        entry.phi_avg = std::stod(value);

        // Skip columns until reaching "fall_cross_section" (Column index 45)
        // Currently at column 12, need to skip 32 columns
        for (int i = 0; i < 55; ++i) std::getline(ss, value, ',');

        // Column 45: fall_cross_section
        std::getline(ss, value, ',');
        entry.cross_section = std::stod(value);

        // Column 46: fall_cross_section_stat_uncertainty
        std::getline(ss, value, ',');
        entry.stat_uncertainty = std::stod(value);

        // Column 47: fall_cross_section_sys_uncertainty
        std::getline(ss, value, ',');
        // entry.sys_uncertainty = std::stod(value);
        // entry.sys_uncertainty = 0;

        entry.cross_section=0.98*entry.cross_section/0.80/0.93;
        entry.stat_uncertainty=1.5*0.98*entry.stat_uncertainty/0.80/0.93;
        entry.sys_uncertainty=0.204*entry.cross_section*3;

        data.push_back(entry);
    }

    file.close();
    return data;
}

// Function to find unique xB bins in the data
std::vector<std::pair<double, double>> find_unique_xB_bins(const CrossSectionData &data) {
    std::set<std::pair<double, double>> unique_bins;
    for (const auto &bin : data) {
        unique_bins.emplace(bin.xB_min, bin.xB_max);
    }
    return std::vector<std::pair<double, double>>(unique_bins.begin(), unique_bins.end());
}

// Function to filter data for a specific xB bin
CrossSectionData filter_data_by_xB(const CrossSectionData &data, const std::pair<double, double> &xB_range) {
    CrossSectionData filtered_data;
    for (const auto &bin : data) {
        if (bin.xB_min == xB_range.first && bin.xB_max == xB_range.second) {
            filtered_data.push_back(bin);
        }
    }
    return filtered_data;
}

// Plotting function for each xB bin with comparison between datasets
void plot_for_xB_bin(const CrossSectionData &data_first, const CrossSectionData &data_second, int xB_index) {
    // Step 1: Identify unique (Q2min, Q2max, tmin, tmax) bins in the first dataset
    std::map<std::tuple<double, double, double, double>, CrossSectionData> qt_bins_first, qt_bins_second;

    // Populate qt_bins_first and qt_bins_second
    for (const auto &bin : data_first) {
        auto key = std::make_tuple(bin.Q2_min, bin.Q2_max, bin.t_min, bin.t_max);
        qt_bins_first[key].push_back(bin);
    }
    for (const auto &bin : data_second) {
        auto key = std::make_tuple(bin.Q2_min, bin.Q2_max, bin.t_min, bin.t_max);
        qt_bins_second[key].push_back(bin);
    }

    // Collect all cross-section values (excluding zeros) to determine global y-axis limits
    double min_cs = std::numeric_limits<double>::max();
    double max_cs = std::numeric_limits<double>::lowest();

    // Iterate over all bins in both datasets to find global min and max
    for (const auto &bins_first : qt_bins_first) {
        for (const auto &bin : bins_first.second) {
            if (bin.cross_section > 0) {
                double cs_min = bin.cross_section;
                double cs_max = bin.cross_section;
                if (cs_min < min_cs && cs_min > 0) min_cs = cs_min;
                if (cs_max > max_cs) max_cs = cs_max;
            }
        }
    }
    for (const auto &bins_second : qt_bins_second) {
        for (const auto &bin : bins_second.second) {
            if (bin.cross_section > 0) {
                double cs_min = bin.cross_section;
                double cs_max = bin.cross_section;
                if (cs_min < min_cs && cs_min > 0) min_cs = cs_min;
                if (cs_max > max_cs) max_cs = cs_max;
            }
        }
    }

    // Adjust min_cs and max_cs to be positive and exclude zeros
    if (min_cs <= 0 || min_cs == std::numeric_limits<double>::max()) min_cs = 1e-4;
    if (max_cs <= 0 || max_cs == std::numeric_limits<double>::lowest()) max_cs = 1.0;

    // Calculate y_min and y_max as powers of 10
    double log_min_cs = std::floor(std::log10(min_cs));
    double log_max_cs = std::ceil(std::log10(max_cs));

    double y_min = std::pow(10, log_min_cs);
    double y_max = std::pow(10, log_max_cs);

    // Calculate grid_size before using it
    int num_plots = qt_bins_first.size();
    int grid_size = std::ceil(std::sqrt(num_plots));

    // Adjust grid size for canvas layout (if needed)
    if ((xB_index == 3 || xB_index == 4) && grid_size * (grid_size - 1) >= num_plots) {
        grid_size -= 1;
    }

    // Step 3: Set up the canvas
    TCanvas canvas("canvas", "Cross Section Cross Check", 1200, 1200);
    canvas.Divide(grid_size, grid_size, 0.02, 0.02); // Small padding between pads

    int plot_index = 1;
    for (const auto &[qt_key, bins_first] : qt_bins_first) {
        canvas.cd(plot_index);
        gPad->SetLeftMargin(0.15);  // Adds padding to the left
        gPad->SetBottomMargin(0.15); // Adds padding to the bottom
        gPad->SetLogy();  // Set log scale for the y-axis

        // Prepare data for the first dataset
        struct DataPoint {
            double phi;
            double cs;
            double phi_err;
            double stat_err;
            double sys_err;
            double phi_min;
            double phi_max;
        };
        std::vector<DataPoint> data_points_first;
        for (const auto &bin : bins_first) {
            // Skip data points with zero or negative cross-section
            if (bin.cross_section <= 0) continue;

            DataPoint dp;
            dp.phi = bin.phi_avg;
            dp.cs = bin.cross_section;
            dp.phi_err = (bin.phi_max - bin.phi_min) / 2.0;
            dp.stat_err = bin.stat_uncertainty;
            dp.sys_err = bin.sys_uncertainty;
            dp.phi_min = bin.phi_min;
            dp.phi_max = bin.phi_max;
            data_points_first.push_back(dp);
        }
        // Sort data points by phi
        std::sort(data_points_first.begin(), data_points_first.end(), [](const DataPoint &a, const DataPoint &b) {
            return a.phi < b.phi;
        });

        // Prepare data for the second dataset if available
        auto it_second = qt_bins_second.find(qt_key);
        std::vector<DataPoint> data_points_second;
        if (it_second != qt_bins_second.end()) {
            const auto &bins_second = it_second->second;
            for (const auto &bin : bins_second) {
                // Skip data points with zero or negative cross-section
                if (bin.cross_section <= 0) continue;

                DataPoint dp;
                dp.phi = bin.phi_avg;
                dp.cs = bin.cross_section;
                dp.phi_err = (bin.phi_max - bin.phi_min) / 2.0;
                dp.stat_err = bin.stat_uncertainty;
                dp.sys_err = bin.sys_uncertainty;
                dp.phi_min = bin.phi_min;
                dp.phi_max = bin.phi_max;
                data_points_second.push_back(dp);
            }
            // Sort data points by phi
            std::sort(data_points_second.begin(), data_points_second.end(), [](const DataPoint &a, const DataPoint &b) {
                return a.phi < b.phi;
            });
        }

        // Create a frame for the plot with global y-axis limits
        TH1F *frame = gPad->DrawFrame(0, y_min, 360, y_max);
        frame->GetXaxis()->SetTitle("#phi [deg]");
        frame->GetYaxis()->SetTitle("d#sigma/dx_{B}dQ^{2}d|t|d#phi (nb/GeV^{4})");

        // Declare graph pointers before the if blocks
        TGraphErrors *graph_first_stat = nullptr;
        TGraphErrors *graph_second_stat = nullptr;

        // Proceed only if there are data points to plot
        if (!data_points_first.empty() || !data_points_second.empty()) {
            // Plot first dataset if it has data points
            if (!data_points_first.empty()) {
                int n_points_first = data_points_first.size();
                std::vector<double> phi_values_first(n_points_first);
                std::vector<double> cs_first(n_points_first);
                std::vector<double> phi_errors_first(n_points_first);
                std::vector<double> stat_err_first(n_points_first);
                std::vector<double> sys_err_first(n_points_first);
                std::vector<double> phi_min_first(n_points_first);
                std::vector<double> phi_max_first(n_points_first);

                for (int i = 0; i < n_points_first; ++i) {
                    phi_values_first[i] = data_points_first[i].phi;
                    cs_first[i] = data_points_first[i].cs;
                    phi_errors_first[i] = data_points_first[i].phi_err;
                    stat_err_first[i] = data_points_first[i].stat_err;
                    sys_err_first[i] = data_points_first[i].sys_err;
                    phi_min_first[i] = data_points_first[i].phi_min;
                    phi_max_first[i] = data_points_first[i].phi_max;
                }

                // Create TGraphErrors for statistical uncertainties (first dataset)
                graph_first_stat = new TGraphErrors(n_points_first,
                                                    &phi_values_first[0], &cs_first[0],
                                                    &phi_errors_first[0], &stat_err_first[0]);
                graph_first_stat->SetMarkerStyle(20);
                graph_first_stat->SetMarkerColor(kBlue);
                graph_first_stat->SetLineColor(kBlue);

                // Draw the systematic uncertainty bands as TBox for each bin
                for (int i = 0; i < n_points_first; ++i) {
                    double x_min = phi_min_first[i];
                    double x_max = phi_max_first[i];
                    double y_min = cs_first[i] - sys_err_first[i];
                    double y_max = cs_first[i] + sys_err_first[i];
                    TBox *box = new TBox(x_min, y_min, x_max, y_max);
                    box->SetFillColorAlpha(kBlue, 0.35);
                    box->SetLineColor(kBlue);
                    box->Draw("same");
                }

                // Draw the statistical uncertainties
                graph_first_stat->Draw("P SAME");
            }

            // Plot second dataset if it has data points
            if (!data_points_second.empty()) {
                int n_points_second = data_points_second.size();
                std::vector<double> phi_values_second(n_points_second);
                std::vector<double> cs_second(n_points_second);
                std::vector<double> phi_errors_second(n_points_second);
                std::vector<double> stat_err_second(n_points_second);
                std::vector<double> sys_err_second(n_points_second);
                std::vector<double> phi_min_second(n_points_second);
                std::vector<double> phi_max_second(n_points_second);

                for (int i = 0; i < n_points_second; ++i) {
                    phi_values_second[i] = data_points_second[i].phi;
                    cs_second[i] = data_points_second[i].cs;
                    phi_errors_second[i] = data_points_second[i].phi_err;
                    stat_err_second[i] = data_points_second[i].stat_err;
                    sys_err_second[i] = data_points_second[i].sys_err;
                    phi_min_second[i] = data_points_second[i].phi_min;
                    phi_max_second[i] = data_points_second[i].phi_max;
                }

                // Create TGraphErrors for statistical uncertainties (second dataset)
                graph_second_stat = new TGraphErrors(n_points_second,
                                                     &phi_values_second[0], &cs_second[0],
                                                     &phi_errors_second[0], &stat_err_second[0]);
                graph_second_stat->SetMarkerStyle(21);
                graph_second_stat->SetMarkerColor(kRed);
                graph_second_stat->SetLineColor(kRed);

                // Draw the systematic uncertainty bands as TBox for each bin
                for (int i = 0; i < n_points_second; ++i) {
                    double x_min = phi_min_second[i];
                    double x_max = phi_max_second[i];
                    double y_min = cs_second[i] - sys_err_second[i];
                    double y_max = cs_second[i] + sys_err_second[i];
                    TBox *box = new TBox(x_min, y_min, x_max, y_max);
                    box->SetFillColorAlpha(kRed, 0.35);
                    box->SetLineColor(kRed);
                    box->Draw("same");
                }

                // Draw the statistical uncertainties
                graph_second_stat->Draw("P SAME");
            }

            // Create and draw legend
            TLegend* legend = new TLegend(0.55, 0.65, 0.9, 0.85);
            legend->SetTextSize(0.03);
            if (graph_first_stat) {
                legend->AddEntry(graph_first_stat, "pass-1, Lee", "lep");
            }
            if (graph_second_stat) {
                legend->AddEntry(graph_second_stat, "pass-2, Hayward", "lep");
            }
            legend->Draw();
        } else {
            // If no data points to plot, display a message
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.05);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.5, "No data to display");
        }

        // Add title with bin information
        double xB_avg = bins_first[0].xB_avg;
        double Q2_avg = bins_first[0].Q2_avg;
        double t_avg = bins_first[0].t_avg;
        std::string title = Form("x_{B} = %.2f, Q^{2} = %.2f, -t = %.2f", xB_avg, Q2_avg, t_avg);
        TLatex latex_title;
        latex_title.SetNDC();
        latex_title.SetTextSize(0.04);
        latex_title.SetTextAlign(22); // Center alignment
        latex_title.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

        plot_index++;
    }

    // Save the canvas
    std::string save_path = Form("output/cross_section_cross_check/cross_section_cross_check_xB_%d.pdf", xB_index);
    ensure_directory_exists("output/cross_section_cross_check/");
    canvas.SaveAs(save_path.c_str());
}

// Main function to control the import, processing, and plotting
void plot_cross_section_comparison(const std::string &first_csv_file,
                                   const std::string &second_csv_file) {
    ensure_directory_exists("output/cross_section_cross_check/");

    // Read data from CSV files
    CrossSectionData data_first_csv = read_first_csv(first_csv_file);
    CrossSectionData data_second_csv = read_second_csv(second_csv_file);

    // Find unique xB bins
    auto unique_xB_bins = find_unique_xB_bins(data_first_csv);

    int xB_index = 0;
    for (const auto &xB_range : unique_xB_bins) {
        // Filter data for the current xB bin range in both datasets
        auto filtered_data_first = filter_data_by_xB(data_first_csv, xB_range);
        auto filtered_data_second = filter_data_by_xB(data_second_csv, xB_range);

        // Pass both datasets to the plotting function
        plot_for_xB_bin(filtered_data_first, filtered_data_second, xB_index);
        xB_index++;
    }
}