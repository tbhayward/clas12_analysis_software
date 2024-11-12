// plot_unfolding.cpp

#include <TH1D.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <TSystem.h>   // For gSystem->mkdir
#include <cmath>
#include <string>
#include <vector>
#include <map>
#include <iostream>
#include <algorithm>
#include <numeric>

#include "bin_boundaries.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include "unfolding_data.h"

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Helper function to check if phi is within bin, accounting for wrap-around
static bool phi_in_bin(double phi_deg, double phi_low, double phi_high) {
    if (phi_low <= phi_high) {
        return phi_deg >= phi_low && phi_deg <= phi_high;
    } else {
        // Wrap-around bin (e.g., from 330 to 30 degrees)
        return phi_deg >= phi_low || phi_deg <= phi_high;
    }
}

// Function implementation
std::map<std::string, std::vector<UnfoldingData>> plot_unfolding(
    const std::string& output_dir,
    int xB_bin,
    const std::vector<BinBoundary>& bin_boundaries,
    std::vector<TTreeReader>& data_readers,
    std::vector<TTreeReader>& mc_gen_readers,
    std::vector<TTreeReader>& mc_rec_readers,
    std::vector<TTreeReader>& eppi0_readers,
    std::vector<TTreeReader>& mc_gen_aaogen_readers,
    std::vector<TTreeReader>& mc_rec_aaogen_readers
) {
    // Define topologies
    std::vector<std::string> topologies = {"CD_FD", "CD_FT", "FD_FD"};
    std::string combined_topology = "combined";

    // Names of the periods
    std::vector<std::string> period_names = {"Fa18Inb", "Fa18Out", "Sp19Inb"};
    const int n_periods = 3; // Number of periods for data and MC

    // Initialize a map to hold all unfolding data for each topology
    std::map<std::string, std::vector<UnfoldingData>> topology_unfolding_data;

    // Precompute relevant bins for the given xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);
    if (relevant_bins.empty()) {
        std::cerr << "Error: No relevant bins found for xB_bin: " << xB_bin << std::endl;
        return topology_unfolding_data; // Return empty map
    }

    // Group the bins by (Q2_low, t_low)
    std::map<std::pair<double, double>, std::vector<int>> bin_groups;
    for (int idx : relevant_bins) {
        const auto& bin = bin_boundaries[idx];
        auto key = std::make_pair(bin.Q2_low, bin.t_low);
        bin_groups[key].push_back(idx);
    }

    if (bin_groups.empty()) {
        std::cerr << "Error: No bin groups created. Exiting function." << std::endl;
        return topology_unfolding_data; // Return empty map
    }

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    gStyle->SetOptStat(0);

    // Map to hold combined unfolding data
    std::vector<UnfoldingData> combined_unfolding_data;

    // For each topology
    for (const auto& topology : topologies) {
        // Vector to hold all unfolding data for this topology
        std::vector<UnfoldingData> all_unfolding_data;

        int n_bins = bin_groups.size();
        int n_columns = std::ceil(std::sqrt(n_bins));
        int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

        // Initialize UnfoldingData structures and histograms
        for (const auto& group : bin_groups) {
            const auto& idx_list = group.second; // List of indices in bin_boundaries

            // Create an UnfoldingData instance
            UnfoldingData unfolding_data;

            // Initialize bin ranges and averages
            const auto& bin_example = bin_boundaries[idx_list[0]];
            unfolding_data.bin_number = xB_bin;
            unfolding_data.xB_min = bin_example.xB_low;
            unfolding_data.xB_max = bin_example.xB_high;
            unfolding_data.xB_avg = bin_example.xB_avg;
            unfolding_data.Q2_min = bin_example.Q2_low;
            unfolding_data.Q2_max = bin_example.Q2_high;
            unfolding_data.Q2_avg = bin_example.Q2_avg;
            unfolding_data.t_min = bin_example.t_low;
            unfolding_data.t_max = bin_example.t_high;
            unfolding_data.t_avg = bin_example.t_avg;

            // Collect phi bins
            for (int idx : idx_list) {
                const auto& bin = bin_boundaries[idx];
                unfolding_data.phi_min.push_back(bin.phi_low);
                unfolding_data.phi_max.push_back(bin.phi_high);
            }

            size_t n_phi_bins = unfolding_data.phi_min.size();
            size_t total_periods = 6; // 0-2: DVCS, 3-5: eppi0

            // Initialize counts and acceptance vectors
            unfolding_data.raw_yields.resize(total_periods, std::vector<int>(n_phi_bins, 0));
            unfolding_data.acceptance.resize(total_periods, std::vector<double>(n_phi_bins, 0.0));
            unfolding_data.unfolded_yields.resize(total_periods, std::vector<double>(n_phi_bins, 0.0));

            // Store the unfolding_data instance
            all_unfolding_data.push_back(unfolding_data);
        }

        // Initialize local variables to store mc_gen_counts and mc_rec_counts
        size_t total_periods = 6; // 0-2: DVCS, 3-5: eppi0
        size_t n_groups = all_unfolding_data.size();
        // mc_gen_counts[period][group_idx][phi_idx]
        std::vector<std::vector<std::vector<double>>> mc_gen_counts(total_periods, std::vector<std::vector<double>>(n_groups));
        std::vector<std::vector<std::vector<double>>> mc_rec_counts(total_periods, std::vector<std::vector<double>>(n_groups));

        // Initialize mc_gen_counts and mc_rec_counts vectors
        for (size_t period = 0; period < total_periods; ++period) {
            for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
                size_t n_phi_bins = all_unfolding_data[group_idx].phi_min.size();
                mc_gen_counts[period][group_idx].resize(n_phi_bins, 0.0);
                mc_rec_counts[period][group_idx].resize(n_phi_bins, 0.0);
            }
        }

        // The rest of the code remains the same for processing data and MC
        // (This includes processing data readers, eppi0 readers, mc_gen_readers, mc_rec_readers, etc.)
        // ...

        // [Process data readers, eppi0 readers, mc_gen_readers, mc_rec_readers as in your existing code]

        // Compute acceptances and unfolded yields
        // [Compute acceptances and unfolded yields as in your existing code]

        // Add the unfolding data for this topology to the combined data
        if (combined_unfolding_data.empty()) {
            combined_unfolding_data = all_unfolding_data;
        } else {
            // Sum the data across topologies
            for (size_t group_idx = 0; group_idx < n_groups; ++group_idx) {
                UnfoldingData& combined_data = combined_unfolding_data[group_idx];
                UnfoldingData& current_data = all_unfolding_data[group_idx];

                for (size_t period = 0; period < total_periods; ++period) {
                    for (size_t phi_idx = 0; phi_idx < combined_data.phi_min.size(); ++phi_idx) {
                        combined_data.raw_yields[period][phi_idx] += current_data.raw_yields[period][phi_idx];
                        combined_data.unfolded_yields[period][phi_idx] += current_data.unfolded_yields[period][phi_idx];
                        combined_data.acceptance[period][phi_idx] = (combined_data.acceptance[period][phi_idx] + current_data.acceptance[period][phi_idx]) / 2.0;
                    }
                }
            }
        }

        // Plotting code for DVCS data (periods 0-2)
        for (int period = 0; period < n_periods; ++period) {
            // Create output directories
            std::string analysis_type = "dvcs";
            std::string output_subdir = output_dir + "/unfolded/" + analysis_type + "/" + period_names[period];
            gSystem->mkdir(output_subdir.c_str(), true);

            size_t n_bins = all_unfolding_data.size();
            int n_columns = std::ceil(std::sqrt(n_bins));
            int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

            // Create the canvas outside the loop
            TCanvas* canvas = new TCanvas(Form("c_unfolding_dvcs_period_%d_%s", period, topology.c_str()),
                                          "Unfolding Results DVCS", canvas_width, canvas_height);
            canvas->Divide(n_columns, n_rows);

            for (size_t group_idx = 0; group_idx < all_unfolding_data.size(); ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];
                size_t n_phi_bins = unfolding_data.phi_min.size();

                // Prepare data for plotting
                std::vector<double> phi_centers(n_phi_bins);
                std::vector<double> phi_widths(n_phi_bins);
                std::vector<double> unfolded_counts(n_phi_bins);
                std::vector<double> unfolded_errors(n_phi_bins);

                for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                    phi_centers[phi_idx] = (unfolding_data.phi_min[phi_idx] + unfolding_data.phi_max[phi_idx]) / 2.0;
                    phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                    unfolded_counts[phi_idx] = unfolding_data.unfolded_yields[period][phi_idx];
                    unfolded_errors[phi_idx] = std::sqrt(unfolding_data.raw_yields[period][phi_idx]) /
                                               (unfolding_data.acceptance[period][phi_idx] > 0 ? unfolding_data.acceptance[period][phi_idx] : 1);
                }

                // Move to the appropriate pad
                canvas->cd(group_idx + 1);

                // Create TGraphErrors
                TGraphErrors* graph_unfolded = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts[0], &phi_widths[0], &unfolded_errors[0]);

                // Set styles
                graph_unfolded->SetMarkerColor(kBlack);
                graph_unfolded->SetMarkerStyle(20);
                graph_unfolded->SetLineColor(kBlack);

                // Draw graph
                TPad* pad = (TPad*)gPad;
                pad->SetLeftMargin(0.15);
                pad->SetBottomMargin(0.15);

                // Create a frame histogram for axes
                double phi_min = 0.0;
                double phi_max = 360.0;
                double count_max = *std::max_element(unfolded_counts.begin(), unfolded_counts.end()) * 1.2;
                double count_min = 0.0;

                TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
                frame->GetXaxis()->SetTitle("#phi [deg]");
                frame->GetYaxis()->SetTitle("Unfolded Counts");

                graph_unfolded->Draw("P SAME");

                // Add title with averages and topology
                std::string title = Form("DVCS, %s, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                         period_names[period].c_str(),
                                         topology.c_str(),
                                         unfolding_data.xB_avg,
                                         unfolding_data.Q2_avg,
                                         unfolding_data.t_avg);
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.04);
                latex.SetTextAlign(22); // Center alignment
                latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

                // Clean up graph (if needed)
                // delete graph_unfolded;
            }

            // Save the canvas after all subplots have been added
            std::string filename = output_subdir + "/unfolding_results_DVCS_" + period_names[period] +
                                   "_xB_bin_" + std::to_string(xB_bin) + "_" + topology + ".pdf";
            canvas->SaveAs(filename.c_str());

            // Clean up
            delete canvas;
        }

        // Plotting code for eppi0 data (periods 3-5)
        for (int period = 0; period < n_periods; ++period) {
            int eppi0_period = period + 3; // Adjust the period index for eppi0 data

            // Create output directories
            std::string analysis_type = "eppi0";
            std::string output_subdir = output_dir + "/unfolded/" + analysis_type + "/" + period_names[period];
            gSystem->mkdir(output_subdir.c_str(), true);

            size_t n_bins = all_unfolding_data.size();
            int n_columns = std::ceil(std::sqrt(n_bins));
            int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

            // Create the canvas outside the loop
            TCanvas* canvas = new TCanvas(Form("c_unfolding_eppi0_period_%d_%s", eppi0_period, topology.c_str()),
                                          "Unfolding Results eppi0", canvas_width, canvas_height);
            canvas->Divide(n_columns, n_rows);

            for (size_t group_idx = 0; group_idx < all_unfolding_data.size(); ++group_idx) {
                UnfoldingData& unfolding_data = all_unfolding_data[group_idx];
                size_t n_phi_bins = unfolding_data.phi_min.size();

                // Prepare data for plotting
                std::vector<double> phi_centers(n_phi_bins);
                std::vector<double> phi_widths(n_phi_bins);
                std::vector<double> unfolded_counts(n_phi_bins);
                std::vector<double> unfolded_errors(n_phi_bins);

                for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                    phi_centers[phi_idx] = (unfolding_data.phi_min[phi_idx] + unfolding_data.phi_max[phi_idx]) / 2.0;
                    phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                    unfolded_counts[phi_idx] = unfolding_data.unfolded_yields[eppi0_period][phi_idx];
                    unfolded_errors[phi_idx] = std::sqrt(unfolding_data.raw_yields[eppi0_period][phi_idx]) /
                                               (unfolding_data.acceptance[eppi0_period][phi_idx] > 0 ? unfolding_data.acceptance[eppi0_period][phi_idx] : 1);
                }

                // Move to the appropriate pad
                canvas->cd(group_idx + 1);

                // Create TGraphErrors
                TGraphErrors* graph_unfolded = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts[0], &phi_widths[0], &unfolded_errors[0]);

                // Set styles
                graph_unfolded->SetMarkerColor(kBlack);
                graph_unfolded->SetMarkerStyle(20);
                graph_unfolded->SetLineColor(kBlack);

                // Draw graph
                TPad* pad = (TPad*)gPad;
                pad->SetLeftMargin(0.15);
                pad->SetBottomMargin(0.15);

                // Create a frame histogram for axes
                double phi_min = 0.0;
                double phi_max = 360.0;
                double count_max = *std::max_element(unfolded_counts.begin(), unfolded_counts.end()) * 1.2;
                double count_min = 0.0;

                TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
                frame->GetXaxis()->SetTitle("#phi [deg]");
                frame->GetYaxis()->SetTitle("Unfolded Counts");

                graph_unfolded->Draw("P SAME");

                // Add title with averages and topology
                std::string title = Form("eppi0, %s, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                         period_names[period].c_str(),
                                         topology.c_str(),
                                         unfolding_data.xB_avg,
                                         unfolding_data.Q2_avg,
                                         unfolding_data.t_avg);
                TLatex latex;
                latex.SetNDC();
                latex.SetTextSize(0.04);
                latex.SetTextAlign(22); // Center alignment
                latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5

                // Clean up graph (if needed)
                // delete graph_unfolded;
            }

            // Save the canvas after all subplots have been added
            std::string filename = output_subdir + "/unfolding_results_eppi0_" + period_names[period] +
                                   "_xB_bin_" + std::to_string(xB_bin) + "_" + topology + ".pdf";
            canvas->SaveAs(filename.c_str());

            // Clean up
            delete canvas;
        }

        // Store the unfolding data for this topology
        topology_unfolding_data[topology] = all_unfolding_data;
    }

    // Add the combined topology data
    topology_unfolding_data[combined_topology] = combined_unfolding_data;

    // Plotting code for the combined topology
    // DVCS data (periods 0-2)
    for (int period = 0; period < n_periods; ++period) {
        // Create output directories
        std::string analysis_type = "dvcs";
        std::string output_subdir = output_dir + "/unfolded/" + analysis_type + "/" + period_names[period];
        gSystem->mkdir(output_subdir.c_str(), true);

        size_t n_bins = combined_unfolding_data.size();
        int n_columns = std::ceil(std::sqrt(n_bins));
        int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

        // Create the canvas
        TCanvas* canvas = new TCanvas(Form("c_unfolding_dvcs_period_%d_%s", period, combined_topology.c_str()),
                                      "Unfolding Results DVCS Combined", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        for (size_t group_idx = 0; group_idx < combined_unfolding_data.size(); ++group_idx) {
            UnfoldingData& unfolding_data = combined_unfolding_data[group_idx];
            size_t n_phi_bins = unfolding_data.phi_min.size();

            // Prepare data for plotting
            std::vector<double> phi_centers(n_phi_bins);
            std::vector<double> phi_widths(n_phi_bins);
            std::vector<double> unfolded_counts(n_phi_bins);
            std::vector<double> unfolded_errors(n_phi_bins);

            for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                phi_centers[phi_idx] = (unfolding_data.phi_min[phi_idx] + unfolding_data.phi_max[phi_idx]) / 2.0;
                phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                unfolded_counts[phi_idx] = unfolding_data.unfolded_yields[period][phi_idx];
                unfolded_errors[phi_idx] = std::sqrt(unfolding_data.raw_yields[period][phi_idx]) /
                                           (unfolding_data.acceptance[period][phi_idx] > 0 ? unfolding_data.acceptance[period][phi_idx] : 1);
            }

            // Move to the appropriate pad
            canvas->cd(group_idx + 1);

            // Create TGraphErrors
            TGraphErrors* graph_unfolded = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts[0], &phi_widths[0], &unfolded_errors[0]);

            // Set styles
            graph_unfolded->SetMarkerColor(kBlack);
            graph_unfolded->SetMarkerStyle(20);
            graph_unfolded->SetLineColor(kBlack);

            // Draw graph
            TPad* pad = (TPad*)gPad;
            pad->SetLeftMargin(0.15);
            pad->SetBottomMargin(0.15);

            // Create a frame histogram for axes
            double phi_min = 0.0;
            double phi_max = 360.0;
            double count_max = *std::max_element(unfolded_counts.begin(), unfolded_counts.end()) * 1.2;
            double count_min = 0.0;

            TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
            frame->GetXaxis()->SetTitle("#phi [deg]");
            frame->GetYaxis()->SetTitle("Unfolded Counts");

            graph_unfolded->Draw("P SAME");

            // Add title with averages
            std::string title = Form("DVCS Combined, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                     period_names[period].c_str(),
                                     unfolding_data.xB_avg,
                                     unfolding_data.Q2_avg,
                                     unfolding_data.t_avg);
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5
        }

        // Save the canvas
        std::string filename = output_subdir + "/unfolding_results_DVCS_" + period_names[period] +
                               "_xB_bin_" + std::to_string(xB_bin) + "_combined.pdf";
        canvas->SaveAs(filename.c_str());

        // Clean up
        delete canvas;
    }

    // Similar plotting code for eppi0 data (periods 3-5)
    for (int period = 0; period < n_periods; ++period) {
        int eppi0_period = period + 3;

        // Create output directories
        std::string analysis_type = "eppi0";
        std::string output_subdir = output_dir + "/unfolded/" + analysis_type + "/" + period_names[period];
        gSystem->mkdir(output_subdir.c_str(), true);

        size_t n_bins = combined_unfolding_data.size();
        int n_columns = std::ceil(std::sqrt(n_bins));
        int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

        // Create the canvas
        TCanvas* canvas = new TCanvas(Form("c_unfolding_eppi0_period_%d_%s", eppi0_period, combined_topology.c_str()),
                                      "Unfolding Results eppi0 Combined", canvas_width, canvas_height);
        canvas->Divide(n_columns, n_rows);

        for (size_t group_idx = 0; group_idx < combined_unfolding_data.size(); ++group_idx) {
            UnfoldingData& unfolding_data = combined_unfolding_data[group_idx];
            size_t n_phi_bins = unfolding_data.phi_min.size();

            // Prepare data for plotting
            std::vector<double> phi_centers(n_phi_bins);
            std::vector<double> phi_widths(n_phi_bins);
            std::vector<double> unfolded_counts(n_phi_bins);
            std::vector<double> unfolded_errors(n_phi_bins);

            for (size_t phi_idx = 0; phi_idx < n_phi_bins; ++phi_idx) {
                phi_centers[phi_idx] = (unfolding_data.phi_min[phi_idx] + unfolding_data.phi_max[phi_idx]) / 2.0;
                phi_widths[phi_idx] = (unfolding_data.phi_max[phi_idx] - unfolding_data.phi_min[phi_idx]) / 2.0;

                unfolded_counts[phi_idx] = unfolding_data.unfolded_yields[eppi0_period][phi_idx];
                unfolded_errors[phi_idx] = std::sqrt(unfolding_data.raw_yields[eppi0_period][phi_idx]) /
                                           (unfolding_data.acceptance[eppi0_period][phi_idx] > 0 ? unfolding_data.acceptance[eppi0_period][phi_idx] : 1);
            }

            // Move to the appropriate pad
            canvas->cd(group_idx + 1);

            // Create TGraphErrors
            TGraphErrors* graph_unfolded = new TGraphErrors(n_phi_bins, &phi_centers[0], &unfolded_counts[0], &phi_widths[0], &unfolded_errors[0]);

            // Set styles
            graph_unfolded->SetMarkerColor(kBlack);
            graph_unfolded->SetMarkerStyle(20);
            graph_unfolded->SetLineColor(kBlack);

            // Draw graph
            TPad* pad = (TPad*)gPad;
            pad->SetLeftMargin(0.15);
            pad->SetBottomMargin(0.15);

            // Create a frame histogram for axes
            double phi_min = 0.0;
            double phi_max = 360.0;
            double count_max = *std::max_element(unfolded_counts.begin(), unfolded_counts.end()) * 1.2;
            double count_min = 0.0;

            TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
            frame->GetXaxis()->SetTitle("#phi [deg]");
            frame->GetYaxis()->SetTitle("Unfolded Counts");

            graph_unfolded->Draw("P SAME");

            // Add title with averages
            std::string title = Form("eppi0 Combined, %s, x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f",
                                     period_names[period].c_str(),
                                     unfolding_data.xB_avg,
                                     unfolding_data.Q2_avg,
                                     unfolding_data.t_avg);
            TLatex latex;
            latex.SetNDC();
            latex.SetTextSize(0.04);
            latex.SetTextAlign(22); // Center alignment
            latex.DrawLatex(0.5, 0.95, title.c_str()); // Centered on x=0.5
        }

        // Save the canvas
        std::string filename = output_subdir + "/unfolding_results_eppi0_" + period_names[period] +
                               "_xB_bin_" + std::to_string(xB_bin) + "_combined.pdf";
        canvas->SaveAs(filename.c_str());

        // Clean up
        delete canvas;
    }

    // Return the unfolding data for all topologies
    return topology_unfolding_data;
}