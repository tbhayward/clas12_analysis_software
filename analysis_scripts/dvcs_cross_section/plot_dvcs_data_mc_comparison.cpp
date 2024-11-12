// plot_dvcs_data_mc_comparison.cpp

#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TLatex.h>
#include <TGraphErrors.h>
#include <cmath>
#include <string>
#include <vector>
#include "bin_boundaries.h"
#include "kinematic_cuts.h"
#include "bin_helpers.h"
#include <algorithm>
#include <map>
#include <numeric>
#include <iostream>

// Constant to convert radians to degrees
constexpr double RAD_TO_DEG = 180.0 / M_PI;

// Structure to hold counts and errors for each phi bin
struct BinCounts {
    double phi_center;
    double phi_width;
    double count_data;
    double count_mc_gen;
    double count_mc_rec;
};

// Helper function to check if phi is within bin, accounting for wrap-around
bool phi_in_bin(double phi_deg, double phi_low, double phi_high) {
    if (phi_low <= phi_high) {
        return phi_deg >= phi_low && phi_deg <= phi_high;
    } else {
        // Wrap-around bin
        return phi_deg >= phi_low || phi_deg <= phi_high;
    }
}

void plot_dvcs_data_mc_comparison(const std::string& output_dir, 
                                  const std::string& analysisType, 
                                  const std::string& dataset, 
                                  int xB_bin,
                                  const std::vector<BinBoundary>& bin_boundaries, 
                                  TTreeReader& data_reader, 
                                  TTreeReader& mc_gen_reader, 
                                  TTreeReader& mc_rec_reader) {

    // Precompute the relevant bins for the xB_bin
    std::vector<int> relevant_bins = precompute_relevant_bins(xB_bin, bin_boundaries);

    // Debugging: Print xB_bin and number of bin boundaries
    std::cout << "xB_bin: " << xB_bin << std::endl;
    std::cout << "Number of bin boundaries: " << bin_boundaries.size() << std::endl;

    // Check if relevant_bins is empty
    if (relevant_bins.empty()) {
        std::cerr << "Error: No relevant bins found for xB_bin: " << xB_bin << std::endl;
        return;
    } else {
        std::cout << "Number of relevant bins: " << relevant_bins.size() << std::endl;
        // Optionally, print the relevant bin indices and xB ranges
        for (int idx : relevant_bins) {
            const auto& bin = bin_boundaries[idx];
            std::cout << "Relevant bin idx: " << idx << ", xB: [" << bin.xB_low << ", " << bin.xB_high << "]" << std::endl;
        }
    }

    // Group the bins by (Q2_low, t_low)
    // Since xB is fixed for all relevant bins, we group by Q2_low and t_low
    std::map<std::pair<double, double>, std::vector<int>> bin_groups;

    for (int idx : relevant_bins) {
        const auto& bin = bin_boundaries[idx];
        auto key = std::make_pair(bin.Q2_low, bin.t_low);
        bin_groups[key].push_back(idx);
    }

    // Debugging: Print number of bin groups
    std::cout << "Number of bin groups: " << bin_groups.size() << std::endl;
    if (bin_groups.empty()) {
        std::cerr << "Error: No bin groups created. Exiting function." << std::endl;
        return;
    }

    // Map to store counts for each bin group
    std::map<std::pair<double, double>, std::vector<BinCounts>> bin_counts_map;

    // Initialize bin_counts_map
    for (const auto& group : bin_groups) {
        const auto& key = group.first;
        const auto& idx_list = group.second;

        std::vector<BinCounts> bin_counts_vector;

        for (int idx : idx_list) {
            const auto& bin = bin_boundaries[idx];
            double phi_center = bin.phi_avg; // Use phi_avg instead of (phi_low + phi_high)/2.0
            double phi_width = bin.phi_high - bin.phi_low;

            BinCounts bin_counts;
            bin_counts.phi_center = phi_center;
            bin_counts.phi_width = phi_width;
            bin_counts.count_data = 0.0;
            bin_counts.count_mc_gen = 0.0;
            bin_counts.count_mc_rec = 0.0;

            bin_counts_vector.push_back(bin_counts);
        }

        bin_counts_map[key] = bin_counts_vector;
    }

    // Adjust canvas size
    const int base_canvas_width = 1200;
    const int base_canvas_height = 800;
    int canvas_width = static_cast<int>(1.5 * base_canvas_width);
    int canvas_height = static_cast<int>(1.5 * base_canvas_height);

    // Determine number of columns and rows for the canvas
    int n_bins = bin_groups.size();
    int n_columns = std::ceil(std::sqrt(n_bins));
    int n_rows = std::ceil(static_cast<double>(n_bins) / n_columns);

    TCanvas* canvas = new TCanvas("c1", "Data vs MC", canvas_width, canvas_height);
    canvas->Divide(n_columns, n_rows);

    gStyle->SetOptStat(0);

    int pad_idx = 1;  // Pad index for canvas

    // Restart the readers before declaring TTreeReaderValue objects
    data_reader.Restart();
    mc_gen_reader.Restart();
    mc_rec_reader.Restart();

    // Handle theta_neutral_neutral based on analysis type (dvcs or eppi0)
    std::string theta_variable_name;
    if (analysisType == "dvcs") {
        theta_variable_name = "theta_gamma_gamma";
    } else if (analysisType == "eppi0") {
        theta_variable_name = "theta_pi0_pi0";
    } else {
        std::cerr << "Error: Unknown analysisType '" << analysisType << "'" << std::endl;
        return;
    }

    // Readers for data
    TTreeReaderValue<double> phi_data(data_reader, "phi2");
    TTreeReaderValue<double> xB_data(data_reader, "x");
    TTreeReaderValue<double> Q2_data(data_reader, "Q2");
    TTreeReaderValue<double> t1_data(data_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_data(data_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_data(data_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_data(data_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_data(data_reader, "pTmiss");
    TTreeReaderValue<double> theta_neutral_neutral_data(data_reader, theta_variable_name.c_str());

    // Readers for MC-generated data
    TTreeReaderValue<double> phi_mc_gen(mc_gen_reader, "phi2");
    TTreeReaderValue<double> xB_mc_gen(mc_gen_reader, "x");
    TTreeReaderValue<double> Q2_mc_gen(mc_gen_reader, "Q2");
    TTreeReaderValue<double> t1_mc_gen(mc_gen_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_gen(mc_gen_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_gen(mc_gen_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_gen(mc_gen_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_gen(mc_gen_reader, "pTmiss");
    TTreeReaderValue<double> theta_neutral_neutral_mc_gen(mc_gen_reader, theta_variable_name.c_str());

    // Readers for MC-reconstructed data
    TTreeReaderValue<double> phi_mc_rec(mc_rec_reader, "phi2");
    TTreeReaderValue<double> xB_mc_rec(mc_rec_reader, "x");
    TTreeReaderValue<double> Q2_mc_rec(mc_rec_reader, "Q2");
    TTreeReaderValue<double> t1_mc_rec(mc_rec_reader, "t1");
    TTreeReaderValue<double> open_angle_ep2_mc_rec(mc_rec_reader, "open_angle_ep2");
    TTreeReaderValue<double> Emiss2_mc_rec(mc_rec_reader, "Emiss2");
    TTreeReaderValue<double> Mx2_1_mc_rec(mc_rec_reader, "Mx2_1");
    TTreeReaderValue<double> pTmiss_mc_rec(mc_rec_reader, "pTmiss");
    TTreeReaderValue<double> theta_neutral_neutral_mc_rec(mc_rec_reader, theta_variable_name.c_str());

    // Initialize event counters
    int data_event_count = 0;
    int data_events_assigned = 0;

    // Fill data counts
    while (data_reader.Next()) {
        data_event_count++;
        double phi_deg = *phi_data * RAD_TO_DEG;  
        double xB = *xB_data;
        double Q2 = *Q2_data;
        double t_abs = std::abs(*t1_data);

        // Adjust phi_deg to be within [0, 360)
        phi_deg = std::fmod(phi_deg + 360.0, 360.0);

        // Apply kinematic cuts
        if (!apply_kinematic_cuts(*t1_data, *open_angle_ep2_data, *theta_neutral_neutral_data, *Emiss2_data, *Mx2_1_data, *pTmiss_data)) {
            continue;
        }

        // Check if xB is within the xB_bin range
        const auto& bin_example = bin_boundaries[relevant_bins[0]]; // Representative bin
        if (xB < bin_example.xB_low || xB > bin_example.xB_high) {
            continue;
        }

        // Find the bin group
        bool event_assigned = false;
        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            double Q2_low = key.first;
            double Q2_high = bin_boundaries[idx_list[0]].Q2_high;
            double t_low = key.second;
            double t_high = bin_boundaries[idx_list[0]].t_high;

            if (Q2 >= Q2_low && Q2 <= Q2_high && t_abs >= t_low && t_abs <= t_high) {
                // Now find the phi bin within this group
                std::vector<BinCounts>& bin_counts_vector = bin_counts_map[key];

                for (size_t i = 0; i < idx_list.size(); ++i) {
                    int idx = idx_list[i];
                    const auto& bin = bin_boundaries[idx];

                    if (phi_in_bin(phi_deg, bin.phi_low, bin.phi_high)) {
                        // Increment the count
                        bin_counts_vector[i].count_data += 1.0;
                        data_events_assigned++;
                        event_assigned = true;
                        break;
                    }
                }
                break; // Exit loop after assigning to a bin group
            }
        }
    }
    std::cout << "Total data events processed: " << data_event_count << std::endl;
    std::cout << "Total data events assigned to bins: " << data_events_assigned << std::endl;

    // Similar code for mc_gen and mc_rec
    int mc_gen_event_count = 0;
    int mc_gen_events_assigned = 0;

    // Fill mc_gen counts
    while (mc_gen_reader.Next()) {
        mc_gen_event_count++;
        double phi_deg = *phi_mc_gen * RAD_TO_DEG;  
        double xB = *xB_mc_gen;
        double Q2 = *Q2_mc_gen;
        double t_abs = std::abs(*t1_mc_gen);

        // Adjust phi_deg to be within [0, 360)
        phi_deg = std::fmod(phi_deg + 360.0, 360.0);

        // Apply kinematic cuts
        if (!apply_kinematic_cuts(*t1_mc_gen, *open_angle_ep2_mc_gen, *theta_neutral_neutral_mc_gen, *Emiss2_mc_gen, *Mx2_1_mc_gen, *pTmiss_mc_gen)) {
            continue;
        }

        // Check if xB is within the xB_bin range
        const auto& bin_example = bin_boundaries[relevant_bins[0]]; // Representative bin
        if (xB < bin_example.xB_low || xB > bin_example.xB_high) {
            continue;
        }

        // Find the bin group
        bool event_assigned = false;
        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            double Q2_low = key.first;
            double Q2_high = bin_boundaries[idx_list[0]].Q2_high;
            double t_low = key.second;
            double t_high = bin_boundaries[idx_list[0]].t_high;

            if (Q2 >= Q2_low && Q2 <= Q2_high && t_abs >= t_low && t_abs <= t_high) {
                // Now find the phi bin within this group
                std::vector<BinCounts>& bin_counts_vector = bin_counts_map[key];

                for (size_t i = 0; i < idx_list.size(); ++i) {
                    int idx = idx_list[i];
                    const auto& bin = bin_boundaries[idx];

                    if (phi_in_bin(phi_deg, bin.phi_low, bin.phi_high)) {
                        // Increment the count
                        bin_counts_vector[i].count_mc_gen += 1.0;
                        mc_gen_events_assigned++;
                        event_assigned = true;
                        break;
                    }
                }
                break; // Exit loop after assigning to a bin group
            }
        }
    }
    std::cout << "Total mc_gen events processed: " << mc_gen_event_count << std::endl;
    std::cout << "Total mc_gen events assigned to bins: " << mc_gen_events_assigned << std::endl;

    int mc_rec_event_count = 0;
    int mc_rec_events_assigned = 0;

    // Fill mc_rec counts
    while (mc_rec_reader.Next()) {
        mc_rec_event_count++;
        double phi_deg = *phi_mc_rec * RAD_TO_DEG;  
        double xB = *xB_mc_rec;
        double Q2 = *Q2_mc_rec;
        double t_abs = std::abs(*t1_mc_rec);

        // Adjust phi_deg to be within [0, 360)
        phi_deg = std::fmod(phi_deg + 360.0, 360.0);

        // Apply kinematic cuts
        if (!apply_kinematic_cuts(*t1_mc_rec, *open_angle_ep2_mc_rec, *theta_neutral_neutral_mc_rec, *Emiss2_mc_rec, *Mx2_1_mc_rec, *pTmiss_mc_rec)) {
            continue;
        }

        // Check if xB is within the xB_bin range
        const auto& bin_example = bin_boundaries[relevant_bins[0]]; // Representative bin
        if (xB < bin_example.xB_low || xB > bin_example.xB_high) {
            continue;
        }

        // Find the bin group
        bool event_assigned = false;
        for (const auto& group : bin_groups) {
            const auto& key = group.first;
            const auto& idx_list = group.second;

            double Q2_low = key.first;
            double Q2_high = bin_boundaries[idx_list[0]].Q2_high;
            double t_low = key.second;
            double t_high = bin_boundaries[idx_list[0]].t_high;

            if (Q2 >= Q2_low && Q2 <= Q2_high && t_abs >= t_low && t_abs <= t_high) {
                // Now find the phi bin within this group
                std::vector<BinCounts>& bin_counts_vector = bin_counts_map[key];

                for (size_t i = 0; i < idx_list.size(); ++i) {
                    int idx = idx_list[i];
                    const auto& bin = bin_boundaries[idx];

                    if (phi_in_bin(phi_deg, bin.phi_low, bin.phi_high)) {
                        // Increment the count
                        bin_counts_vector[i].count_mc_rec += 1.0;
                        mc_rec_events_assigned++;
                        event_assigned = true;
                        break;
                    }
                }
                break; // Exit loop after assigning to a bin group
            }
        }
    }
    std::cout << "Total mc_rec events processed: " << mc_rec_event_count << std::endl;
    std::cout << "Total mc_rec events assigned to bins: " << mc_rec_events_assigned << std::endl;

    // Now create TGraphErrors for each bin group
    pad_idx = 1;
    for (const auto& group : bin_groups) {
        const auto& key = group.first;
        const auto& bin_counts_vector = bin_counts_map[key];

        size_t n_points = bin_counts_vector.size();

        if (n_points == 0) {
            std::cerr << "Warning: No phi bins for bin group. Skipping plotting for this bin group." << std::endl;
            ++pad_idx;
            continue;
        }

        // Arrays to hold phi centers and counts
        std::vector<double> phi_centers(n_points);
        std::vector<double> phi_errors(n_points, 0.0); // Set x-errors to zero or phi_width/2 if desired
        std::vector<double> counts_data(n_points);
        std::vector<double> errors_data(n_points);
        std::vector<double> counts_mc_gen(n_points);
        std::vector<double> errors_mc_gen(n_points);
        std::vector<double> counts_mc_rec(n_points);
        std::vector<double> errors_mc_rec(n_points);

        // Fill the arrays
        for (size_t i = 0; i < n_points; ++i) {
            const auto& bin_counts = bin_counts_vector[i];
            phi_centers[i] = bin_counts.phi_center;
            phi_errors[i] = 0.0; // bin_counts.phi_width / 2.0; // Optional

            counts_data[i] = bin_counts.count_data;
            errors_data[i] = std::sqrt(bin_counts.count_data); // Before normalization

            counts_mc_gen[i] = bin_counts.count_mc_gen;
            errors_mc_gen[i] = std::sqrt(bin_counts.count_mc_gen); // Before normalization

            counts_mc_rec[i] = bin_counts.count_mc_rec;
            errors_mc_rec[i] = std::sqrt(bin_counts.count_mc_rec); // Before normalization
        }

        // Sort the data by phi_centers to ensure proper plotting
        std::vector<size_t> indices(n_points);
        std::iota(indices.begin(), indices.end(), 0);
        std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
            return phi_centers[a] < phi_centers[b];
        });

        // Reorder all arrays according to sorted indices
        auto reorder = [&](std::vector<double>& vec) {
            std::vector<double> temp(n_points);
            for (size_t i = 0; i < n_points; ++i) {
                temp[i] = vec[indices[i]];
            }
            vec = std::move(temp);
        };

        reorder(phi_centers);
        reorder(phi_errors);
        reorder(counts_data);
        reorder(errors_data);
        reorder(counts_mc_gen);
        reorder(errors_mc_gen);
        reorder(counts_mc_rec);
        reorder(errors_mc_rec);

        // Total counts for scaling
        double sum_data = std::accumulate(counts_data.begin(), counts_data.end(), 0.0);
        double sum_mc_rec = std::accumulate(counts_mc_rec.begin(), counts_mc_rec.end(), 0.0);
        double sum_mc_gen = std::accumulate(counts_mc_gen.begin(), counts_mc_gen.end(), 0.0);

        // Scale reconstructed MC to have the same total counts as data
        if (sum_mc_rec > 0 && sum_data > 0) {
            double scale_factor = sum_data / sum_mc_rec;
            for (size_t i = 0; i < n_points; ++i) {
                counts_mc_rec[i] *= scale_factor;
                errors_mc_rec[i] *= scale_factor;
            }
        }

        // Now scale generated MC using the ratio of its original integral to the original reconstructed MC integral
        if (sum_mc_gen > 0 && sum_mc_rec > 0 && sum_data > 0) {
            double scale_factor = (sum_mc_gen / sum_mc_rec) * (sum_data / sum_mc_rec);
            for (size_t i = 0; i < n_points; ++i) {
                counts_mc_gen[i] *= scale_factor;
                errors_mc_gen[i] *= scale_factor;
            }
        }

        // Create TGraphErrors
        TGraphErrors* graph_data = new TGraphErrors(n_points, &phi_centers[0], &counts_data[0], &phi_errors[0], &errors_data[0]);
        TGraphErrors* graph_mc_gen = new TGraphErrors(n_points, &phi_centers[0], &counts_mc_gen[0], &phi_errors[0], &errors_mc_gen[0]);
        TGraphErrors* graph_mc_rec = new TGraphErrors(n_points, &phi_centers[0], &counts_mc_rec[0], &phi_errors[0], &errors_mc_rec[0]);

        // Set styles
        graph_data->SetMarkerColor(kBlue);
        graph_data->SetMarkerStyle(20);
        graph_data->SetLineColor(kBlue);

        graph_mc_gen->SetLineColor(kBlack);
        graph_mc_gen->SetLineStyle(2);

        graph_mc_rec->SetMarkerColor(kRed);
        graph_mc_rec->SetMarkerStyle(22);
        graph_mc_rec->SetLineColor(kRed);

        // Draw graphs
        canvas->cd(pad_idx);
        TPad* pad = (TPad*)gPad;
        pad->SetLeftMargin(0.15);
        pad->SetBottomMargin(0.15);
        pad->SetLogy(); // Set y-axis to logarithmic scale

        // Create a frame histogram for axes
        double phi_min = *std::min_element(phi_centers.begin(), phi_centers.end()) - 5.0;
        double phi_max = *std::max_element(phi_centers.begin(), phi_centers.end()) + 5.0;
        double count_max = std::max(*std::max_element(counts_data.begin(), counts_data.end()),
                                    std::max(*std::max_element(counts_mc_rec.begin(), counts_mc_rec.end()),
                                             *std::max_element(counts_mc_gen.begin(), counts_mc_gen.end()))) * 5.0; // Adjusted for log scale
        double count_min = 0.1; // Set minimum count for log scale

        TH1F* frame = pad->DrawFrame(phi_min, count_min, phi_max, count_max);
        frame->GetXaxis()->SetTitle("#phi [deg]");
        frame->GetYaxis()->SetTitle("Counts");

        graph_data->Draw("P SAME");
        graph_mc_rec->Draw("P SAME");
        graph_mc_gen->Draw("L SAME");

        // Add legend
        TLegend* legend = new TLegend(0.575, 0.45, 0.9, 0.75);
        legend->AddEntry(graph_data, "Data", "lep");
        legend->AddEntry(graph_mc_rec, "Reconstructed MC", "lep");
        legend->AddEntry(graph_mc_gen, "Generated MC", "l");
        legend->SetTextSize(0.04);
        legend->Draw();

        // Add title with averages instead of ranges
        const auto& bin = bin_boundaries[bin_groups[key][0]];
        std::string title = Form("%s, %s: x_{B}=%.3f, Q^{2}=%.3f, -t=%.3f", 
                                 analysisType.c_str(), 
                                 dataset.c_str(),  
                                 bin.xB_avg,
                                 bin.Q2_avg,
                                 bin.t_avg);
        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.15, 0.85, title.c_str());

        ++pad_idx;
    }

    // Save canvas
    std::string filename = output_dir + "/phi_data_mc_comparison_" + analysisType + "_" + dataset + "_xB_bin_" + std::to_string(xB_bin) + ".pdf";
    canvas->SaveAs(filename.c_str());

    // Clean up
    delete canvas;
}