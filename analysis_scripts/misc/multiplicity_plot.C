#include <iostream>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <cmath>

// Include ROOT headers
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TApplication.h"
#include "TROOT.h"

void multiplicity_plot(const char* filename)
{
    // Open the ROOT file
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file " << filename << std::endl;
        return;
    }

    // Get the TTree named "PhysicsEvents"
    TTree *tree = (TTree*)file->Get("PhysicsEvents");
    if (!tree) {
        std::cerr << "Error: TTree 'PhysicsEvents' not found in file " << filename << std::endl;
        file->Close();
        return;
    }

    // Enable only the necessary branches
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("runnum", 1);
    tree->SetBranchStatus("num_pos", 1);
    tree->SetBranchStatus("num_neg", 1);
    tree->SetBranchStatus("num_neutral", 1);
    // Enable branches for sector definitions
    tree->SetBranchStatus("Q2", 1);
    tree->SetBranchStatus("W", 1);
    tree->SetBranchStatus("y", 1);
    tree->SetBranchStatus("Mx2", 1);
    tree->SetBranchStatus("e_phi", 1);

    // Set branch addresses
    int runnum, num_pos, num_neg, num_neutral;
    double Q2, W, y, Mx2, e_phi;

    tree->SetBranchAddress("runnum", &runnum);
    tree->SetBranchAddress("num_pos", &num_pos);
    tree->SetBranchAddress("num_neg", &num_neg);
    tree->SetBranchAddress("num_neutral", &num_neutral);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W", &W);
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("Mx2", &Mx2);
    tree->SetBranchAddress("e_phi", &e_phi);

    // Define data structures to accumulate sums per run number
    struct Accumulator {
        double sum_pos = 0;
        double sum_pos2 = 0;
        double sum_neg = 0;
        double sum_neg2 = 0;
        double sum_neutral = 0;
        double sum_neutral2 = 0;
        int count = 0;
    };

    std::map<int, Accumulator> run_data;

    // Data structures for per-sector accumulation
    std::map<int, std::map<int, Accumulator>> sector_run_data; // sector -> runnum -> Accumulator

    // Loop over the tree once
    Long64_t nentries = tree->GetEntries();
    for (Long64_t i = 0; i < nentries; i++) {
        tree->GetEntry(i);

        // Accumulate per run number as before
        Accumulator &acc = run_data[runnum];
        acc.sum_pos += num_pos;
        acc.sum_pos2 += num_pos * num_pos;
        acc.sum_neg += num_neg;
        acc.sum_neg2 += num_neg * num_neg;
        acc.sum_neutral += num_neutral;
        acc.sum_neutral2 += num_neutral * num_neutral;
        acc.count += 1;

        // Determine sector
        int sector = 0;
        if (Q2 > 1 && W > 2 && y < 0.80 && Mx2 > 1.8225) {
            if (e_phi < 0.7 || e_phi > 5.9) {
                sector = 1;
            } else if (e_phi > 0.7 && e_phi < 1.8) {
                sector = 2;
            } else if (e_phi > 1.8 && e_phi < 2.8) {
                sector = 3;
            } else if (e_phi > 2.8 && e_phi < 3.8) {
                sector = 4;
            } else if (e_phi > 3.8 && e_phi < 4.8) {
                sector = 5;
            } else if (e_phi > 4.8 && e_phi < 5.85) {
                sector = 6;
            }
        }

        // If sector is non-zero, accumulate per sector per run number
        if (sector != 0) {
            Accumulator &acc_sector = sector_run_data[sector][runnum];
            acc_sector.sum_pos += num_pos;
            acc_sector.sum_pos2 += num_pos * num_pos;
            acc_sector.sum_neg += num_neg;
            acc_sector.sum_neg2 += num_neg * num_neg;
            acc_sector.sum_neutral += num_neutral;
            acc_sector.sum_neutral2 += num_neutral * num_neutral;
            acc_sector.count += 1;
        }
    }

    // Prepare vectors for plotting (for the overall data)
    std::vector<double> runnum_vals, mean_pos, mean_neg, mean_neutral;
    std::vector<double> err_pos, err_neg, err_neutral;

    // Calculate means and errors for each run number
    for (const auto& kv : run_data) {
        int curr_runnum = kv.first;
        const Accumulator &acc = kv.second;
        int N = acc.count;

        double mean_p = acc.sum_pos / N;
        double mean_n = acc.sum_neg / N;
        double mean_neu = acc.sum_neutral / N;

        // Variance calculation
        double var_p = (acc.sum_pos2 / N) - (mean_p * mean_p);
        double var_n = (acc.sum_neg2 / N) - (mean_n * mean_n);
        double var_neu = (acc.sum_neutral2 / N) - (mean_neu * mean_neu);

        // Statistical uncertainty (standard error)
        double err_p = sqrt(var_p / N);
        double err_n = sqrt(var_n / N);
        double err_neu = sqrt(var_neu / N);

        runnum_vals.push_back(curr_runnum);
        mean_pos.push_back(mean_p);
        mean_neg.push_back(mean_n - 1);
        mean_neutral.push_back(mean_neu);

        err_pos.push_back(err_p);
        err_neg.push_back(err_n);
        err_neutral.push_back(err_neu);
    }

    // Sort data by run number
    std::vector<size_t> indices(runnum_vals.size());
    std::iota(indices.begin(), indices.end(), 0);

    std::sort(indices.begin(), indices.end(),
        [&runnum_vals](size_t i1, size_t i2) { return runnum_vals[i1] < runnum_vals[i2]; });

    // Create sorted vectors
    std::vector<double> runnum_vals_sorted, mean_pos_sorted, mean_neg_sorted, mean_neutral_sorted;
    std::vector<double> err_pos_sorted, err_neg_sorted, err_neutral_sorted;

    for (size_t idx : indices) {
        runnum_vals_sorted.push_back(runnum_vals[idx]);
        mean_pos_sorted.push_back(mean_pos[idx]);
        mean_neg_sorted.push_back(mean_neg[idx]);
        mean_neutral_sorted.push_back(mean_neutral[idx]);
        err_pos_sorted.push_back(err_pos[idx]);
        err_neg_sorted.push_back(err_neg[idx]);
        err_neutral_sorted.push_back(err_neutral[idx]);
    }

    // Create TGraphErrors for plotting
    TGraphErrors *gr_pos = new TGraphErrors(runnum_vals_sorted.size(), &runnum_vals_sorted[0], &mean_pos_sorted[0], nullptr, &err_pos_sorted[0]);
    TGraphErrors *gr_neg = new TGraphErrors(runnum_vals_sorted.size(), &runnum_vals_sorted[0], &mean_neg_sorted[0], nullptr, &err_neg_sorted[0]);
    TGraphErrors *gr_neutral = new TGraphErrors(runnum_vals_sorted.size(), &runnum_vals_sorted[0], &mean_neutral_sorted[0], nullptr, &err_neutral_sorted[0]);

    // Customize graph appearance
    gr_pos->SetMarkerStyle(20);
    gr_pos->SetMarkerColor(kRed);
    gr_pos->SetLineColor(kRed);

    gr_neg->SetMarkerStyle(21);
    gr_neg->SetMarkerColor(kBlue);
    gr_neg->SetLineColor(kBlue);

    gr_neutral->SetMarkerStyle(22);
    gr_neutral->SetMarkerColor(kGreen+2);
    gr_neutral->SetLineColor(kGreen+2);

    // Create and customize canvas
    TCanvas *c1 = new TCanvas("c1", "Multiplicity vs Run Number", 800, 600);
    c1->SetGrid();

    // Draw graphs
    gr_pos->SetTitle("");
    gr_pos->GetXaxis()->SetTitle("runnum");
    gr_pos->GetYaxis()->SetTitle("Multiplicity");
    gr_pos->GetYaxis()->SetRangeUser(0.0, 1.3);  // Set y-axis range from 0 to 1.3
    gr_pos->Draw("AP");
    gr_neg->Draw("P SAME");
    gr_neutral->Draw("P SAME");

    // Add legend
    TLegend *leg = new TLegend(0.7, 0.75, 0.9, 0.9); // Position in top right
    leg->SetBorderSize(1);           // Draw border around legend
    leg->SetFillColor(kWhite);       // Solid white background
    leg->SetFillStyle(1001);         // Solid fill
    leg->AddEntry(gr_pos, "positives", "P");
    leg->AddEntry(gr_neg, "negatives", "P");
    leg->AddEntry(gr_neutral, "neutrals", "P");
    leg->Draw();

    // Save canvas to file
    c1->SaveAs("nh3_multiplicity.pdf");

    // Now process per-sector data

    // Create new canvas for sector plots
    TCanvas *c2 = new TCanvas("c2", "Multiplicity vs Run Number per Sector", 1200, 800);
    c2->Divide(3,2); // 2 rows, 3 columns

    for (int sector = 1; sector <=6; sector++) {
        c2->cd(sector);

        // Get the data for this sector
        const auto &run_data_sector = sector_run_data[sector];

        // Prepare vectors
        std::vector<double> runnum_vals_s, mean_pos_s, mean_neg_s, mean_neutral_s;
        std::vector<double> err_pos_s, err_neg_s, err_neutral_s;

        // Calculate means and errors for each run number
        for (const auto& kv : run_data_sector) {
            int curr_runnum = kv.first;
            const Accumulator &acc = kv.second;
            int N = acc.count;

            double mean_p = acc.sum_pos / N;
            double mean_n = acc.sum_neg / N;
            double mean_neu = acc.sum_neutral / N;

            // Variance calculation
            double var_p = (acc.sum_pos2 / N) - (mean_p * mean_p);
            double var_n = (acc.sum_neg2 / N) - (mean_n * mean_n);
            double var_neu = (acc.sum_neutral2 / N) - (mean_neu * mean_neu);

            // Statistical uncertainty (standard error)
            double err_p = sqrt(var_p / N);
            double err_n = sqrt(var_n / N);
            double err_neu = sqrt(var_neu / N);

            runnum_vals_s.push_back(curr_runnum);
            mean_pos_s.push_back(mean_p);
            mean_neg_s.push_back(mean_n - 1);
            mean_neutral_s.push_back(mean_neu);

            err_pos_s.push_back(err_p);
            err_neg_s.push_back(err_n);
            err_neutral_s.push_back(err_neu);
        }

        // Sort data by run number
        std::vector<size_t> indices_s(runnum_vals_s.size());
        std::iota(indices_s.begin(), indices_s.end(), 0);

        std::sort(indices_s.begin(), indices_s.end(),
            [&runnum_vals_s](size_t i1, size_t i2) { return runnum_vals_s[i1] < runnum_vals_s[i2]; });

        // Create sorted vectors
        std::vector<double> runnum_vals_sorted_s, mean_pos_sorted_s, mean_neg_sorted_s, mean_neutral_sorted_s;
        std::vector<double> err_pos_sorted_s, err_neg_sorted_s, err_neutral_sorted_s;

        for (size_t idx : indices_s) {
            runnum_vals_sorted_s.push_back(runnum_vals_s[idx]);
            mean_pos_sorted_s.push_back(mean_pos_s[idx]);
            mean_neg_sorted_s.push_back(mean_neg_s[idx]);
            mean_neutral_sorted_s.push_back(mean_neutral_s[idx]);
            err_pos_sorted_s.push_back(err_pos_s[idx]);
            err_neg_sorted_s.push_back(err_neg_s[idx]);
            err_neutral_sorted_s.push_back(err_neutral_s[idx]);
        }

        // Create TGraphErrors for plotting
        TGraphErrors *gr_pos_s = new TGraphErrors(runnum_vals_sorted_s.size(), &runnum_vals_sorted_s[0], &mean_pos_sorted_s[0], nullptr, &err_pos_sorted_s[0]);
        TGraphErrors *gr_neg_s = new TGraphErrors(runnum_vals_sorted_s.size(), &runnum_vals_sorted_s[0], &mean_neg_sorted_s[0], nullptr, &err_neg_sorted_s[0]);
        TGraphErrors *gr_neutral_s = new TGraphErrors(runnum_vals_sorted_s.size(), &runnum_vals_sorted_s[0], &mean_neutral_sorted_s[0], nullptr, &err_neutral_sorted_s[0]);

        // Customize graph appearance
        gr_pos_s->SetMarkerStyle(20);
        gr_pos_s->SetMarkerColor(kRed);
        gr_pos_s->SetLineColor(kRed);

        gr_neg_s->SetMarkerStyle(21);
        gr_neg_s->SetMarkerColor(kBlue);
        gr_neg_s->SetLineColor(kBlue);

        gr_neutral_s->SetMarkerStyle(22);
        gr_neutral_s->SetMarkerColor(kGreen+2);
        gr_neutral_s->SetLineColor(kGreen+2);

        // Draw graphs
        gr_pos_s->SetTitle(Form("Sector %d", sector));
        gr_pos_s->GetXaxis()->SetTitle("runnum");
        gr_pos_s->GetYaxis()->SetTitle("Multiplicity");
        gr_pos_s->GetYaxis()->SetRangeUser(0.0, 1.3);  // Adjust y-axis range as needed
        gr_pos_s->Draw("AP");
        gr_neg_s->Draw("P SAME");
        gr_neutral_s->Draw("P SAME");

        // Add legend
        TLegend *leg_s = new TLegend(0.7, 0.75, 0.9, 0.9); // Position in top right
        leg_s->SetBorderSize(1);           // Draw border around legend
        leg_s->SetFillColor(kWhite);       // Solid white background
        leg_s->SetFillStyle(1001);         // Solid fill
        leg_s->AddEntry(gr_pos_s, "positives", "P");
        leg_s->AddEntry(gr_neg_s, "negatives", "P");
        leg_s->AddEntry(gr_neutral_s, "neutrals", "P");
        leg_s->Draw();
    }

    // Save canvas to file
    c2->SaveAs("nh3_multiplicity_per_sector.pdf");

    // Clean up
    delete gr_pos;
    delete gr_neg;
    delete gr_neutral;
    delete c1;
    delete c2;
    file->Close();
    delete file;
}

// Main function to run the script from command line
int main(int argc, char** argv) {
    // Initialize ROOT application to handle graphics
    TApplication app("app", &argc, argv);

    if (argc != 2) {
        std::cerr << "Usage: ./multiplicity_plot <root_file>" << std::endl;
        return 1;
    }

    multiplicity_plot(argv[1]);

    // If you want to display the canvases, uncomment the next line
    // app.Run();

    return 0;
}