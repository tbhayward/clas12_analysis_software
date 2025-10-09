#include "uncorrected_cross_section.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1.h>
#include <TAxis.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <nlohmann/json.hpp>

using json = nlohmann::json;
namespace fs = std::filesystem;

// ------------------------------------------------------------
// Helper: Read total integrated luminosity from text file
// ------------------------------------------------------------
// ------------------------------------------------------------
// Helper: Read total integrated luminosity from CSV text file
// (sum of 3rd and 4th columns = pos + neg)
// ------------------------------------------------------------
double read_total_luminosity(const std::string& filepath) {
    std::ifstream in(filepath);
    if (!in.is_open()) {
        std::cerr << "[luminosity] Cannot open " << filepath << "\n";
        return 0.0;
    }

    auto trim = [](std::string s){
        size_t a = s.find_first_not_of(" \t\r\n");
        size_t b = s.find_last_not_of(" \t\r\n");
        if (a == std::string::npos) return std::string();
        return s.substr(a, b - a + 1);
    };

    double sum_pos = 0.0, sum_neg = 0.0;
    std::string line;

    while (std::getline(in, line)) {
        line = trim(line);
        if (line.empty()) continue;
        if (line[0] == '#') continue;

        // CSV: run,total,pos,neg,*,*
        // Example: 5036,169255.790,87120.751,87500.855,0,0
        std::stringstream ss(line);
        std::string tok;
        std::vector<std::string> cols;
        while (std::getline(ss, tok, ',')) cols.push_back(trim(tok));
        if (cols.size() < 4) continue;

        // cols[2] = pos, cols[3] = neg
        try {
            double pos = std::stod(cols[2]);
            double neg = std::stod(cols[3]);
            if (pos > 0) sum_pos += pos;
            if (neg > 0) sum_neg += neg;
        } catch (...) {
            // ignore malformed lines
        }
    }

    in.close();
    std::cout << "[luminosity] " << filepath
              << "  total_pos=" << sum_pos
              << "  total_neg=" << sum_neg << "\n";
    return sum_pos + sum_neg;
}

// ------------------------------------------------------------
// Load bin volume JSON for energy
// ------------------------------------------------------------
json load_json(const std::string& path) {
    std::ifstream f(path);
    if (!f.is_open()) {
        std::cerr << "[json] Failed to open " << path << "\n";
        return json();
    }
    json j;
    f >> j;
    return j;
}

// ------------------------------------------------------------
// Calculate uncorrected cross section
// ------------------------------------------------------------
void compute_uncorrected_cross_sections(
    const std::vector<Binning>& binning_scheme,
    const std::string& bin_volume_json_dir,
    const std::string& unfolded_counts_dir,
    const std::string& luminosity_dir,
    const std::string& output_dir)
{
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetLineWidth(2);
    gStyle->SetFrameLineWidth(2);

    // ----------------------------
    // Load luminosities
    // ----------------------------
    const double L_fa18_inb = read_total_luminosity(luminosity_dir + "/rga_fa18_inb.txt");
    const double L_fa18_out = read_total_luminosity(luminosity_dir + "/rga_fa18_out.txt");
    const double L_sp19_inb = read_total_luminosity(luminosity_dir + "/rga_sp19_inb.txt");

    const double L_10p6 = L_fa18_inb + L_fa18_out;
    const double L_10p2 = L_sp19_inb;

    const std::map<std::string, double> energy_luminosities = {
        {"10.60", L_10p6},
        {"10.2",  L_10p2}
    };

    // ----------------------------
    // Load bin volume JSONs
    // ----------------------------
    const std::vector<std::string> energies = {"10.60", "10.2"};
    std::map<std::string, json> binvolumes;
    for (const auto& e : energies) {
        std::string path = bin_volume_json_dir + "/bin_volume_" + e + ".json";
        binvolumes[e] = load_json(path);
        if (binvolumes[e].empty()) {
            std::cerr << "[binvol] Missing " << path << "\n";
        } else {
            std::cout << "[binvol] Loaded " << path << "\n";
        }
    }

    // ----------------------------
    // Prepare output
    // ----------------------------
    fs::create_directories(output_dir);
    fs::create_directories(output_dir + "/jsons");
    fs::create_directories(output_dir + "/plots");

    // ----------------------------
    // Loop over energies
    // ----------------------------
    for (const auto& e : energies) {
        double L_total = energy_luminosities.at(e);
        json& jvol = binvolumes[e];

        if (jvol.empty() || L_total == 0) {
            std::cerr << "[xsec] Skipping " << e << " due to missing data.\n";
            continue;
        }

        // Load unfolded counts JSON (placeholder paths for now)
        std::string counts_file_pos = unfolded_counts_dir + "/unfolded_counts_pos_" + e + ".json";
        std::string counts_file_neg = unfolded_counts_dir + "/unfolded_counts_neg_" + e + ".json";

        json counts_pos = load_json(counts_file_pos);
        json counts_neg = load_json(counts_file_neg);

        // Construct output json
        json jxsec;
        jxsec["energy"] = e;
        jxsec["luminosity_total"] = L_total;

        // Each helicity handled separately
        for (const std::string& helicity : {"pos", "neg"}) {
            const json& counts = (helicity == "pos" ? counts_pos : counts_neg);

            if (counts.empty()) {
                std::cerr << "[xsec] No unfolded counts for " << e << " helicity=" << helicity << "\n";
                continue;
            }

            json helicity_bins;

            // Loop through bins
            for (auto it = jvol["bins"].begin(); it != jvol["bins"].end(); ++it) {
                std::string bin_key = it.key();
                const auto& vol_entry = it.value();

                double bin_volume = 0.0;
                if (vol_entry.contains("vol")) {
                    double sum_vol = 0.0;
                    for (double v : vol_entry["vol"]) sum_vol += v;
                    bin_volume = sum_vol;
                }

                double counts_bin = 0.0;
                if (counts.contains(bin_key)) counts_bin = counts[bin_key].get<double>();

                double xsec = (bin_volume > 0.0 && L_total > 0.0)
                              ? counts_bin / (bin_volume * L_total)
                              : 0.0;

                helicity_bins[bin_key] = {
                    {"counts", counts_bin},
                    {"bin_volume", bin_volume},
                    {"xsec_uncorr", xsec}
                };
            }

            jxsec["helicity_" + helicity] = helicity_bins;
        }

        // Save JSON
        std::string out_json = output_dir + "/jsons/uncorrected_xsec_" + e + ".json";
        std::ofstream ofs(out_json);
        ofs << std::setw(2) << jxsec;
        ofs.close();
        std::cout << "[xsec] Wrote " << out_json << "\n";

        // --------------------------
        // Plot cross section vs φ
        // --------------------------
        TCanvas* c = new TCanvas(("c_xsec_" + e).c_str(), ("c_xsec_" + e).c_str(), 1000, 600);
        c->Divide(1,2);

        for (int h=0; h<2; ++h) {
            c->cd(h+1);
            gPad->SetGrid(1,1);
            gPad->SetTopMargin(0.08);
            gPad->SetBottomMargin(0.15);
            gPad->SetLeftMargin(0.15);
            gPad->SetRightMargin(0.10);
            gPad->SetLogy();

            TH1* frame = gPad->DrawFrame(0, 1e-4, 360, 1.0);
            frame->GetXaxis()->SetTitle("#phi (deg)");
            frame->GetYaxis()->SetTitle("Uncorrected dσ/d#phi (arb. units)");
            frame->GetYaxis()->CenterTitle();
            frame->GetXaxis()->CenterTitle();

            TGraphErrors* gr = new TGraphErrors();
            int p=0;
            for (int phi_i=0; phi_i<12; ++phi_i) {
                double phi_deg = (phi_i+0.5)*30.0;
                double val = 1e-3 * (1 + 0.5*sin(phi_deg*M_PI/180.0)); // placeholder
                gr->SetPoint(p++, phi_deg, val);
            }
            gr->SetMarkerStyle(20);
            gr->SetMarkerSize(0.9);
            gr->SetLineWidth(2);
            gr->SetLineColor(h==0 ? kRed : kBlue);
            gr->Draw("P SAME");

            TLatex lab;
            lab.SetNDC();
            lab.SetTextSize(0.05);
            lab.DrawLatex(0.15,0.93,(h==0?"Positive helicity":"Negative helicity"));
        }

        std::string out_plot = output_dir + "/plots/uncorrected_xsec_" + e + ".png";
        c->SaveAs(out_plot.c_str());
        delete c;
    }

    std::cout << "[xsec] Finished uncorrected cross-section generation.\n";
}