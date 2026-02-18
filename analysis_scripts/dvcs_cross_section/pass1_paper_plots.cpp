#include "pass1_paper_plots.h"
#include "model_predictions.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TBox.h>
#include <TPad.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <filesystem>

namespace {

struct Row {
    int bin;
    double xB;
    double Q2;
    double t;
    double phi;
    double val;
    double stat;
    double syst;
};

std::vector<std::string> split(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while(std::getline(ss, item, ',')) {
        out.push_back(item);
    } //endwhile
    return out;
}

} // anonymous namespace


bool make_pass1_phi_panels(const std::string &csv_path,
                           double Ebeam,
                           int helicity)
{
    std::ifstream file(csv_path.c_str());
    if(!file) {
        std::cerr << "Cannot open CSV: " << csv_path << "\n";
        return false;
    }

    std::vector<Row> rows;

    std::string line;
    std::getline(file, line); // skip header

    while(std::getline(file, line)) {
        if(line.empty()) continue;

        auto cols = split(line);
        if(cols.size() < 8) continue;

        Row r;
        r.bin  = std::stoi(cols[0]);
        r.xB   = std::stod(cols[1]);
        r.Q2   = std::stod(cols[2]);
        r.t    = std::stod(cols[3]);
        r.phi  = std::stod(cols[4]);
        r.val  = std::stod(cols[5]);
        r.stat = std::stod(cols[6]);
        r.syst = std::stod(cols[7]);

        rows.push_back(r);
    } //endwhile

    std::vector<int> wanted_bins = {63,65,66,67,68};

    std::map<int, std::vector<Row>> grouped;

    for(const auto &r : rows) {
        for(int b : wanted_bins) {
            if(r.bin == b) {
                grouped[b].push_back(r);
            }
        } //endfor
    } //endfor

    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    std::filesystem::create_directories("output/pass1_paper_plots");

    TCanvas *c = new TCanvas("c_pass1", "pass1", 1500, 900);
    c->Divide(3,2,0.0001,0.0001);

    int pad_index = 1;

    for(int b : wanted_bins) {

        c->cd(pad_index);
        TPad *p = (TPad*)gPad;
        p->SetLeftMargin(0.15);
        p->SetBottomMargin(0.15);
        p->SetRightMargin(0.05);
        p->SetTopMargin(0.08);

        auto &vec = grouped[b];
        if(vec.empty()) continue;

        double tval = vec.front().t;

        TH1F *frame = (TH1F*)p->DrawFrame(0, 1e-3, 360, 1.2);
        frame->GetYaxis()->SetTitleOffset(1.6);
        frame->GetYaxis()->SetTitle(
            "#frac{d#sigma_{ep#rightarrow e'p'#gamma}}{dx_{B} dQ^{2} d|t| d#phi} (nb/GeV^{4})"
        );

        if(pad_index >= 4) {
            frame->GetXaxis()->SetTitle("#phi (deg)");
        }

        p->SetLogy();

        int N = vec.size();

        TGraphErrors *g_stat = new TGraphErrors(N);
        TGraphErrors *g_syst = new TGraphErrors(N);

        for(int i=0;i<N;i++) {
            g_stat->SetPoint(i, vec[i].phi, vec[i].val);
            g_stat->SetPointError(i, 0, vec[i].stat);

            g_syst->SetPoint(i, vec[i].phi, vec[i].val);
            g_syst->SetPointError(i, 0, vec[i].syst);
        } //endfor

        g_syst->SetFillColor(kGray);
        g_syst->SetLineColor(kGray);
        g_syst->Draw("E2 SAME");

        g_stat->SetMarkerStyle(20);
        g_stat->SetMarkerColor(kBlack);
        g_stat->SetLineColor(kBlack);
        g_stat->Draw("PE SAME");

        // -------- Models --------

        const int nphi = 200;
        TGraph *g_bh   = new TGraph(nphi);
        TGraph *g_vgg  = new TGraph(nphi);
        TGraph *g_km15 = new TGraph(nphi);

        for(int i=0;i<nphi;i++) {
            double phi = 360.0 * i / (nphi-1);

            double bh  = vgg_bh_only(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam);
            double vgg = vgg_xs(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam, helicity);
            double km  = km15_xs(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam, helicity);

            g_bh->SetPoint(i, phi, bh);
            g_vgg->SetPoint(i, phi, vgg);
            g_km15->SetPoint(i, phi, km);
        } //endfor

        g_bh->SetLineColor(kRed);
        g_bh->SetLineWidth(2);

        g_vgg->SetLineColor(kOrange+7);
        g_vgg->SetLineWidth(2);

        g_km15->SetLineColor(kCyan+2);
        g_km15->SetLineWidth(2);

        g_bh->Draw("L SAME");
        g_vgg->Draw("L SAME");
        g_km15->Draw("L SAME");

        TLatex tl;
        tl.SetNDC();
        tl.SetTextSize(0.05);
        tl.DrawLatex(0.18,0.88,
            Form("|t| = %.2f GeV^{2}", tval));

        if(pad_index == 6) {
            TLatex info;
            info.SetNDC();
            info.SetTextSize(0.045);
            info.DrawLatex(0.55,0.75,
                "<x_{B}> = 0.17");
            info.DrawLatex(0.55,0.68,
                "<Q^{2}> = 2.24 GeV^{2}");
        }

        pad_index++;

    } //endfor bins

    c->SaveAs("output/pass1_paper_plots/pass1_phi_panels.png");

    delete c;

    std::cout << "Saved: output/pass1_paper_plots/pass1_phi_panels.png\n";

    return true;
}