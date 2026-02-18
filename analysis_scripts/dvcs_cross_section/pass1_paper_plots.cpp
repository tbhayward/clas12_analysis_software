#include "pass1_paper_plots.h"
#include "model_predictions.h"

#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TPad.h>
#include <TH1F.h>
#include <TROOT.h>

#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cmath>
#include <iostream>
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

} // namespace


bool make_pass1_phi_panels(const std::string &csv_path,
                           double Ebeam,
                           int helicity_input)
{
    Helicity helicity;
    if(helicity_input == 1) helicity = Helicity::Plus;
    else if(helicity_input == -1) helicity = Helicity::Minus;
    else helicity = Helicity::Unpol;

    std::ifstream file(csv_path.c_str());
    if(!file) {
        std::cerr << "Cannot open CSV\n";
        return false;
    }

    std::vector<Row> rows;
    std::string line;
    std::getline(file, line);

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
    }

    std::vector<int> wanted_bins = {63,65,66,67,68};
    std::map<int, std::vector<Row>> grouped;

    for(const auto &r : rows) {
        for(int b : wanted_bins) {
            if(r.bin == b) grouped[b].push_back(r);
        }
    }

    gStyle->SetOptStat(0);
    gROOT->SetBatch(kTRUE);

    std::filesystem::create_directories("output/pass1_paper_plots");

    TCanvas *c = new TCanvas("c_pass1","",1500,900);
    c->Divide(3,2,0.0,0.0);

    int pad_index = 1;

    for(int row=0; row<2; row++) {
        for(int col=0; col<3; col++) {

            if(pad_index > 5) break;

            int bin = wanted_bins[pad_index-1];

            c->cd(pad_index);
            TPad *p = (TPad*)gPad;

            // Tight margins
            p->SetLeftMargin(col==0 ? 0.18 : 0.05);
            p->SetRightMargin(0.02);
            p->SetBottomMargin(row==1 ? 0.18 : 0.05);
            p->SetTopMargin(0.05);

            auto &vec = grouped[bin];
            if(vec.empty()) continue;

            TH1F *frame = new TH1F("frame","",100,0,360);
            frame->SetMinimum(1e-3);
            frame->SetMaximum(1.0);
            frame->Draw();

            p->SetLogy();

            frame->GetXaxis()->SetTitleSize(0.07);
            frame->GetYaxis()->SetTitleSize(0.07);
            frame->GetXaxis()->CenterTitle();
            frame->GetYaxis()->CenterTitle();

            if(row==1)
                frame->GetXaxis()->SetTitle("#phi (^{#circ})");

            if(col==0)
                frame->GetYaxis()->SetTitle(
                    "#frac{d#sigma_{ep#rightarrow e'p'#gamma}}{dx_{B} dQ^{2} d|t| d#phi} (nb/GeV^{4})"
                );

            int N = vec.size();
            TGraphErrors *g_stat = new TGraphErrors(N);
            TGraphErrors *g_syst = new TGraphErrors(N);

            for(int i=0;i<N;i++) {
                g_stat->SetPoint(i, vec[i].phi, vec[i].val);
                g_stat->SetPointError(i, 0, vec[i].stat);

                g_syst->SetPoint(i, vec[i].phi, vec[i].val);
                g_syst->SetPointError(i, 0, vec[i].syst);
            }

            g_syst->SetFillColor(kGray);
            g_syst->SetLineColor(kGray);
            g_syst->Draw("E2 SAME");

            g_stat->SetMarkerStyle(20);
            g_stat->SetMarkerColor(kBlack);
            g_stat->SetLineColor(kBlack);
            g_stat->Draw("PE SAME");

            const int nphi=8;
            TGraph *g_bh   = new TGraph(nphi);
            TGraph *g_vgg  = new TGraph(nphi);
            TGraph *g_km15 = new TGraph(nphi);

            for(int i=0;i<nphi;i++) {
                double phi = 360.0*i/(nphi-1);
                double bh  = vgg_bh_only(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam);
                double vgg = vgg_xs(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam, helicity);
                double km  = km15_xs(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam, helicity);

                g_bh->SetPoint(i,phi,bh);
                g_vgg->SetPoint(i,phi,vgg);
                g_km15->SetPoint(i,phi,km);
            }

            g_bh->SetLineColor(kRed);
            g_vgg->SetLineColor(kOrange+7);
            g_km15->SetLineColor(kCyan+2);

            g_bh->SetLineWidth(2);
            g_vgg->SetLineWidth(2);
            g_km15->SetLineWidth(2);

            g_bh->Draw("L SAME");
            g_vgg->Draw("L SAME");
            g_km15->Draw("L SAME");

            TLatex tl;
            tl.SetNDC();
            tl.SetTextAlign(22);
            tl.SetTextSize(0.06);
            tl.DrawLatex(0.5,0.93,
                Form("|t| = %.2f GeV^{2}", vec[0].t));

            if(pad_index==5) {
                TLatex info;
                info.SetNDC();
                info.SetTextSize(0.055);

                info.DrawLatex(0.55,0.75,"<x_{B}> = 0.17");
                info.DrawLatex(0.55,0.68,"<Q^{2}> = 2.24 GeV^{2}");

                TLegend *leg = new TLegend(0.55,0.45,0.95,0.65);
                leg->SetBorderSize(0);
                leg->SetTextSize(0.05);

                leg->AddEntry(g_stat,"Data","p");
                leg->AddEntry(g_syst,"Syst. Unc.","f");
                leg->AddEntry(g_bh,"BH","l");
                leg->AddEntry(g_vgg,"VGG","l");
                leg->AddEntry(g_km15,"KM15","l");

                leg->Draw();
            }

            pad_index++;
        }
    }

    c->SaveAs("output/pass1_paper_plots/pass1_phi_panels.png");
    delete c;

    return true;
}