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
#include <TF1.h>
#include <TMath.h>

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
    }
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

    // Updated bins: 64 instead of 63
    std::vector<int> wanted_bins = {64,65,66,67,68};
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
    c->Divide(3,2,0.0,0.0);  // zero spacing

    int pad_index = 1;

    for(int row=0; row<2; row++) {
        for(int col=0; col<3; col++) {

            c->cd(pad_index);
            TPad *p = (TPad*)gPad;

            // Tight margins so pads touch
            p->SetLeftMargin(col==0 ? 0.17 : 0.01);
            p->SetRightMargin(0.01);
            p->SetTopMargin(0.01);
            p->SetBottomMargin(row==1 ? 0.16 : 0.01);

            if(pad_index <= 5) {

                int bin = wanted_bins[pad_index-1];
                auto &vec = grouped[bin];
                if(vec.empty()) {
                    pad_index++;
                    continue;
                }

                TH1F *frame = new TH1F(Form("frame_%d",bin),"",100,0,360);
                frame->SetMinimum(1e-3);
                frame->SetMaximum(1.0);
                frame->Draw();

                p->SetLogy();

                frame->GetXaxis()->SetTitleSize(0.075);
                frame->GetYaxis()->SetTitleSize(0.075);
                frame->GetXaxis()->CenterTitle();
                frame->GetYaxis()->CenterTitle();

                // X-axis labels only on bottom row of each column
                if(row==1 || pad_index==3)
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
                g_stat->Draw("PE SAME");

                // Fit
                TF1 *fit = new TF1(Form("fit_%d",bin),
                    "[0]*(1 + [1]*cos(x*TMath::DegToRad()))",
                    0.0, 360.0);

                fit->SetParameters(vec[0].val, 0.1);
                fit->SetLineColor(kBlack);
                fit->SetLineStyle(2);
                fit->SetLineWidth(2);

                g_stat->Fit(fit,"Q");
                fit->Draw("SAME");

                // Models
                const int nphi=20;
                TGraph *g_bh   = new TGraph(nphi);
                TGraph *g_vgg  = new TGraph(nphi);
                TGraph *g_km15 = new TGraph(nphi);

                for(int i=0;i<nphi;i++) {
                    double phi = 360.0*i/(nphi-1);
                    g_bh->SetPoint(i,phi,
                        vgg_bh_only(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam));
                    g_vgg->SetPoint(i,phi,
                        vgg_xs(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam, helicity));
                    g_km15->SetPoint(i,phi,
                        km15_xs(vec[0].xB, vec[0].Q2, vec[0].t, phi, Ebeam, helicity));
                }

                g_bh->SetLineColor(kRed);
                g_vgg->SetLineColor(kOrange+7);
                g_km15->SetLineColor(kCyan+2);

                g_bh->Draw("L SAME");
                g_vgg->Draw("L SAME");
                g_km15->Draw("L SAME");

                TLatex tl;
                tl.SetNDC();
                tl.SetTextAlign(22);
                tl.SetTextSize(0.06);

                // Shifted slightly down
                tl.DrawLatex(0.5,0.88,
                    Form("|t| = %.2f GeV^{2}", vec[0].t));
            }

            // Sixth panel — annotation block
            if(pad_index == 6) {

                TLatex text;
                text.SetNDC();
                text.SetTextSize(0.07);

                text.DrawLatex(0.15,0.80,"<x_{B}> = 0.17");
                text.DrawLatex(0.15,0.70,"<Q^{2}> = 2.24 GeV^{2}");

                TLegend *leg = new TLegend(0.15,0.40,0.85,0.65);
                leg->SetBorderSize(0);
                leg->SetTextSize(0.06);

                leg->AddEntry((TObject*)0,"Data","p");
                leg->AddEntry((TObject*)0,"Syst. Unc.","f");
                leg->AddEntry((TObject*)0,"Fit (c_{0}(1+c_{1}cos#phi))","l");
                leg->AddEntry((TObject*)0,"BH","l");
                leg->AddEntry((TObject*)0,"VGG","l");
                leg->AddEntry((TObject*)0,"KM15","l");

                leg->Draw();
            }

            pad_index++;
        }
    }

    c->SaveAs("output/pass1_paper_plots/pass1_phi_panels.png");
    delete c;

    return true;
}