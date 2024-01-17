#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TObjArray.h>
#include <TObjString.h>

void vertex_study() {
    // Define run periods and channels
    const char* run_periods[] = {"RGA_Fa18_Inb", "RGA_Fa18_Out", "RGA_Sp19_Inb", "RGB_Sp19_Inb", "RGB_Fa19_Out", "RGB_Sp20_Inb"};
    const char* channels[] = {"eX", "epi+X", "epi-X", "epX", "ek-X"};

    // Create a canvas with six subplots
    TCanvas* canvas = new TCanvas("canvas", "Vertex Study", 1200, 800);
    canvas->Divide(3, 2);

    // Loop over each run period
    for (int i = 0; i < 6; i++) {
        canvas->cd(i + 1);
        TH1F* hist_eX = new TH1F("hist_eX", run_periods[i], 100, -10, 15);
        TH1F* hist_epiX = new TH1F("hist_epiX", run_periods[i], 100, -10, 15);
        TH1F* hist_ekX = new TH1F("hist_ekX", run_periods[i], 100, -10, 15);

        // Loop over channels for each run period
        for (int j = 0; j < 5; j++) {
            TObjArray* tokens = TString(run_periods[i]).Tokenize("_");
            TString group = ((TObjString*)tokens->At(0))->GetString();
            TString season_bending = ((TObjString*)tokens->At(1))->GetString() + ((TObjString*)tokens->At(2))->GetString();
            delete tokens; // Free memory

            // Construct file path
            TString file_name = Form("/volatile/clas12/thayward/vertex_studies/%s/%s/%s_%s.root", 
                                     group.ToLower().Data(), season_bending.ToLower().Data(), TString(run_periods[i]).ToLower().Data(), channels[j]);

            TFile* file = new TFile(file_name);
            TTree* tree = (TTree*)file->Get("PhysicsEvents");

            TString branch_name = (j == 0) ? "vz_e" : "vz_p";
            float vz;
            tree->SetBranchAddress(branch_name, &vz);

            // Fill histograms
            for (int k = 0; k < tree->GetEntries(); k++) {
                tree->GetEntry(k);
                if (j == 0) hist_eX->Fill(vz);
                else if (j == 2) hist_epiX->Fill(vz);
                else if (j == 4) hist_ekX->Fill(vz);
            }

            file->Close();
        }

        // Draw histograms
        hist_eX->SetLineColor(kBlack);
        hist_eX->Draw();
        hist_epiX->SetLineColor(kBlue);
        hist_epiX->Draw("same");
        hist_ekX->SetLineColor(kRed);
        hist_ekX->Draw("same");

        // Add legend
        if (i == 0) {
            TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);
            leg->AddEntry(hist_eX, "e^{-}", "l");
            leg->AddEntry(hist_epiX, "#pi^{-}", "l");
            leg->AddEntry(hist_ekX, "k^{-}", "l");
            leg->Draw();
        }

        // Customize plot
        gPad->SetLeftMargin(0.15);
        hist_eX->GetYaxis()->SetTitle("Counts");
        hist_eX->GetYaxis()->CenterTitle();
        hist_eX->GetYaxis()->SetTitleSize(0.05);
        hist_eX->GetXaxis()->SetTitle("v_{z}");
        hist_eX->GetXaxis()->CenterTitle();
        hist_eX->GetXaxis()->SetTitleSize(0.05);
    }

    // Save canvas
    canvas->SaveAs("output/negative_vz.png");
}
