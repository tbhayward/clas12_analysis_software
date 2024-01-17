#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>

void vertex_study() {
    const char* run_periods[] = {"Fa18_Inb", "Fa18_Out", "Sp19_Inb", "Sp19_Inb", "Fa19_Out", "Sp20_Inb"};
    const char* neg_channels[] = {"eX", "epi-X", "ek-X"};  // Only negative tracks
    const char* pos_channels[] = {"epi+X", "epX"};  // Only positive tracks
    const char* colors[] = {"black", "blue", "red"};

    TCanvas* c1 = new TCanvas("c1", "Vertex Position Study", 2000, 1200);
    c1->Divide(3, 2);

    for (int i = 0; i < 6; i++) {
        c1->cd(i+1);
        TLegend* leg = new TLegend(0.1, 0.7, 0.3, 0.9);
        for (int j = 0; j < 3; j++) {
            TString file_path = Form("/volatile/clas12/thayward/vertex_studies/rg%c/%s/rg%c_%s_%s.root", 
                                     (i < 3 ? 'a' : 'b'), run_periods[i], (i < 3 ? 'a' : 'b'), run_periods[i], neg_channels[j]);
            TFile* file = new TFile(file_path);
            if (!file || file->IsZombie()) {
                std::cerr << "Error opening file: " << file_path << std::endl;
                continue;
            }

            TTree* tree = (TTree*)file->Get("PhysicsEvents");
            if (!tree) {
                std::cerr << "Tree PhysicsEvents not found in file: " << file_path << std::endl;
                file->Close();
                continue;
            }
            TTree* tree = (TTree*)file->Get("PhysicsEvents");
            TString hist_name = Form("hist_%d_%d", i, j);
            TH1F* hist = new TH1F(hist_name, run_periods[i], 100, -10, 15);
            TString var_name = (j == 0 ? "vz_e" : "vz_p");
            tree->Draw(Form("%s>>%s", var_name.Data(), hist_name.Data()));
            hist->SetLineColor(j + 1);  // Color codes: 1-Black, 2-Red, 3-Green, 4-Blue, etc.
            if (j == 0) hist->Draw();
            else hist->Draw("SAME");
            leg->AddEntry(hist, neg_channels[j], "l");
        }
        leg->Draw();
        gPad->Modified();
    }

    c1->SaveAs("output/negative_vz.png");
}
